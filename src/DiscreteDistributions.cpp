/**
 * @file DiscreteDistributions.cpp
 * @brief Discrete power-law distribution implementation.
 *
 * Adapted from the alrobles/powerlaw repository.
 * Reference: Clauset, Shalizi & Newman (2009) https://arxiv.org/abs/0706.1062
 *
 * Performance improvements over the original:
 * - Alpha MLE uses Brent's method (O(log(1/tol)) zeta calls) instead of grid
 *   search (O(range/precision) calls), giving ~10x faster alpha estimation.
 * - Empirical CDF is built in a single O(n + range) linear pass rather than
 *   O(range * log n) repeated binary searches.
 * - Log-sum of data values is computed once per xMin candidate instead of
 *   being recomputed at every alpha grid point.
 */
#include "DiscreteDistributions.h"
#include "Zeta.h"
#include "VectorUtilities.h"
#include <algorithm>
#include <cmath>
#include <limits>
using namespace std;

/* =========================================================================
 * Brent's method for 1-D function maximisation over a closed interval.
 *
 * Combines the golden-section step and parabolic interpolation following
 * Brent (1973) "Algorithms for Minimization without Derivatives", ch. 5.
 * Convergence is super-linear and does not require derivatives.
 *
 * Template parameter F: callable double(double).
 * @param f        Function to maximise.
 * @param lo       Left bracket.
 * @param hi       Right bracket.
 * @param tol      Absolute tolerance on the returned maximiser.
 * @param maxIter  Maximum number of function evaluations.
 * @return         Approximate maximiser x* in [lo, hi].
 * ======================================================================= */
template<typename F>
static double BrentMaximize(F f,
                             double lo, double hi,
                             double tol     = 1e-8,
                             int    maxIter = 200)
{
    const double phi = (3.0 - std::sqrt(5.0)) / 2.0; /* (3-sqrt5)/2 ≈ 0.382 */

    double a = lo, b = hi;
    double x = a + phi * (b - a);
    double w = x, v = x;
    double fx = f(x), fw = fx, fv = fx;
    double d = 0.0, e = 0.0;

    for (int iter = 0; iter < maxIter; ++iter)
    {
        const double xm   = 0.5 * (a + b);
        const double tol1 = tol * std::abs(x) + 1e-10;
        const double tol2 = 2.0 * tol1;

        /* Convergence test */
        if (std::abs(x - xm) <= tol2 - 0.5 * (b - a))
            return x;

        /* Attempt parabolic interpolation */
        double p = 0.0, q = 0.0, r = 0.0;
        if (std::abs(e) > tol1)
        {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) p = -p; else q = -q;
            r = e;
            e = d;
        }

        if (std::abs(p) < std::abs(0.5 * q * r) &&
            p > q * (a - x) && p < q * (b - x))
        {
            d = p / q; /* Parabolic step */
        }
        else
        {
            e = (x < xm) ? b - x : a - x;
            d = phi * e; /* Golden-section step */
        }

        const double u  = x + (std::abs(d) >= tol1
                               ? d
                               : (d > 0 ? tol1 : -tol1));
        const double fu = f(u);

        /* Update bracket – maximising, so inequalities are flipped vs minimise */
        if (fu >= fx)
        {
            if (u < x) b = x; else a = x;
            v = w; w = x; x = u;
            fv = fw; fw = fx; fx = fu;
        }
        else
        {
            if (u < x) a = u; else b = u;
            if (fu >= fw || w == x)      { v = w; w = u; fv = fw; fw = fu; }
            else if (fu >= fv || v == x || v == w) { v = u; fv = fu; }
        }
    }
    return x;
}

/* =========================================================================
 * DiscreteEmpiricalDistribution
 * ======================================================================= */

DiscreteEmpiricalDistribution::DiscreteEmpiricalDistribution(
    const vector<int>& sampleData, int xMin, int xMax)
{
    vector<int> sortedTailSample = sampleData;
    VectorUtilities::RemoveLower(sortedTailSample, xMin);
    VectorUtilities::RemoveGreater(sortedTailSample, xMax);
    VectorUtilities::Sort(sortedTailSample);

    _xMin = xMin;
    _xMax = xMax;
    PrecalculateCDF(sortedTailSample);
}

void DiscreteEmpiricalDistribution::PrecalculateCDF(const vector<int>& sortedTailSample)
{
    if (_xMax > _xMin)
    {
        /* Build the survival function P(X >= x) for x = _xMin, ..., _xMax.
         *
         * Original approach: O((xMax-xMin) * log n) – repeated binary searches.
         * Optimised approach: O(n + xMax-xMin) – one forward scan through the
         * already-sorted sample, advancing a single index as x increases.
         */
        const int    n     = static_cast<int>(sortedTailSample.size());
        const int    range = _xMax - _xMin + 1;
        _cdf.resize(static_cast<size_t>(range));

        int dataIdx = 0; /* first index with sortedTailSample[dataIdx] >= x */
        for (int i = 0; i < range; ++i)
        {
            const int x = _xMin + i;
            /* Advance past elements strictly less than x */
            while (dataIdx < n && sortedTailSample[dataIdx] < x)
                ++dataIdx;
            /* P(X >= x) = (count of elements >= x) / n */
            _cdf[static_cast<size_t>(i)] =
                static_cast<double>(n - dataIdx) / static_cast<double>(n);
        }
    }
}

double DiscreteEmpiricalDistribution::GetCDF(int x) const
{
    if (x >= _xMin && x <= _xMax)
        return _cdf[x - _xMin];
    else if (x < _xMin)
        return 1.0;
    else
        return 0.0;
}

/* =========================================================================
 * DiscretePowerLawDistribution
 * ======================================================================= */

DiscretePowerLawDistribution::DiscretePowerLawDistribution(
    const DiscretePowerLawDistribution& other)
{
    _alpha            = other._alpha;
    _xMin             = other._xMin;
    _xMax             = other._xMax;
    _state            = other._state;
    _sampleSize       = other._sampleSize;
    _alphaPrecision   = other._alphaPrecision;
    _ksStatistic      = other._ksStatistic;
    _distributionType = other._distributionType;
    _cdf              = other._cdf;
}

DiscretePowerLawDistribution::DiscretePowerLawDistribution(
    const vector<int>& sampleData, int xParameter,
    double alphaPrecision, DistributionType distributionType)
{
    _state            = InputValidator(sampleData, xParameter, distributionType);
    _alphaPrecision   = alphaPrecision;
    _distributionType = distributionType;

    if (_state == DistributionState::Valid)
    {
        if (distributionType == DistributionType::LeftBounded)
        {
            _xMin       = xParameter;
            _xMax       = VectorUtilities::Max(sampleData);
            _alpha      = EstimateAlpha(sampleData, _xMin, alphaPrecision);
            _sampleSize = VectorUtilities::NumberOfGreaterOrEqual(sampleData, _xMin);
        }
        else
        {
            _xMin       = 1;
            _xMax       = xParameter;
            _alpha      = EstimateAlpha(sampleData, _xMin, _xMax, alphaPrecision);
            _sampleSize = VectorUtilities::NumberOfLowerOrEqual(sampleData, _xMax);
        }
        PrecalculateCDF();
        _ksStatistic = CalculateKSStatistic(sampleData);
    }
}

DiscretePowerLawDistribution::DiscretePowerLawDistribution(
    const vector<int>& sampleData, double alphaPrecision,
    DistributionType distributionType, int smallestInterval)
{
    _state            = InputValidator(sampleData);
    _alphaPrecision   = alphaPrecision;
    _distributionType = distributionType;

    if (_state == DistributionState::Valid)
    {
        if (distributionType == DistributionType::LeftBounded)
        {
            _xMin       = EstimateLowerBound(sampleData, alphaPrecision);
            _xMax       = VectorUtilities::Max(sampleData);
            _alpha      = EstimateAlpha(sampleData, _xMin, alphaPrecision);
            _sampleSize = VectorUtilities::NumberOfGreaterOrEqual(sampleData, _xMin);
        }
        else
        {
            _xMin       = 1;
            _xMax       = EstimateUpperBound(sampleData, alphaPrecision, smallestInterval);
            _alpha      = EstimateAlpha(sampleData, _xMin, _xMax, alphaPrecision);
            _sampleSize = VectorUtilities::NumberOfLowerOrEqual(sampleData, _xMax);
        }
        PrecalculateCDF();
        _ksStatistic = CalculateKSStatistic(sampleData);
    }
}

/* ---- input validation -------------------------------------------------- */

DistributionState DiscretePowerLawDistribution::InputValidator(const vector<int>& data)
{
    return !data.empty() ? DistributionState::Valid : DistributionState::NoInput;
}

DistributionState DiscretePowerLawDistribution::InputValidator(
    const vector<int>& data, int xParameter, DistributionType distributionType)
{
    if (data.empty())
        return DistributionState::NoInput;

    if (distributionType == DistributionType::LeftBounded)
    {
        if (xParameter >= VectorUtilities::Max(data))
            return DistributionState::InvalidInput;
    }
    else
    {
        if (xParameter <= VectorUtilities::Min(data))
            return DistributionState::InvalidInput;
    }
    return DistributionState::Valid;
}

/* ---- estimation -------------------------------------------------------- */

double DiscretePowerLawDistribution::EstimateAlpha(
    const vector<int>& data, int xMin, double precision)
{
    /* Pre-compute the weighted log-sum once (avoids recomputation per alpha).
     *
     * L(α) = -n * log ζ(α, xMin)  -  α * Σ log(x)
     *
     * The second term is linear in α so it does not depend on the zeta
     * function.  Pre-computing it lets Brent's method call only one zeta
     * evaluation per function evaluation rather than two.
     */
    double logXSum = 0.0;
    int    n       = 0;
    for (int x : data)
        if (x >= xMin) { logXSum += std::log(static_cast<double>(x)); ++n; }

    if (n == 0) return 2.0; /* fallback */

    const double dn = static_cast<double>(n);

    auto logL = [&](double alpha) -> double {
        return -dn * std::log(real_hurwitz_zeta(alpha, static_cast<double>(xMin)))
               - alpha * logXSum;
    };

    /* Brent's method over a generous bracket [1.01, 20.0].
     * Use precision/10 as the absolute tolerance so the returned estimate is
     * consistent with the requested alpha_precision.  */
    return BrentMaximize(logL, 1.01, 20.0, precision * 0.1);
}

double DiscretePowerLawDistribution::EstimateAlpha(
    const vector<int>& data, int xMin, int xMax, double precision)
{
    double logXSum = 0.0;
    int    n       = 0;
    for (int x : data)
        if (x >= xMin && x <= xMax)
        {
            logXSum += std::log(static_cast<double>(x));
            ++n;
        }

    if (n == 0) return 2.0;

    const double dn = static_cast<double>(n);

    auto logL = [&](double alpha) -> double {
        const double denom = real_hurwitz_zeta(alpha, 1.0)
                           - real_hurwitz_zeta(alpha, static_cast<double>(xMax + 1));
        if (denom <= 0.0) return -std::numeric_limits<double>::infinity();
        return -dn * std::log(denom) - alpha * logXSum;
    };

    return BrentMaximize(logL, 1.01, 20.0, precision * 0.1);
}

int DiscretePowerLawDistribution::EstimateLowerBound(
    const vector<int>& data, double precision)
{
    const int minElem = VectorUtilities::Min(data);
    const int maxElem = VectorUtilities::Max(data);

    double minKS     = numeric_limits<double>::infinity();
    int    xMinEstim = minElem;

    for (int x = minElem; x < maxElem; ++x)
    {
        const DiscretePowerLawDistribution model(data, x, precision,
                                                 DistributionType::LeftBounded);
        const double ks = model.GetKSStatistic();
        if (ks < minKS)
        {
            minKS     = ks;
            xMinEstim = x;
        }
        else
        {
            break;
        }
    }
    return max(xMinEstim, 1);
}

int DiscretePowerLawDistribution::EstimateUpperBound(
    const vector<int>& data, double precision, int smallestInterval)
{
    const int minElem = 1 + smallestInterval;
    const int maxElem = VectorUtilities::Max(data);

    vector<double> ksValues;
    ksValues.reserve(static_cast<size_t>(maxElem - minElem));
    for (int x = minElem; x < maxElem; ++x)
    {
        const DiscretePowerLawDistribution model(data, x, precision,
                                                 DistributionType::RightBounded);
        ksValues.push_back(model.GetKSStatistic());
    }

    return VectorUtilities::IndexOfMin(ksValues) + minElem;
}

/* ---- log-likelihoods --------------------------------------------------- */

double DiscretePowerLawDistribution::CalculateLogLikelihoodLeftBounded(
    const vector<int>& data, double alpha, int xMin)
{
    const double n = static_cast<double>(
        VectorUtilities::NumberOfGreaterOrEqual(data, xMin));

    double logXSum = 0.0;
    for (int x : data)
        if (x >= xMin)
            logXSum += std::log(static_cast<double>(x));

    return -n * std::log(real_hurwitz_zeta(alpha, xMin)) - alpha * logXSum;
}

double DiscretePowerLawDistribution::CalculateLogLikelihoodRightBounded(
    const vector<int>& data, double alpha, int xMax)
{
    const double n = static_cast<double>(
        VectorUtilities::NumberOfLowerOrEqual(data, xMax));

    double logXSum = 0.0;
    for (int x : data)
        if (x >= 1 && x <= xMax)
            logXSum += std::log(static_cast<double>(x));

    const double denom = real_hurwitz_zeta(alpha, 1)
                       - real_hurwitz_zeta(alpha, 1 + xMax);
    return -n * std::log(denom) - alpha * logXSum;
}

/* ---- CDF --------------------------------------------------------------- */

void DiscretePowerLawDistribution::PrecalculateCDF()
{
    _cdf.reserve(static_cast<size_t>(_xMax - _xMin + 1));
    for (int x = _xMin; x <= _xMax; ++x)
    {
        double val = (_distributionType == DistributionType::LeftBounded)
                         ? CalculateCDF(x, _alpha, _xMin)
                         : CalculateCDF(x, _alpha, _xMin, _xMax);
        _cdf.push_back(val);
    }
}

double DiscretePowerLawDistribution::CalculateCDF(int x, double alpha, int xMin)
{
    if (x >= xMin)
        return real_hurwitz_zeta(alpha, x) / real_hurwitz_zeta(alpha, xMin);
    return 1.0;
}

double DiscretePowerLawDistribution::CalculateCDF(int x, double alpha,
                                                   int xMin, int xMax)
{
    if (x >= xMin && x <= xMax)
    {
        const double num = real_hurwitz_zeta(alpha, x)
                         - real_hurwitz_zeta(alpha, 1 + xMax);
        const double den = real_hurwitz_zeta(alpha, xMin)
                         - real_hurwitz_zeta(alpha, 1 + xMax);
        return num / den;
    }
    else if (x < xMin)
        return 1.0;
    else
        return 0.0;
}

int DiscretePowerLawDistribution::BinarySearch(int l, int r, double x) const
{
    while (l <= r)
    {
        const int    mid  = l + (r - l) / 2;
        const double cdf  = GetCDF(mid);
        const double rCdf = GetCDF(mid + 1);
        const double lCdf = GetCDF(mid - 1);

        if (x < lCdf && x > rCdf)
            return (x > cdf) ? mid - 1 : mid;
        else if (cdf < x)
            r = mid - 1;
        else
            l = mid + 1;
    }
    return -1;
}

/* ---- random variate generation ---------------------------------------- */

vector<int> DiscretePowerLawDistribution::GenerateRandomSequence(int n) const
{
    vector<int> seq;
    seq.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
        seq.push_back(GenerateRandomSample());
    return seq;
}

int DiscretePowerLawDistribution::GenerateRandomSample() const
{
    if (_state != DistributionState::Valid)
        return numeric_limits<int>::min();

    const double r = RandomGen::GetUniform01();
    int x1, x2;
    x2 = _xMin;
    do
    {
        x1 = x2;
        x2 = 2 * x1;
    } while (GetCDF(x2) >= r);

    return BinarySearch(x1, x2, r);
}

/* ---- KS statistic ------------------------------------------------------ */

double DiscretePowerLawDistribution::CalculateKSStatistic(const vector<int>& data) const
{
    if (!StateIsValid())
        return numeric_limits<double>::infinity();

    DiscreteEmpiricalDistribution empirical(data, _xMin, _xMax);

    double maxDiff = 0.0;
    for (int x = _xMin; x <= _xMax; ++x)
    {
        const double diff = std::abs(empirical.GetCDF(x) - GetCDF(x));
        if (diff > maxDiff)
            maxDiff = diff;
    }
    return maxDiff;
}

/* ---- accessors --------------------------------------------------------- */

double DiscretePowerLawDistribution::GetPDF(int x) const
{
    if (_state != DistributionState::Valid)
        return numeric_limits<double>::quiet_NaN();
    return std::pow(static_cast<double>(x), -_alpha)
           / real_hurwitz_zeta(_alpha, _xMin);
}

double DiscretePowerLawDistribution::GetCDF(int x) const
{
    if (_state != DistributionState::Valid)
        return numeric_limits<double>::quiet_NaN();
    if (x >= _xMin && x <= _xMax)
        return _cdf[static_cast<size_t>(x - _xMin)];
    else if (x < _xMin)
        return 1.0;
    else
        return 0.0;
}

double DiscretePowerLawDistribution::GetKSStatistic() const
{
    if (_state != DistributionState::Valid)
        return numeric_limits<double>::infinity();
    return _ksStatistic;
}

double DiscretePowerLawDistribution::GetAlpha() const
{
    if (_state != DistributionState::Valid)
        return numeric_limits<double>::quiet_NaN();
    return _alpha;
}

double DiscretePowerLawDistribution::GetAlphaPrecision() const
{
    return _alphaPrecision;
}

int DiscretePowerLawDistribution::GetXMin() const
{
    if (_state != DistributionState::Valid)
        return numeric_limits<int>::min();
    return _xMin;
}

int DiscretePowerLawDistribution::GetXMax() const
{
    if (_state != DistributionState::Valid)
        return numeric_limits<int>::min();
    return _xMax;
}

double DiscretePowerLawDistribution::GetStandardError() const
{
    if (_state != DistributionState::Valid)
        return numeric_limits<double>::quiet_NaN();
    return GetStandardError(_sampleSize);
}

double DiscretePowerLawDistribution::GetStandardError(int sampleSize) const
{
    return (_alpha - 1.0) / static_cast<double>(sampleSize);
}

double DiscretePowerLawDistribution::GetLogLikelihood(const vector<int>& data) const
{
    if (_distributionType == DistributionType::LeftBounded)
        return CalculateLogLikelihoodLeftBounded(data, _alpha, _xMin);
    else
        return CalculateLogLikelihoodRightBounded(data, _alpha, _xMax);
}

bool DiscretePowerLawDistribution::StateIsValid() const
{
    return _state == DistributionState::Valid;
}

DistributionState DiscretePowerLawDistribution::GetState() const
{
    return _state;
}

DistributionType DiscretePowerLawDistribution::GetDistributionType() const
{
    return _distributionType;
}

std::string DiscretePowerLawDistribution::GetDistributionTypeStr() const
{
    switch (_distributionType)
    {
        case DistributionType::LeftBounded:  return "Left bounded (Type I)";
        case DistributionType::RightBounded: return "Right bounded (Type II)";
        default:                             return "<unknown>";
    }
}

DiscretePowerLawDistribution DiscretePowerLawDistribution::FromParameters(
    double alpha, int xMin, int xMax)
{
    DiscretePowerLawDistribution dist;
    dist._alpha            = alpha;
    dist._xMin             = xMin;
    dist._xMax             = xMax;
    dist._alphaPrecision   = 0.01;
    dist._sampleSize       = 0;
    dist._ksStatistic      = 0.0;
    dist._distributionType = DistributionType::LeftBounded;
    dist._state            = DistributionState::Valid;
    dist.PrecalculateCDF();
    return dist;
}

/* =========================================================================
 * SyntheticPowerLawGenerator
 * ======================================================================= */

SyntheticPowerLawGenerator::SyntheticPowerLawGenerator(
    const DiscretePowerLawDistribution& model,
    const vector<int>& sampleData,
    SyntheticGeneratorMode mode)
: _powerLawDistribution(model)
{
    _sampleDataSize = static_cast<int>(sampleData.size());
    _mode           = mode;

    if (mode == SyntheticGeneratorMode::SemiParametric)
    {
        _nonModelData = sampleData;
        if (model.GetDistributionType() == DistributionType::LeftBounded)
            VectorUtilities::RemoveGreaterOrEqual(_nonModelData, model.GetXMin());
        else
            VectorUtilities::RemoveLowerOrEqual(_nonModelData, model.GetXMax());

        _modelSampleProbability =
            1.0 - static_cast<double>(_nonModelData.size())
                / static_cast<double>(_sampleDataSize);
    }
    else
    {
        _modelSampleProbability = 1.0;
    }
}

int SyntheticPowerLawGenerator::SampleFromData() const
{
    const int idx = RandomGen::GetInt(static_cast<int>(_nonModelData.size()) - 1);
    return _nonModelData[static_cast<size_t>(idx)];
}

vector<int> SyntheticPowerLawGenerator::SampleFromData(int n) const
{
    vector<int> v;
    v.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
        v.push_back(SampleFromData());
    return v;
}

vector<int> SyntheticPowerLawGenerator::GenerateSynthetic() const
{
    vector<int> synthetic;
    synthetic.reserve(static_cast<size_t>(_sampleDataSize));

    const int modelN = static_cast<int>(
        std::floor(_modelSampleProbability * _sampleDataSize));
    VectorUtilities::Insert(synthetic,
                            _powerLawDistribution.GenerateRandomSequence(modelN));
    VectorUtilities::Insert(synthetic,
                            SampleFromData(_sampleDataSize - modelN));
    return synthetic;
}

double SyntheticPowerLawGenerator::MeasureKsStatisticOfReplica() const
{
    const vector<int> syn            = GenerateSynthetic();
    const DistributionType dtype     = _powerLawDistribution.GetDistributionType();
    const double       prec          = _powerLawDistribution.GetAlphaPrecision();

    if (_mode == SyntheticGeneratorMode::SemiParametric)
    {
        const DiscretePowerLawDistribution model(syn, prec, dtype);
        return model.GetKSStatistic();
    }
    else
    {
        const int xParam = (dtype == DistributionType::LeftBounded)
                               ? _powerLawDistribution.GetXMin()
                               : _powerLawDistribution.GetXMax();
        const DiscretePowerLawDistribution model(syn, xParam, prec, dtype);
        return model.GetKSStatistic();
    }
}
