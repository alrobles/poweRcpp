/**
 * @file DiscreteDistributions.h
 * @brief Discrete power-law distribution classes.
 *
 * Adapted from the alrobles/powerlaw repository.
 * Reference: Clauset, Shalizi & Newman (2009) https://arxiv.org/abs/0706.1062
 *
 * Alpha MLE uses Brent's method (super-linear convergence, ~50 zeta function
 * evaluations) instead of a fixed grid search, giving a ~10× speed-up for the
 * default alpha_precision = 0.01 case.
 */
#pragma once
#include <vector>
#include <string>
#include <limits>
#include "RandomGen.h"

/**
 * Stores an empirical CDF for a discrete distribution with truncated tail.
 */
class DiscreteEmpiricalDistribution
{
private:
    int _xMin, _xMax;
    std::vector<double> _cdf;

    void PrecalculateCDF(const std::vector<int>& sortedTailSample);

public:
    DiscreteEmpiricalDistribution(const std::vector<int>& sampleData,
                                  int xMin, int xMax);

    /// Returns the survival function (complementary CDF) at x.
    double GetCDF(int x) const;
};

enum class DistributionType
{
    LeftBounded,  ///< Type I: lower-bounded distribution (xMin estimated or fixed)
    RightBounded  ///< Type II: upper-bounded distribution (xMax estimated or fixed)
};

enum class DistributionState
{
    Valid,        ///< Distribution was successfully fitted
    NoInput,      ///< Empty data vector supplied
    InvalidInput  ///< Parameter inconsistent with data
};

enum class SyntheticGeneratorMode
{
    SemiParametric, ///< Mix: parametric tail + non-parametric body
    FullParametric  ///< Fully parametric generation
};

/**
 * Discrete power-law distribution: parameter estimation, PDF/CDF, and
 * random variate generation.
 *
 * Implements the maximum-likelihood estimator described in:
 * Clauset, Shalizi & Newman (2009), SIAM Review 51(4):661–703.
 */
class DiscretePowerLawDistribution
{
private:
    DistributionType  _distributionType;
    DistributionState _state;
    double _alpha;
    double _ksStatistic;
    double _alphaPrecision;
    int _xMin, _xMax;
    int _sampleSize;
    std::vector<double> _cdf;

    static DistributionState InputValidator(const std::vector<int>& data);
    static DistributionState InputValidator(const std::vector<int>& data,
                                            int xParameter,
                                            DistributionType distributionType);

    static double EstimateAlpha(const std::vector<int>& data,
                                int xMin,
                                double precision = 0.01);

    static double EstimateAlpha(const std::vector<int>& data,
                                int xMin, int xMax,
                                double precision = 0.01);

    static int EstimateLowerBound(const std::vector<int>& data,
                                  double precision = 0.01);

    static int EstimateUpperBound(const std::vector<int>& data,
                                  double precision = 0.01,
                                  int smallestInterval = 20);

    static double CalculateLogLikelihoodLeftBounded(const std::vector<int>& data,
                                                    double alpha, int xMin);

    static double CalculateLogLikelihoodRightBounded(const std::vector<int>& data,
                                                     double alpha, int xMax);

    static double CalculateCDF(int x, double alpha, int xMin);
    static double CalculateCDF(int x, double alpha, int xMin, int xMax);

    int    BinarySearch(int l, int r, double x) const;
    double GetStandardError(int sampleSize) const;
    void   PrecalculateCDF();

    /* Private default constructor – used only by FromParameters. */
    DiscretePowerLawDistribution() = default;

public:
    DiscretePowerLawDistribution(const DiscretePowerLawDistribution& other);

    /**
     * Fit a power-law with a fixed xParameter.
     *
     * @param sampleData  Input data.
     * @param xParameter  Fixed xMin (LeftBounded) or xMax (RightBounded).
     * @param alphaPrecision  Grid precision for alpha search.
     * @param distributionType  Type I or Type II distribution.
     */
    DiscretePowerLawDistribution(const std::vector<int>& sampleData,
                                 int xParameter,
                                 double alphaPrecision = 0.01,
                                 DistributionType distributionType = DistributionType::LeftBounded);

    /**
     * Fit a power-law with all parameters estimated from data.
     *
     * @param sampleData  Input data.
     * @param alphaPrecision  Grid precision for alpha search.
     * @param distributionType  Type I or Type II distribution.
     * @param smallestInterval  Minimum allowed xMax − xMin for Type II.
     */
    explicit DiscretePowerLawDistribution(const std::vector<int>& sampleData,
                                          double alphaPrecision = 0.01,
                                          DistributionType distributionType = DistributionType::LeftBounded,
                                          int smallestInterval = 20);

    std::vector<int> GenerateRandomSequence(int n) const;
    int GenerateRandomSample() const;

    double GetPDF(int x) const;
    double GetCDF(int x) const;
    double GetKSStatistic() const;
    double GetAlpha() const;
    double GetAlphaPrecision() const;
    double GetStandardError() const;
    double GetLogLikelihood(const std::vector<int>& data) const;
    int    GetXMin() const;
    int    GetXMax() const;
    bool   StateIsValid() const;
    DistributionState  GetState() const;
    DistributionType   GetDistributionType() const;
    std::string        GetDistributionTypeStr() const;

    /**
     * Compute the KS statistic between the fitted model and the empirical
     * distribution of @p data.  This method is public so that callers who
     * build a model via FromParameters can still obtain the KS statistic.
     */
    double CalculateKSStatistic(const std::vector<int>& data) const;

    /**
     * Construct with explicitly supplied alpha, xMin, and xMax (no fitting).
     * Intended for PDF/CDF evaluation and random generation when the
     * parameters are already known.  Call CalculateKSStatistic() separately
     * if the KS statistic is required.
     *
     * @param alpha            Power-law exponent (> 1).
     * @param xMin             Lower cut-off (>= 1).
     * @param xMax             Upper cut-off.
     * @param distributionType Left- or right-bounded (default: LeftBounded).
     */
    static DiscretePowerLawDistribution FromParameters(
        double alpha, int xMin, int xMax,
        DistributionType distributionType = DistributionType::LeftBounded);
};

/**
 * Generates synthetic power-law replicas for goodness-of-fit testing.
 */
class SyntheticPowerLawGenerator
{
private:
    DiscretePowerLawDistribution _powerLawDistribution;
    SyntheticGeneratorMode       _mode;
    std::vector<int>             _nonModelData;
    double _modelSampleProbability;
    int    _sampleDataSize;

    int              SampleFromData() const;
    std::vector<int> SampleFromData(int n) const;

public:
    SyntheticPowerLawGenerator(const DiscretePowerLawDistribution& model,
                               const std::vector<int>& sampleData,
                               SyntheticGeneratorMode mode = SyntheticGeneratorMode::SemiParametric);

    std::vector<int> GenerateSynthetic() const;
    double           MeasureKsStatisticOfReplica() const;
};
