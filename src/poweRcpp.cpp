/**
 * @file poweRcpp.cpp
 * @brief Rcpp/RcppParallel interface for discrete power-law fitting.
 *
 * Exposes fitting, PDF/CDF evaluation, random generation, and parallel
 * goodness-of-fit testing to R via Rcpp export attributes.
 *
 * The goodness-of-fit bootstrap (powerlaw_gof_cpp) uses RcppParallel's
 * parallelReduce so that each synthetic replica is generated and evaluated
 * on a separate thread, dramatically reducing computation time for large
 * numbers of replicas.
 */

// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "DiscreteDistributions.h"
#include "Zeta.h"
using namespace Rcpp;
using namespace RcppParallel;

/* -----------------------------------------------------------------------
 * Helpers
 * --------------------------------------------------------------------- */

/// Convert an R integer vector to std::vector<int>, discarding NAs.
static std::vector<int> to_int_vector(const IntegerVector& x)
{
    std::vector<int> v;
    v.reserve(static_cast<size_t>(x.size()));
    for (int i = 0; i < x.size(); ++i)
        if (x[i] != NA_INTEGER)
            v.push_back(x[i]);
    return v;
}

/* -----------------------------------------------------------------------
 * Exported functions
 * --------------------------------------------------------------------- */

//' Fit a discrete power-law distribution
//'
//' Estimates the exponent (alpha) and, optionally, the lower cut-off (xMin)
//' of a discrete power-law distribution using maximum likelihood and the
//' Kolmogorov-Smirnov criterion described in Clauset, Shalizi & Newman (2009).
//'
//' @param data     Integer vector of observations (positive integers).
//' @param x_min    Integer lower cut-off. If \code{NA} (default) it is
//'                 estimated from the data via KS minimisation.
//' @param alpha_precision Numeric grid step for the alpha search
//'                 (default 0.01).
//' @param type     Distribution type: \code{"left"} (Type I, default) or
//'                 \code{"right"} (Type II, right-bounded).
//'
//' @return A named list with elements:
//'   \describe{
//'     \item{alpha}{Estimated power-law exponent.}
//'     \item{x_min}{Estimated (or supplied) lower cut-off.}
//'     \item{x_max}{Upper cut-off (maximum observed value for Type I).}
//'     \item{ks_statistic}{Kolmogorov-Smirnov goodness-of-fit statistic.}
//'     \item{standard_error}{Standard error of the alpha estimate.}
//'     \item{log_likelihood}{Log-likelihood at the estimated parameters.}
//'     \item{n_tail}{Number of observations in the power-law tail.}
//'     \item{type}{Distribution type string.}
//'     \item{valid}{Logical: \code{TRUE} if fitting succeeded.}
//'   }
//'
//' @examples
//' set.seed(42)
//' data <- sample(1:500, 200, replace = TRUE, prob = (1:500)^(-2.5))
//' result <- fit_powerlaw(data)
//' print(result$alpha)
//'
//' @export
// [[Rcpp::export]]
List fit_powerlaw_cpp(IntegerVector data,
                      Nullable<int>  x_min        = R_NilValue,
                      double         alpha_precision = 0.01,
                      std::string    type          = "left")
{
    if (data.size() == 0)
        stop("'data' must be a non-empty integer vector.");
    if (alpha_precision <= 0.0 || alpha_precision > 1.0)
        stop("'alpha_precision' must be in (0, 1].");

    const DistributionType dtype = (type == "right")
                                       ? DistributionType::RightBounded
                                       : DistributionType::LeftBounded;

    std::vector<int> v = to_int_vector(data);

    DiscretePowerLawDistribution model =
        x_min.isNull()
            ? DiscretePowerLawDistribution(v, alpha_precision, dtype)
            : DiscretePowerLawDistribution(v, as<int>(x_min), alpha_precision, dtype);

    const bool valid = model.StateIsValid();

    return List::create(
        Named("alpha")          = valid ? model.GetAlpha()          : NA_REAL,
        Named("x_min")          = valid ? model.GetXMin()           : NA_INTEGER,
        Named("x_max")          = valid ? model.GetXMax()           : NA_INTEGER,
        Named("ks_statistic")   = valid ? model.GetKSStatistic()    : NA_REAL,
        Named("standard_error") = valid ? model.GetStandardError()  : NA_REAL,
        Named("log_likelihood") = valid ? model.GetLogLikelihood(v) : NA_REAL,
        Named("n_tail")         = valid ? (int) v.size()            : NA_INTEGER,
        Named("type")           = model.GetDistributionTypeStr(),
        Named("valid")          = valid
    );
}

//' Discrete power-law probability density function
//'
//' Returns P(X = x) for a discrete power-law distribution with the given
//' parameters.
//'
//' @param x     Integer vector of values at which to evaluate the PMF.
//' @param alpha Power-law exponent (alpha > 1).
//' @param x_min Lower cut-off (positive integer, default 1).
//'
//' @return Numeric vector of probability mass values.
//'
//' @examples
//' powerlaw_pdf(1:10, alpha = 2.5, x_min = 1)
//'
//' @export
// [[Rcpp::export]]
NumericVector powerlaw_pdf_cpp(IntegerVector x, double alpha, int x_min = 1)
{
    if (alpha <= 1.0) stop("'alpha' must be > 1.");
    if (x_min < 1)   stop("'x_min' must be >= 1.");

    /* Compute PDF directly: P(X=x) = x^(-alpha) / zeta(alpha, x_min) */
    const double denom = real_hurwitz_zeta(alpha, static_cast<double>(x_min));

    NumericVector out(x.size());
    for (int i = 0; i < x.size(); ++i)
    {
        if (x[i] == NA_INTEGER || x[i] < x_min)
            out[i] = (x[i] == NA_INTEGER) ? NA_REAL : 0.0;
        else
            out[i] = std::pow(static_cast<double>(x[i]), -alpha) / denom;
    }
    return out;
}

//' Discrete power-law cumulative distribution function
//'
//' Returns the survival function \eqn{P(X \geq x)} for a discrete power-law
//' distribution with the given parameters.
//'
//' @param x     Integer vector of values at which to evaluate the CDF.
//' @param alpha Power-law exponent (alpha > 1).
//' @param x_min Lower cut-off (positive integer, default 1).
//'
//' @return Numeric vector of survival-function values.
//'
//' @examples
//' powerlaw_cdf(1:10, alpha = 2.5, x_min = 1)
//'
//' @export
// [[Rcpp::export]]
NumericVector powerlaw_cdf_cpp(IntegerVector x, double alpha, int x_min = 1)
{
    if (alpha <= 1.0) stop("'alpha' must be > 1.");
    if (x_min < 1)   stop("'x_min' must be >= 1.");

    /* Survival function: P(X >= x) = zeta(alpha, x) / zeta(alpha, x_min) */
    const double denom = real_hurwitz_zeta(alpha, static_cast<double>(x_min));

    NumericVector out(x.size());
    for (int i = 0; i < x.size(); ++i)
    {
        if (x[i] == NA_INTEGER)
            out[i] = NA_REAL;
        else if (x[i] < x_min)
            out[i] = 1.0;
        else
            out[i] = real_hurwitz_zeta(alpha, static_cast<double>(x[i])) / denom;
    }
    return out;
}

//' Generate random samples from a discrete power-law distribution
//'
//' @param n     Number of samples to generate.
//' @param alpha Power-law exponent (alpha > 1).
//' @param x_min Lower cut-off (positive integer, default 1).
//' @param x_max Upper cut-off. If \code{NA} (default), uses \code{x_min * 1000}.
//'
//' @return Integer vector of length \code{n}.
//'
//' @examples
//' set.seed(1)
//' samples <- powerlaw_generate(100, alpha = 2.5, x_min = 1)
//' hist(samples, breaks = 20)
//'
//' @export
// [[Rcpp::export]]
IntegerVector powerlaw_generate_cpp(int n, double alpha, int x_min = 1,
                                    Nullable<int> x_max = R_NilValue)
{
    if (n <= 0)       stop("'n' must be a positive integer.");
    if (alpha <= 1.0) stop("'alpha' must be > 1.");
    if (x_min < 1)   stop("'x_min' must be >= 1.");

    const int xmax_val = x_max.isNull() ? x_min * 1000 : as<int>(x_max);
    if (xmax_val <= x_min)
        stop("'x_max' must be greater than 'x_min'.");

    DiscretePowerLawDistribution model =
        DiscretePowerLawDistribution::FromParameters(alpha, x_min, xmax_val);

    RandomGen::Seed();
    const std::vector<int> seq = model.GenerateRandomSequence(n);

    IntegerVector out(seq.begin(), seq.end());
    return out;
}

//' Goodness-of-fit test for a discrete power-law distribution
//'
//' Computes a p-value for the hypothesis that \code{data} follow a discrete
//' power-law distribution using a semi-parametric bootstrap procedure
//' (Clauset, Shalizi & Newman 2009).  The bootstrap is parallelised with
//' \pkg{RcppParallel}.
//'
//' @param data     Integer vector of observations.
//' @param alpha    Power-law exponent.  If \code{NA} (default) it is estimated
//'                 from \code{data}.
//' @param x_min    Lower cut-off.  If \code{NA} (default) it is estimated
//'                 from \code{data}.
//' @param replicas Number of bootstrap replicas (default 1000).
//' @param type     Distribution type: \code{"left"} (default) or
//'                 \code{"right"}.
//'
//' @return A named list with elements:
//'   \describe{
//'     \item{p_value}{Bootstrap p-value.}
//'     \item{ks_statistic}{Observed KS statistic.}
//'     \item{alpha}{Alpha used for the test.}
//'     \item{x_min}{xMin used for the test.}
//'     \item{replicas}{Number of replicas used.}
//'   }
//'
//' @examples
//' \dontrun{
//' set.seed(42)
//' data <- sample(1:500, 200, replace = TRUE, prob = (1:500)^(-2.5))
//' gof  <- powerlaw_gof(data, replicas = 500)
//' print(gof$p_value)
//' }
//'
//' @export
// [[Rcpp::export]]
List powerlaw_gof_cpp(IntegerVector data,
                      Nullable<double> alpha    = R_NilValue,
                      Nullable<int>    x_min    = R_NilValue,
                      int              replicas = 1000,
                      std::string      type     = "left")
{
    if (data.size() == 0)
        stop("'data' must be a non-empty integer vector.");
    if (replicas <= 0)
        stop("'replicas' must be a positive integer.");

    const DistributionType dtype = (type == "right")
                                       ? DistributionType::RightBounded
                                       : DistributionType::LeftBounded;

    std::vector<int> v = to_int_vector(data);

    /* Fit the model (or use supplied parameters) */
    DiscretePowerLawDistribution model =
        x_min.isNull()
            ? DiscretePowerLawDistribution(v, 0.01, dtype)
            : DiscretePowerLawDistribution(v, as<int>(x_min), 0.01, dtype);

    if (!model.StateIsValid())
        stop("Power-law model fitting failed. Check that data contains positive integers.");

    const double observed_ks = model.GetKSStatistic();

    /* Run parallel bootstrap */
    RandomGen::Seed();
    SyntheticPowerLawGenerator gen(model, v, SyntheticGeneratorMode::SemiParametric);

    /* Use a parallel reduce to collect KS statistics across replicas */
    std::vector<double> ksValues(static_cast<size_t>(replicas));

    /* Custom parallel worker using parallelFor */
    struct KsWorker : public Worker
    {
        const SyntheticPowerLawGenerator& gen;
        std::vector<double>& ks;

        KsWorker(const SyntheticPowerLawGenerator& g, std::vector<double>& k)
            : gen(g), ks(k) {}

        void operator()(std::size_t begin, std::size_t end) override
        {
            for (std::size_t i = begin; i < end; ++i)
                ks[i] = gen.MeasureKsStatisticOfReplica();
        }
    };

    KsWorker worker(gen, ksValues);
    parallelFor(0, static_cast<std::size_t>(replicas), worker);

    /* Count replicas where synthetic KS >= observed KS */
    int larger = 0;
    for (double ks : ksValues)
        if (ks >= observed_ks)
            ++larger;

    const double p_value = static_cast<double>(larger)
                         / static_cast<double>(replicas);

    return List::create(
        Named("p_value")      = p_value,
        Named("ks_statistic") = observed_ks,
        Named("alpha")        = model.GetAlpha(),
        Named("x_min")        = model.GetXMin(),
        Named("replicas")     = replicas
    );
}
