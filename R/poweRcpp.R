#' Fit a discrete power-law distribution
#'
#' Estimates the exponent \eqn{\alpha} and, optionally, the lower cut-off
#' \eqn{x_{\min}} of a discrete power-law distribution using the maximum
#' likelihood estimator and Kolmogorov-Smirnov criterion described in Clauset,
#' Shalizi & Newman (2009).  Alpha is found via Brent's method applied to the
#' log-likelihood, which is faster than a fixed-step grid search.
#'
#' @param data A positive integer vector of observations.
#' @param x_min Integer lower cut-off.  When \code{NULL} (default) it is
#'   estimated from the data by minimising the KS statistic.
#' @param alpha_precision Numeric tolerance used when searching for the optimal
#'   alpha (default \code{0.01}).  Smaller values give higher precision but
#'   require more Hurwitz zeta evaluations.
#' @param type Distribution type.  \code{"left"} (default) fits a left-bounded
#'   (Type I) power law; \code{"right"} fits a right-bounded (Type II) power law.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{alpha}{Estimated power-law exponent.}
#'   \item{x_min}{Estimated (or supplied) lower cut-off.}
#'   \item{x_max}{Upper cut-off.}
#'   \item{ks_statistic}{Kolmogorov-Smirnov goodness-of-fit statistic.}
#'   \item{standard_error}{Standard error of the alpha estimate.}
#'   \item{log_likelihood}{Log-likelihood at the estimated parameters.}
#'   \item{n_tail}{Number of observations in the power-law tail.}
#'   \item{type}{Distribution type string.}
#'   \item{valid}{Logical: \code{TRUE} if fitting succeeded.}
#' }
#'
#' @references
#' Clauset, A., Shalizi, C.R. & Newman, M.E.J. (2009).
#' Power-law distributions in empirical data.
#' \emph{SIAM Review}, 51(4), 661-703.
#' \doi{10.1137/070710111}
#'
#' @examples
#' set.seed(42)
#' data <- sample(1:500, 200, replace = TRUE, prob = (1:500)^(-2.5))
#' result <- fit_powerlaw(data)
#' cat("Estimated alpha:", result$alpha, "\n")
#' cat("Estimated x_min:", result$x_min, "\n")
#'
#' @export
fit_powerlaw <- function(data,
                         x_min           = NULL,
                         alpha_precision = 0.01,
                         type            = c("left", "right"))
{
    type <- match.arg(type)

    if (!is.integer(data))
        data <- as.integer(data)
    if (any(is.na(data)))
    {
        warning("NA values found in 'data'. They will be removed before fitting.")
        data <- data[!is.na(data)]
    }
    if (length(data) == 0L)
        stop("'data' must be a non-empty integer vector.")
    if (any(data <= 0L))
        stop("All values in 'data' must be positive integers.")

    fit_powerlaw_cpp(
        data,
        x_min           = if (is.null(x_min)) NULL else as.integer(x_min),
        alpha_precision = alpha_precision,
        type            = type
    )
}

#' Discrete power-law probability mass function
#'
#' Evaluates the PMF \eqn{P(X = x) = x^{-\alpha} / \zeta(\alpha, x_{\min})}
#' of a left-bounded discrete power-law distribution.
#'
#' @param x     Integer vector of values.
#' @param alpha Power-law exponent (\eqn{\alpha > 1}).
#' @param x_min Lower cut-off (positive integer, default \code{1}).
#'
#' @return Numeric vector of probability mass values.
#'
#' @examples
#' powerlaw_pdf(1:10, alpha = 2.5, x_min = 1)
#'
#' @export
powerlaw_pdf <- function(x, alpha, x_min = 1L)
{
    if (!is.integer(x)) x <- as.integer(x)
    if (alpha <= 1) stop("'alpha' must be > 1.")
    if (x_min < 1)  stop("'x_min' must be >= 1.")
    powerlaw_pdf_cpp(x, alpha = alpha, x_min = as.integer(x_min))
}

#' Discrete power-law survival function (complementary CDF)
#'
#' Evaluates \eqn{P(X \geq x)} for a left-bounded discrete power-law
#' distribution with exponent \code{alpha} and lower cut-off \code{x_min}.
#'
#' @param x     Integer vector of values.
#' @param alpha Power-law exponent (\eqn{\alpha > 1}).
#' @param x_min Lower cut-off (positive integer, default \code{1}).
#'
#' @return Numeric vector of survival-function values in \eqn{[0, 1]}.
#'
#' @examples
#' powerlaw_cdf(1:10, alpha = 2.5, x_min = 1)
#'
#' @export
powerlaw_cdf <- function(x, alpha, x_min = 1L)
{
    if (!is.integer(x)) x <- as.integer(x)
    if (alpha <= 1) stop("'alpha' must be > 1.")
    if (x_min < 1)  stop("'x_min' must be >= 1.")
    powerlaw_cdf_cpp(x, alpha = alpha, x_min = as.integer(x_min))
}

#' Generate random samples from a discrete power-law distribution
#'
#' Generates integer variates from a discrete power-law distribution using
#' inversion sampling on the pre-computed CDF.
#'
#' @param n     Number of samples to generate.
#' @param alpha Power-law exponent (\eqn{\alpha > 1}).
#' @param x_min Lower cut-off (positive integer, default \code{1}).
#' @param x_max Upper cut-off.  When \code{NULL} (default), uses
#'   \code{x_min * 1000}.
#'
#' @return Integer vector of length \code{n}.
#'
#' @examples
#' set.seed(1)
#' samples <- powerlaw_generate(100, alpha = 2.5, x_min = 1)
#' hist(samples, breaks = 20, main = "Power-law samples")
#'
#' @export
powerlaw_generate <- function(n, alpha, x_min = 1L, x_max = NULL)
{
    if (!is.numeric(n) || n <= 0)
        stop("'n' must be a positive integer.")
    if (alpha <= 1)
        stop("'alpha' must be > 1.")
    if (x_min < 1)
        stop("'x_min' must be >= 1.")

    powerlaw_generate_cpp(
        as.integer(n),
        alpha = alpha,
        x_min = as.integer(x_min),
        x_max = if (is.null(x_max)) NULL else as.integer(x_max)
    )
}

#' Goodness-of-fit test for a discrete power-law distribution
#'
#' Computes a p-value for the hypothesis that \code{data} are drawn from a
#' discrete power-law distribution.  A semi-parametric bootstrap procedure
#' (Clauset, Shalizi & Newman 2009) is used: synthetic datasets are generated
#' by combining a parametric power-law tail with non-parametric resampling of
#' the body, and the fraction of replicas whose KS statistic exceeds the
#' observed value is returned as the p-value.  The bootstrap loop is
#' parallelised with \pkg{RcppParallel}.
#'
#' @param data     Positive integer vector of observations.
#' @param alpha    Power-law exponent.  When \code{NULL} (default) it is
#'   estimated from the data.
#' @param x_min    Lower cut-off.  When \code{NULL} (default) it is estimated
#'   from the data.
#' @param replicas Number of bootstrap replicas (default \code{1000}).
#'   Values \eqn{\geq 1000} are recommended for stable p-values.
#' @param type     Distribution type (\code{"left"} or \code{"right"}).
#'
#' @return A named list with elements:
#' \describe{
#'   \item{p_value}{Bootstrap p-value.}
#'   \item{ks_statistic}{Observed KS statistic.}
#'   \item{alpha}{Alpha used for the test.}
#'   \item{x_min}{xMin used for the test.}
#'   \item{replicas}{Number of replicas used.}
#' }
#'
#' @references
#' Clauset, A., Shalizi, C.R. & Newman, M.E.J. (2009).
#' Power-law distributions in empirical data.
#' \emph{SIAM Review}, 51(4), 661-703.
#' \doi{10.1137/070710111}
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' data <- sample(1:500, 200, replace = TRUE, prob = (1:500)^(-2.5))
#' gof  <- powerlaw_gof(data, replicas = 500)
#' cat("p-value:", gof$p_value, "\n")
#' }
#'
#' @export
powerlaw_gof <- function(data,
                         alpha    = NULL,
                         x_min   = NULL,
                         replicas = 1000L,
                         type     = c("left", "right"))
{
    type <- match.arg(type)

    if (!is.integer(data)) data <- as.integer(data)
    if (any(is.na(data)))
    {
        warning("NA values found in 'data'. They will be removed.")
        data <- data[!is.na(data)]
    }
    if (length(data) == 0L)
        stop("'data' must be a non-empty integer vector.")
    if (any(data <= 0L))
        stop("All values in 'data' must be positive integers.")
    if (replicas < 1L)
        stop("'replicas' must be a positive integer.")

    powerlaw_gof_cpp(
        data,
        alpha    = if (is.null(alpha)) NULL else as.double(alpha),
        x_min    = if (is.null(x_min)) NULL else as.integer(x_min),
        replicas = as.integer(replicas),
        type     = type
    )
}
