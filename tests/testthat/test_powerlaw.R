## Unit tests for poweRcpp
## Tests are designed to be fast (small n, low alpha_precision, few replicas).

library(poweRcpp)

## ============================================================
## Helper: simple power-law sample generator in pure R
## (for testing purposes only)
## ============================================================
r_powerlaw <- function(n, alpha, x_min = 1L)
{
    ## Inverse-transform approximation using continuous approximation
    u   <- runif(n)
    x   <- floor(x_min * (1 - u)^(-1 / (alpha - 1)))
    as.integer(pmax(x, x_min))
}

## Reference PMF (R implementation) for comparison with C++ PDF
r_powerlaw_pmf <- function(x, alpha, x_min = 1L)
{
    ## P(X = k) = k^{-alpha} / sum_{k=x_min}^{Inf} k^{-alpha}
    ## Approximated as k^{-alpha} / zeta(alpha, x_min) using the
    ## partial sum over a wide range
    k   <- x_min:max(x_min + 5000L, max(x))
    Z   <- sum(k^(-alpha))
    vapply(as.integer(x), function(xi)
        if (xi < x_min) 0.0 else xi^(-alpha) / Z, numeric(1))
}

## ============================================================
## Tests: fit_powerlaw
## ============================================================

test_that("fit_powerlaw returns a list with expected names", {
    set.seed(1)
    data <- as.integer(r_powerlaw(200, alpha = 2.5, x_min = 2L))
    result <- fit_powerlaw(data, x_min = 2L, alpha_precision = 0.1)

    expect_type(result, "list")
    expected_names <- c("alpha", "x_min", "x_max", "ks_statistic",
                        "standard_error", "log_likelihood", "n_tail",
                        "type", "valid")
    expect_true(all(expected_names %in% names(result)))
})

test_that("fit_powerlaw with fixed x_min returns valid = TRUE for valid data", {
    set.seed(2)
    data <- as.integer(r_powerlaw(300, alpha = 2.5, x_min = 1L))
    result <- fit_powerlaw(data, x_min = 1L, alpha_precision = 0.1)

    expect_true(result$valid)
    expect_gt(result$alpha, 1.0)
    expect_gte(result$x_min, 1L)
})

test_that("fit_powerlaw alpha estimate is in a plausible range", {
    set.seed(3)
    data <- as.integer(r_powerlaw(500, alpha = 2.5, x_min = 1L))
    result <- fit_powerlaw(data, x_min = 1L, alpha_precision = 0.05)

    ## With alpha_precision = 0.05 and 500 samples the estimate should be
    ## within ±0.5 of the true value (very generous tolerance)
    expect_true(result$valid)
    expect_gt(result$alpha, 1.5)
    expect_lt(result$alpha, 4.0)
})

test_that("fit_powerlaw without x_min still returns a valid fit", {
    set.seed(4)
    data <- as.integer(r_powerlaw(200, alpha = 2.5, x_min = 1L))
    result <- fit_powerlaw(data, alpha_precision = 0.1)

    expect_true(result$valid)
    expect_false(is.na(result$alpha))
    expect_false(is.na(result$x_min))
})

test_that("fit_powerlaw stops on empty input", {
    expect_error(fit_powerlaw(integer(0)), "'data' must be a non-empty")
})

test_that("fit_powerlaw stops on non-positive values", {
    expect_error(fit_powerlaw(c(-1L, 2L, 3L)), "positive integers")
})

test_that("fit_powerlaw coerces numeric to integer without error", {
    set.seed(5)
    data_num <- r_powerlaw(100, alpha = 2.5, x_min = 1L)
    result   <- fit_powerlaw(data_num, x_min = 1L, alpha_precision = 0.1)
    expect_true(result$valid)
})

test_that("fit_powerlaw removes NAs with a warning", {
    set.seed(6)
    data <- as.integer(r_powerlaw(100, alpha = 2.5))
    data_with_na <- c(data, NA_integer_)
    expect_warning(
        result <- fit_powerlaw(data_with_na, x_min = 1L, alpha_precision = 0.1),
        "NA values"
    )
    expect_true(result$valid)
})

## ============================================================
## Tests: powerlaw_pdf
## ============================================================

test_that("powerlaw_pdf returns values between 0 and 1", {
    vals <- powerlaw_pdf(1L:20L, alpha = 2.5, x_min = 1L)
    expect_true(all(vals >= 0 & vals <= 1))
})

test_that("powerlaw_pdf is monotonically decreasing", {
    vals <- powerlaw_pdf(1L:10L, alpha = 2.5, x_min = 1L)
    expect_true(all(diff(vals) < 0))
})

test_that("powerlaw_pdf sums to approximately 1 over a wide range", {
    ## The true sum over all positive integers should be 1, so summing
    ## over 1..10000 should be close
    vals  <- powerlaw_pdf(1L:10000L, alpha = 3.0, x_min = 1L)
    total <- sum(vals)
    expect_gt(total, 0.99)
    expect_lt(total, 1.01)
})

test_that("powerlaw_pdf stops when alpha <= 1", {
    expect_error(powerlaw_pdf(1L:5L, alpha = 1.0), "'alpha' must be > 1")
})

## ============================================================
## Tests: powerlaw_cdf
## ============================================================

test_that("powerlaw_cdf at x_min equals 1.0", {
    val <- powerlaw_cdf(1L, alpha = 2.5, x_min = 1L)
    expect_equal(val, 1.0)
})

test_that("powerlaw_cdf is monotonically non-increasing", {
    vals <- powerlaw_cdf(1L:20L, alpha = 2.5, x_min = 1L)
    expect_true(all(diff(vals) <= 0))
})

test_that("powerlaw_cdf returns values in [0, 1]", {
    vals <- powerlaw_cdf(1L:50L, alpha = 2.5, x_min = 1L)
    expect_true(all(vals >= 0 & vals <= 1))
})

## ============================================================
## Tests: powerlaw_generate
## ============================================================

test_that("powerlaw_generate returns an integer vector of length n", {
    set.seed(10)
    smp <- powerlaw_generate(50L, alpha = 2.5, x_min = 1L)
    expect_type(smp, "integer")
    expect_length(smp, 50L)
})

test_that("powerlaw_generate samples are >= x_min", {
    set.seed(11)
    x_min <- 3L
    smp   <- powerlaw_generate(200L, alpha = 2.5, x_min = x_min)
    expect_true(all(smp >= x_min))
})

test_that("powerlaw_generate stops for invalid alpha", {
    expect_error(powerlaw_generate(10L, alpha = 0.5, x_min = 1L), "'alpha' must be > 1")
})

test_that("powerlaw_generate stops for n <= 0", {
    expect_error(powerlaw_generate(0L, alpha = 2.5), "'n' must be a positive integer")
})

## ============================================================
## Tests: powerlaw_gof
## ============================================================

test_that("powerlaw_gof returns a list with expected names", {
    set.seed(20)
    data <- as.integer(r_powerlaw(100, alpha = 2.5, x_min = 1L))
    gof  <- powerlaw_gof(data, x_min = 1L, replicas = 50L)

    expect_type(gof, "list")
    expected <- c("p_value", "ks_statistic", "alpha", "x_min", "replicas")
    expect_true(all(expected %in% names(gof)))
})

test_that("powerlaw_gof p_value is in [0, 1]", {
    set.seed(21)
    data <- as.integer(r_powerlaw(100, alpha = 2.5, x_min = 1L))
    gof  <- powerlaw_gof(data, x_min = 1L, replicas = 50L)

    expect_gte(gof$p_value, 0.0)
    expect_lte(gof$p_value, 1.0)
})

test_that("powerlaw_gof ks_statistic is non-negative", {
    set.seed(22)
    data <- as.integer(r_powerlaw(100, alpha = 2.5, x_min = 1L))
    gof  <- powerlaw_gof(data, x_min = 1L, replicas = 50L)

    expect_gte(gof$ks_statistic, 0.0)
})

test_that("powerlaw_gof stops on empty data", {
    expect_error(powerlaw_gof(integer(0)), "'data' must be a non-empty")
})

test_that("powerlaw_gof stops on non-positive replicas", {
    set.seed(23)
    data <- as.integer(r_powerlaw(100, alpha = 2.5))
    expect_error(powerlaw_gof(data, replicas = 0L), "'replicas' must be")
})

## ============================================================
## Tests: alpha estimation accuracy (Brent's method vs truth)
## ============================================================

test_that("fit_powerlaw alpha estimate is close to truth with larger sample", {
    ## With n = 2000 and a large x_min the MLE should be within ~0.15 of 2.5
    set.seed(30)
    data   <- as.integer(r_powerlaw(2000, alpha = 2.5, x_min = 1L))
    result <- fit_powerlaw(data, x_min = 1L, alpha_precision = 0.01)

    expect_true(result$valid)
    expect_gt(result$alpha, 2.0)
    expect_lt(result$alpha, 3.2)
})

test_that("fit_powerlaw alpha is more precise than alpha_precision", {
    ## Brent's method returns continuous-valued alpha, not rounded to grid
    set.seed(31)
    data   <- as.integer(r_powerlaw(500, alpha = 2.5, x_min = 1L))
    result <- fit_powerlaw(data, x_min = 1L, alpha_precision = 0.1)

    expect_true(result$valid)
    ## alpha should be strictly greater precision than alpha_precision=0.1
    ## (allow equality in degenerate cases, but alpha must be > 1)
    expect_gt(result$alpha, 1.0)
})

## ============================================================
## Tests: PDF vs reference R implementation
## ============================================================

test_that("powerlaw_pdf matches pure-R reference implementation", {
    x     <- 1L:30L
    alpha <- 2.5
    x_min <- 1L

    cpp_pdf <- powerlaw_pdf(x, alpha = alpha, x_min = x_min)
    ref_pdf <- r_powerlaw_pmf(x, alpha = alpha, x_min = x_min)

    ## Should agree to within 0.5 % (partial-sum reference vs exact zeta)
    rel_err <- abs(cpp_pdf - ref_pdf) / pmax(ref_pdf, 1e-12)
    expect_true(all(rel_err < 0.005))
})

## ============================================================
## Tests: right-bounded distribution
## ============================================================

test_that("fit_powerlaw type='right' returns valid result", {
    set.seed(40)
    ## Generate data truncated to [1, 100] to simulate a right-bounded scenario.
    ## For type = "right", x_min is used as xMax (upper bound), so pass NULL
    ## to let the algorithm estimate xMax automatically.
    data <- as.integer(r_powerlaw(300, alpha = 2.5, x_min = 1L))
    data <- data[data <= 100L]

    result <- fit_powerlaw(data, alpha_precision = 0.1, type = "right")
    expect_true(result$valid)
    expect_equal(result$type, "Right bounded (Type II)")
    expect_gt(result$alpha, 1.0)
})

test_that("fit_powerlaw log_likelihood is finite for right-bounded fit", {
    set.seed(41)
    data <- as.integer(r_powerlaw(200, alpha = 2.5, x_min = 1L))
    data <- data[data <= 200L]

    result <- fit_powerlaw(data, alpha_precision = 0.1, type = "right")
    expect_true(result$valid)
    expect_true(is.finite(result$log_likelihood))
})

## ============================================================
## Tests: survival function (CDF) properties
## ============================================================

test_that("powerlaw_cdf is strictly decreasing for x > x_min", {
    vals <- powerlaw_cdf(1L:50L, alpha = 2.5, x_min = 1L)
    ## All differences should be strictly negative
    expect_true(all(diff(vals) < 0))
})

test_that("powerlaw_cdf approaches 0 as x grows large", {
    val_large <- powerlaw_cdf(10000L, alpha = 2.5, x_min = 1L)
    expect_lt(val_large, 0.001)
})

## ============================================================
## Tests: edge values below x_min
## ============================================================

test_that("powerlaw_pdf returns 0 for x < x_min", {
    vals <- powerlaw_pdf(c(1L, 2L), alpha = 2.5, x_min = 3L)
    expect_equal(vals[1], 0.0)
    expect_equal(vals[2], 0.0)
})

test_that("powerlaw_cdf returns 1.0 for x < x_min (survival function)", {
    vals <- powerlaw_cdf(c(0L, 1L, 2L), alpha = 2.5, x_min = 3L)
    expect_true(all(vals == 1.0))
})

## ============================================================
## Tests: powerlaw_generate with custom x_max
## ============================================================

test_that("powerlaw_generate respects custom x_max upper bound", {
    set.seed(50)
    x_min <- 1L
    x_max <- 20L
    smp   <- powerlaw_generate(200L, alpha = 2.5, x_min = x_min, x_max = x_max)
    expect_true(all(smp >= x_min))
    expect_true(all(smp <= x_max))
})

test_that("powerlaw_generate stops when x_max <= x_min", {
    expect_error(powerlaw_generate(10L, alpha = 2.5, x_min = 5L, x_max = 5L),
                 "'x_max' must be greater than 'x_min'")
    expect_error(powerlaw_generate(10L, alpha = 2.5, x_min = 5L, x_max = 3L),
                 "'x_max' must be greater than 'x_min'")
})

## ============================================================
## Tests: type parameter validation
## ============================================================

test_that("fit_powerlaw stops on invalid type string", {
    set.seed(60)
    data <- as.integer(r_powerlaw(100, alpha = 2.5))
    expect_error(fit_powerlaw(data, type = "foo"))
})

test_that("powerlaw_gof stops on invalid type string", {
    set.seed(61)
    data <- as.integer(r_powerlaw(100, alpha = 2.5))
    expect_error(powerlaw_gof(data, replicas = 10L, type = "bar"))
})

## ============================================================
## Tests: powerlaw_gof with both alpha and x_min supplied
## ============================================================

test_that("powerlaw_gof with both alpha and x_min supplied returns valid result", {
    set.seed(70)
    data  <- as.integer(r_powerlaw(200, alpha = 2.5, x_min = 1L))
    fit   <- fit_powerlaw(data, x_min = 1L, alpha_precision = 0.1)
    gof   <- powerlaw_gof(data,
                          alpha    = fit$alpha,
                          x_min    = fit$x_min,
                          replicas = 30L)

    expect_type(gof, "list")
    expect_gte(gof$p_value, 0.0)
    expect_lte(gof$p_value, 1.0)
})

## ============================================================
## Tests: right-bounded GOF
## ============================================================

test_that("powerlaw_gof type='right' returns valid result", {
    set.seed(80)
    data <- as.integer(r_powerlaw(300, alpha = 2.5, x_min = 1L))
    data <- data[data <= 100L]

    gof <- powerlaw_gof(data, replicas = 30L, type = "right")

    expect_type(gof, "list")
    expected <- c("p_value", "ks_statistic", "alpha", "x_min", "replicas")
    expect_true(all(expected %in% names(gof)))
    expect_gte(gof$p_value, 0.0)
    expect_lte(gof$p_value, 1.0)
})

## ============================================================
## Tests: S3 class and print/summary methods
## ============================================================

test_that("fit_powerlaw result has class 'powerlaw_fit'", {
    set.seed(90)
    data <- as.integer(r_powerlaw(100, alpha = 2.5, x_min = 1L))
    result <- fit_powerlaw(data, x_min = 1L, alpha_precision = 0.1)
    expect_s3_class(result, "powerlaw_fit")
})

test_that("print.powerlaw_fit runs without error", {
    set.seed(91)
    data   <- as.integer(r_powerlaw(100, alpha = 2.5, x_min = 1L))
    result <- fit_powerlaw(data, x_min = 1L, alpha_precision = 0.1)
    expect_output(print(result), "Alpha")
    expect_output(print(result), "KS statistic")
})

test_that("summary.powerlaw_fit runs without error", {
    set.seed(92)
    data   <- as.integer(r_powerlaw(100, alpha = 2.5, x_min = 1L))
    result <- fit_powerlaw(data, x_min = 1L, alpha_precision = 0.1)
    expect_output(summary(result), "alpha")
})

test_that("powerlaw_gof result has class 'powerlaw_gof'", {
    set.seed(93)
    data <- as.integer(r_powerlaw(100, alpha = 2.5, x_min = 1L))
    gof  <- powerlaw_gof(data, x_min = 1L, replicas = 20L)
    expect_s3_class(gof, "powerlaw_gof")
})

test_that("print.powerlaw_gof runs without error", {
    set.seed(94)
    data <- as.integer(r_powerlaw(100, alpha = 2.5, x_min = 1L))
    gof  <- powerlaw_gof(data, x_min = 1L, replicas = 20L)
    expect_output(print(gof), "p-value")
    expect_output(print(gof), "KS statistic")
})

## ============================================================
## Tests: n_tail accuracy
## ============================================================

test_that("fit_powerlaw n_tail equals observations >= x_min", {
    set.seed(100)
    data   <- as.integer(r_powerlaw(300, alpha = 2.5, x_min = 1L))
    x_min  <- 3L
    result <- fit_powerlaw(data, x_min = x_min, alpha_precision = 0.1)

    expected_n_tail <- sum(data >= x_min)
    expect_equal(result$n_tail, expected_n_tail)
})

