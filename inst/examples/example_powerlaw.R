## ============================================================
## poweRcpp – usage examples
## ============================================================
## This script demonstrates the main functions provided by the
## poweRcpp package.  Run after installing the package:
##
##   install.packages("devtools")
##   devtools::install_github("alrobles/poweRcpp")
##   library(poweRcpp)
## ============================================================

library(poweRcpp)

## ------------------------------------------------------------
## 1.  Fit a power law with all parameters estimated from data
## ------------------------------------------------------------
set.seed(42)

## Simulate 500 observations from a power-law (alpha = 2.5, x_min = 3)
sim_data <- powerlaw_generate(500, alpha = 2.5, x_min = 3L)

fit <- fit_powerlaw(sim_data)
cat("=== Fitted model (all parameters estimated) ===\n")
cat("Alpha            :", fit$alpha, "\n")
cat("x_min            :", fit$x_min, "\n")
cat("x_max            :", fit$x_max, "\n")
cat("KS statistic     :", round(fit$ks_statistic, 4), "\n")
cat("Standard error   :", round(fit$standard_error, 4), "\n")
cat("Log-likelihood   :", round(fit$log_likelihood, 2), "\n")
cat("Tail observations:", fit$n_tail, "\n")
cat("Type             :", fit$type, "\n")
cat("Valid            :", fit$valid, "\n\n")

## ------------------------------------------------------------
## 2.  Fit with a known x_min
## ------------------------------------------------------------
fit_fixed <- fit_powerlaw(sim_data, x_min = 3L)
cat("=== Fitted model (x_min = 3 fixed) ===\n")
cat("Alpha            :", fit_fixed$alpha, "\n")
cat("x_min            :", fit_fixed$x_min, "\n\n")

## ------------------------------------------------------------
## 3.  Evaluate PDF and survival function
## ------------------------------------------------------------
x_vals <- 1L:20L

pdf_vals <- powerlaw_pdf(x_vals, alpha = 2.5, x_min = 1L)
cdf_vals <- powerlaw_cdf(x_vals, alpha = 2.5, x_min = 1L)

cat("=== PMF and survival function (alpha=2.5, x_min=1) ===\n")
cat(sprintf("x = %2d  PMF = %.4f  Survival = %.4f\n",
            x_vals, pdf_vals, cdf_vals), sep = "")

## Check that PMF sums to ≈ 1 over the tail
cat("Sum of PMF[1..20]:", round(sum(pdf_vals), 4),
    "(not exactly 1 because the tail is infinite)\n\n")

## ------------------------------------------------------------
## 4.  Generate random samples and visualise
## ------------------------------------------------------------
set.seed(123)
rng_samples <- powerlaw_generate(1000L, alpha = 2.5, x_min = 1L)

cat("=== Random sample summary (n=1000, alpha=2.5) ===\n")
cat("Min :", min(rng_samples), "\n")
cat("Mean:", round(mean(rng_samples), 2), "\n")
cat("Max :", max(rng_samples), "\n\n")

## Log-log histogram
breaks <- 10^seq(0, log10(max(rng_samples) + 1), length.out = 20)
hist(rng_samples, breaks = breaks, plot = FALSE)

## ------------------------------------------------------------
## 5.  Goodness-of-fit test (small number of replicas for speed)
## ------------------------------------------------------------
cat("=== Goodness-of-fit test (200 replicas) ===\n")
gof <- powerlaw_gof(sim_data, replicas = 200L)
cat("p-value      :", gof$p_value, "\n")
cat("KS statistic :", round(gof$ks_statistic, 4), "\n")
cat("Alpha used   :", gof$alpha, "\n")
cat("x_min used   :", gof$x_min, "\n")
cat("Replicas     :", gof$replicas, "\n\n")

cat("A p-value >= 0.10 provides no evidence against the power-law hypothesis.\n")
