## ============================================================
## Benchmark: poweRcpp vs poweRlaw – Goodness-of-Fit
## ============================================================
## This script compares the performance of the goodness-of-fit
## (bootstrapping) routines in poweRcpp and poweRlaw.
##
## Requirements (install once before running):
##   install.packages("profvis")
##   install.packages("bench")
##   devtools::install_github("alrobles/poweRcpp")
##   install.packages("poweRlaw")
## ============================================================

library(poweRcpp)
library(poweRlaw)

# Generate test data
data <- rpois(10000, lambda = 10)

# Use all available cores for the poweRlaw bootstrap
n_threads <- parallel::detectCores()

# Goodness-of-Fit Benchmark
profiles <- bench::press(
  backend = c("poweRcpp", "poweRlaw"),
  {
    if (backend == "poweRcpp") {
      # Run poweRcpp goodness-of-fit
      result <- poweRcpp::powerlaw_gof(data, replicas = 10)
    } else {
      # Run poweRlaw goodness-of-fit
      pl <- poweRlaw::displ$new(data)
      result <- bootstrap(pl, no_of_sims = 10, threads = n_threads)
    }
  }
)
print(profiles)
