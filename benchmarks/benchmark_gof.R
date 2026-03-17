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
##
## Usage:
##   Rscript benchmarks/benchmark_gof.R
##   Results are printed to the console and saved to
##   benchmarks/results/benchmark_gof_results.rds
## ============================================================

library(poweRcpp)
library(poweRlaw)

# Ensure results directory exists
results_dir <- file.path("benchmarks", "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Generate test data
set.seed(42)
test_data <- rpois(10000, lambda = 10)

# Use all available cores for the poweRlaw bootstrap
n_threads <- parallel::detectCores()
cat(sprintf("Running benchmark on %d core(s)\n\n", n_threads))

# Goodness-of-Fit Benchmark
profiles <- bench::press(
  backend = c("poweRcpp", "poweRlaw"),
  {
    if (backend == "poweRcpp") {
      # Run poweRcpp goodness-of-fit; result captured to include return value
      # in timing (prevents the compiler from optimizing the call away)
      result <- poweRcpp::powerlaw_gof(test_data, replicas = 10)
    } else {
      # Run poweRlaw goodness-of-fit
      pl <- poweRlaw::displ$new(test_data)
      result <- bootstrap(pl, no_of_sims = 10, threads = n_threads)
    }
  }
)

# Print full bench results
cat("\n=== Full Benchmark Results ===\n")
print(profiles)

# Print a human-readable summary comparing the two backends
cat("\n=== Summary: poweRcpp vs poweRlaw ===\n")
timings <- profiles[, c("backend", "median", "mem_alloc")]
timings$median_sec <- as.numeric(timings$median)

poweRcpp_row <- timings[timings$backend == "poweRcpp", ]
poweRlaw_row <- timings[timings$backend == "poweRlaw", ]

cat(sprintf("  poweRcpp  median time : %s\n", format(poweRcpp_row$median)))
cat(sprintf("  poweRlaw  median time : %s\n", format(poweRlaw_row$median)))

speedup <- poweRlaw_row$median_sec / poweRcpp_row$median_sec
if (speedup > 1) {
  cat(sprintf("\n  Result: poweRcpp is %.2fx FASTER than poweRlaw\n", speedup))
} else if (speedup < 1) {
  cat(sprintf("\n  Result: poweRlaw is %.2fx FASTER than poweRcpp\n", 1 / speedup))
} else {
  cat("\n  Result: poweRcpp and poweRlaw have equal performance\n")
}

# Save results for later inspection
saveRDS(profiles, file = file.path(results_dir, "benchmark_gof_results.rds"))
cat(sprintf("\nResults saved to %s\n",
            file.path(results_dir, "benchmark_gof_results.rds")))
