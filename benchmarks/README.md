# Benchmarks

This directory contains scripts for benchmarking `poweRcpp` against `poweRlaw`
to evaluate performance gaps and guide redesign decisions—particularly around
parallelization of the goodness-of-fit bootstrapping.

## Scripts

| Script | Description |
|--------|-------------|
| `benchmark_gof.R` | Compares `poweRcpp::powerlaw_gof` vs `poweRlaw::bootstrap` using `bench::press` |

## Requirements

Install the required packages once before running any benchmark:

```r
install.packages("bench")
install.packages("poweRlaw")
devtools::install_github("alrobles/poweRcpp")
```

## Running

From the repository root:

```sh
Rscript benchmarks/benchmark_gof.R
```

## Output

The script prints:

1. **Full `bench::press` table** — median time, memory allocation, garbage
   collections, and iteration counts for each backend.
2. **Human-readable summary** — which package is faster and by how much
   (speedup factor).
3. Saves an `.rds` file to `benchmarks/results/benchmark_gof_results.rds`
   for later inspection with `readRDS()`.

### Example output

```
Running benchmark on 8 core(s)

=== Full Benchmark Results ===
# A tibble: 2 × 13
  backend  min   median `itr/sec` mem_alloc `gc/sec` ...
  <chr>    <bch:tm> <bch:tm>  <dbl> <bch:byt>    <dbl>
1 poweRcpp ...      ...        ...   ...           ...
2 poweRlaw ...      ...        ...   ...           ...

=== Summary: poweRcpp vs poweRlaw ===
  poweRcpp  median time : Xs
  poweRlaw  median time : Xs

  Result: poweRcpp is X.XXx FASTER than poweRlaw
```

## Interpreting Results

- **`median` time** is the most reliable metric for comparing throughput.
- **`mem_alloc`** shows total memory allocated during the benchmark; lower
  values indicate less GC pressure.
- A **speedup > 1** means `poweRcpp` is faster than `poweRlaw`; a speedup of
  exactly 1 means equal performance; values < 1 indicate `poweRlaw` is faster
  and `poweRcpp` still needs improvement in that area.

## Context

The goal of `poweRcpp` is to outperform `poweRlaw` by leveraging C++ via Rcpp
and parallelizing both the core distribution routines and the goodness-of-fit
bootstrapping. Use these benchmarks to measure progress as the implementation
is refined.
