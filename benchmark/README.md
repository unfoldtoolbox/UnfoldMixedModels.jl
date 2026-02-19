# UnfoldMixedModels.jl Benchmarks

This directory contains benchmark scripts for comparing different implementations of parameter extraction in UnfoldMixedModels.jl.

## Running Benchmarks

### Setup

First, install BenchmarkTools:

```bash
cd benchmark
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Run Parameter Extraction Benchmarks

```bash
cd benchmark
julia --project=.. parameter_extraction_benchmarks.jl
```

Or from the repository root:

```bash
julia --project=. benchmark/parameter_extraction_benchmarks.jl
```

## What is Benchmarked

The `parameter_extraction_benchmarks.jl` script compares:

1. **Raw extraction**: `allpars` vs `pbstrtbl` for getting all parameters from a fit collection
2. **tidyβ**: Using `allpars` internally vs direct `MixedModels.tidyβ`
3. **tidyσs**: Using `allpars` internally vs direct `MixedModels.tidyσs`
4. **make_estimate**: Using `pbstrtbl` vs calling `coef`/`ranef`/`ranefcorr` separately

## Interpreting Results

The benchmark reports time ratios:
- **Ratio < 1.0**: First approach is faster
- **Ratio ≈ 1.0-1.2**: Comparable performance (within 20%)
- **Ratio > 1.2**: First approach is slower (by more than 20%)

## Example Output

```
================================================================================
Benchmark 1: Raw Parameter Extraction - allpars vs pbstrtbl
================================================================================

Testing allpars()...
allpars: 1.234 ms

Testing pbstrtbl()...
pbstrtbl: 1.456 ms

Speed comparison:
  pbstrtbl / allpars = 1.18x
  ≈ Comparable performance (< 20% difference)
```

## Contributing

When adding new optimizations or changing parameter extraction logic, please:
1. Run these benchmarks before and after your changes
2. Document any significant performance changes in your PR
3. Update this README if you add new benchmarks
