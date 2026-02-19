"""
Benchmark script to compare different parameter extraction implementations:
1. tidyβ/tidyσs with allpars vs direct MixedModels functions
2. pbstrtbl vs allpars for extracting all parameters
3. make_estimate with pbstrtbl vs hypothetical allpars-based approach

Run this script with:
    julia --project=. benchmark/parameter_extraction_benchmarks.jl
"""

using BenchmarkTools
using UnfoldMixedModels
using Unfold
using MixedModels
using StatsModels
using DataFrames
using StableRNGs
using UnfoldSim
using Printf

println("="^80)
println("Parameter Extraction Benchmarks")
println("="^80)
println()

# Create test data for benchmarking
println("Setting up test data...")
data, evts = UnfoldSim.predef_2x2(
    StableRNG(1);
    return_epoched = false,
    n_subjects = 5,
    noiselevel = 1,
    signalsize = 10,
    n_items = 16,
)

subj_idx = [parse(Int, split(string(s), 'S')[2]) for s in evts.subject]
evts.latency .+= size(data, 1) .* (subj_idx .- 1)
data = reshape(data, 1, :)
data = vcat(data, data)  # 2 channels

transform!(evts, :subject => categorical => :subject)

# Cut the data into epochs
data_e, times = Unfold.epoch(data = data, tbl = evts, τ = (-1.0, 1.9), sfreq = 10)
evts_e, data_e = Unfold.drop_missing_epochs(copy(evts), data_e)

# Fit the model
println("Fitting model...")
f = @formula 0 ~ 1 + A + B + (1 + A + B | subject)

m = fit(
    UnfoldModel,
    f,
    evts_e,
    data_e,
    times,
    contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding()),
    show_progress = false,
)

println("Model fitted successfully!")
println("Model type: ", typeof(m))
println("Number of fits: ", length(modelfit(m).fits))
println()

# ============================================================================
# Benchmark 1: allpars vs pbstrtbl for raw extraction
# ============================================================================
println("="^80)
println("Benchmark 1: Raw Parameter Extraction - allpars vs pbstrtbl")
println("="^80)

println("\nTesting allpars()...")
b1_allpars = @benchmark MixedModels.allpars($(modelfit(m))) samples=100
println("allpars: ", BenchmarkTools.prettytime(median(b1_allpars).time))

println("\nTesting pbstrtbl()...")
b1_pbstrtbl = @benchmark MixedModels.pbstrtbl($(modelfit(m))) samples=100
println("pbstrtbl: ", BenchmarkTools.prettytime(median(b1_pbstrtbl).time))

ratio1 = median(b1_pbstrtbl).time / median(b1_allpars).time
println("\nSpeed comparison:")
println("  pbstrtbl / allpars = ", @sprintf("%.2fx", ratio1))
if ratio1 < 1.0
    println("  ✓ pbstrtbl is ", @sprintf("%.1f%%", (1-ratio1)*100), " faster")
else
    println("  ✗ pbstrtbl is ", @sprintf("%.1f%%", (ratio1-1)*100), " slower")
end
println()

# ============================================================================
# Benchmark 2: tidyβ implementations
# ============================================================================
println("="^80)
println("Benchmark 2: tidyβ - allpars-based vs direct MixedModels.tidyβ")
println("="^80)

println("\nTesting current tidyβ (using allpars)...")
b2_tidy_allpars = @benchmark MixedModels.tidyβ($m) samples=50
println("tidyβ with allpars: ", BenchmarkTools.prettytime(median(b2_tidy_allpars).time))

println("\nTesting direct MixedModels.tidyβ on collection...")
b2_tidy_direct = @benchmark MixedModels.tidyβ($(modelfit(m))) samples=50
println("MixedModels.tidyβ direct: ", BenchmarkTools.prettytime(median(b2_tidy_direct).time))

ratio2 = median(b2_tidy_allpars).time / median(b2_tidy_direct).time
println("\nSpeed comparison:")
println("  allpars-based / direct = ", @sprintf("%.2fx", ratio2))
if ratio2 < 1.2
    println("  ≈ Comparable performance (< 20% difference)")
elseif ratio2 < 1.0
    println("  ✓ allpars-based is ", @sprintf("%.1f%%", (1-ratio2)*100), " faster")
else
    println("  ✗ allpars-based is ", @sprintf("%.1f%%", (ratio2-1)*100), " slower")
end
println()

# ============================================================================
# Benchmark 3: tidyσs implementations
# ============================================================================
println("="^80)
println("Benchmark 3: tidyσs - allpars-based vs direct MixedModels.tidyσs")
println("="^80)

println("\nTesting current tidyσs (using allpars)...")
b3_tidy_allpars = @benchmark MixedModels.tidyσs($m) samples=50
println("tidyσs with allpars: ", BenchmarkTools.prettytime(median(b3_tidy_allpars).time))

println("\nTesting direct MixedModels.tidyσs on collection...")
b3_tidy_direct = @benchmark MixedModels.tidyσs($(modelfit(m))) samples=50
println("MixedModels.tidyσs direct: ", BenchmarkTools.prettytime(median(b3_tidy_direct).time))

ratio3 = median(b3_tidy_allpars).time / median(b3_tidy_direct).time
println("\nSpeed comparison:")
println("  allpars-based / direct = ", @sprintf("%.2fx", ratio3))
if ratio3 < 1.2
    println("  ≈ Comparable performance (< 20% difference)")
elseif ratio3 < 1.0
    println("  ✓ allpars-based is ", @sprintf("%.1f%%", (1-ratio3)*100), " faster")
else
    println("  ✗ allpars-based is ", @sprintf("%.1f%%", (ratio3-1)*100), " slower")
end
println()

# ============================================================================
# Benchmark 4: make_estimate - pbstrtbl vs allpars-based approach
# ============================================================================
println("="^80)
println("Benchmark 4: make_estimate - pbstrtbl vs calling individual functions")
println("="^80)

println("\nTesting current make_estimate (using pbstrtbl)...")
b4_make_pbstrtbl = @benchmark Unfold.make_estimate($m) samples=50
println("make_estimate with pbstrtbl: ", BenchmarkTools.prettytime(median(b4_make_pbstrtbl).time))

# Alternative: manually call coef, ranef, ranefcorr (what make_estimate used to do)
function make_estimate_old_style(m)
    coefs = StatsModels.coef(m)
    ranef_sigma = MixedModels.ranef(m)
    corrs = UnfoldMixedModels.ranefcorr(m)
    # Just return for timing purposes
    return (coefs, ranef_sigma, corrs)
end

println("\nTesting alternative (calling coef/ranef/ranefcorr separately)...")
b4_make_separate = @benchmark make_estimate_old_style($m) samples=50
println("make_estimate with separate calls: ", BenchmarkTools.prettytime(median(b4_make_separate).time))

ratio4 = median(b4_make_pbstrtbl).time / median(b4_make_separate).time
println("\nSpeed comparison:")
println("  pbstrtbl / separate calls = ", @sprintf("%.2fx", ratio4))
if ratio4 < 1.0
    println("  ✓ pbstrtbl-based is ", @sprintf("%.1f%%", (1-ratio4)*100), " faster")
elseif ratio4 < 1.2
    println("  ≈ Comparable performance (< 20% difference)")
else
    println("  ✗ pbstrtbl-based is ", @sprintf("%.1f%%", (ratio4-1)*100), " slower")
end
println()

# ============================================================================
# Summary
# ============================================================================
println("="^80)
println("SUMMARY")
println("="^80)
println()
println("1. Raw extraction (allpars vs pbstrtbl):")
println("   pbstrtbl / allpars = ", @sprintf("%.2fx", ratio1))
println()
println("2. tidyβ (allpars-based vs direct):")
println("   allpars-based / direct = ", @sprintf("%.2fx", ratio2))
println()
println("3. tidyσs (allpars-based vs direct):")
println("   allpars-based / direct = ", @sprintf("%.2fx", ratio3))
println()
println("4. make_estimate (pbstrtbl vs separate calls):")
println("   pbstrtbl / separate calls = ", @sprintf("%.2fx", ratio4))
println()

println("Interpretation:")
println("- Ratios < 1.0: First approach is faster")
println("- Ratios ≈ 1.0-1.2: Comparable performance")
println("- Ratios > 1.2: First approach is slower")
println()
println("="^80)
