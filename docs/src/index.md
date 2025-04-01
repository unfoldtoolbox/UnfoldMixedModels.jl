```@meta
CurrentModule = UnfoldMixedModels
```

# UnfoldMixedModels.jl Documentation

Welcome to [UnfoldMixedModels.jl](https://github.com/unfoldtoolbox/UnfoldMixedModels.jl): a Julia package to analyse timeseries with Linear Mixed Models.
This is an standalone-addon to `Unfold.jl` with similar syntax, but optimized for LMMs / Hierarchical Models.

```@raw html
<div style="width:60%; margin: auto;">
</div>
```

## Key features

- **MixedModels.jl:** Powered by the up-to-100x faster implementation of MixedModels.jl
- **Subject and Item effects:** Fit the full spectrum of LMMs, random slopes and all!
- **Beta: clusterpermutation:** Combine with `UnfoldStats.jl` and fit LMM clusterpermutation tests
- **Alpha: Overlap:** Model overlap and LMMs (experimental!)

## Installation

```julia-repl
julia> using Pkg; Pkg.add("UnfoldMixedModels")
```

## Usage example

```julia
using UnfoldMixedModels
using UnfoldSim
data, evts = UnfoldSim.predef_eeg(10;return_epoched=true) # 10 subjects
data = reshape(data,size(data,1),:) # concatenate subjects

times = range(-0.1,0.5,size(data,1)) # arbitrary time-vector

fLMM = @formula 0 ~ 1 + condition + (1 + condition|subject) + (1|item)
fit(UnfoldModel, [Any=>(f, times)], evts, data)
nothing #hide
```

## Where to start: Learning roadmap

### 0. First first steps

📌 Goal: Familiarize yourself with rERP Unfold.jl fitting & MixedModels.jl /
🔗 [Unfold.jl Quickstart](https://unfoldtoolbox.github.io/Unfold.jl/stable/tutorials/lm_mu/) | [MixedModels.jl](https://juliastats.org/MixedModels.jl/dev/constructors/)

### 1. First steps

📌 Goal: Learn to fit a mass univariate Linear Mixed Model /
🔗 [Mass Univariate Linear Mixed Models](@ref lmm_massunivariate)

### 2. Intermediate topics

📌 Goal: Extract statistics and cluster permutation tests /
🔗 [How To get P-Values for Mass-Univariate LMM](@ref lmm_overlap) | [LMM Cluster Permutation tests](https://github.com/unfoldtoolbox/UnfoldStats.jl)

We further recommend to skim this online book (from the authors of MixedModels.jl): [embraceuncertaintybook.com/](https://embraceuncertaintybook.com/)

You should also learn about model simplification (keep it maximal ..?) and contrast codings. Enjoy!

## Statement of need

EEG researchers often analyse complex experimental procedures and want to generalize them to populations of subjects, items, schools etc. ([Generalization Crisis - Yarkoni 2020](https://doi.org/10.1017/S0140525X20001685 )). In case of hierarchical structures (e.g. repeated trials in subjects, different stimuli used) the Linear Mixed Model has become very popular. Unfortunately, fitting such models can be quite involved, especially for EEG data which require massive-modelfitting for each sensor and channel. MixedModels.jl provides a fast way for fitting, and UnfoldMixedModels.jl provides the bookkeeping to do so in a massive way. Users can easily extract fixed and random effects over time and sensors, do statistical testing and even correction for multiple comparisons (via UnfoldStats.jl / MixedModelsPermutations.jl).

```@raw html
<!---
Note: The statement of need is also used in the `README.md`. Make sure that they are synchronized.
-->
```
