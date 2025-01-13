```@meta
CurrentModule = UnfoldMixedModels
```

# UnfoldMixedModels

Documentation for [UnfoldMixedModels](https://github.com/unfoldtoolbox/UnfoldMixedModels.jl).

Using this package, one can fit Linear Mixed Models in a mass-univariate way (for every time-point and channel); but also combined with overlap correction (experimental!)

As `UnfoldMixedModels.jl` is like an *addon* to `Unfold.jl`, we recommend checking out these tutorials first.

```julia
using UnfoldMixedModels
using UnfoldSim
data, evts = UnfoldSim.predef_eeg(10;return_epoched=true) # 10 subjects
data = reshape(data,size(data,1),:) # concatenate subjects

times = range(-0.1,0.5,size(data,1)) # arbitrary time-vector

fLMM = @formula 0 ~ 1 + condition + (1|subject) + (1|item)
fit(UnfoldModel, [Any=>(f, times)], evts, data)
nothing #hide
```
