# ![UnfoldMixedModels - MixedModels in EEG](https://github.com/user-attachments/assets/a2fd2d6b-4d9c-4d23-b3ca-936fe4055a99)



[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://unfoldtoolbox.github.io/UnfoldMixedModels.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://unfoldtoolbox.github.io/UnfoldMixedModels.jl/dev)
[![Build Status](https://github.com/unfoldtoolbox/UnfoldMixedModels.jl/workflows/Test/badge.svg)](https://github.com/unfoldtoolbox/UnfoldMixedModels.jl/actions)
[![Test workflow status](https://github.com/unfoldtoolbox/UnfoldMixedModels.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/unfoldtoolbox/UnfoldMixedModels.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/unfoldtoolbox/UnfoldMixedModels.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/unfoldtoolbox/UnfoldMixedModels.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/unfoldtoolbox/UnfoldMixedModels.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/unfoldtoolbox/UnfoldMixedModels.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/unfoldtoolbox/UnfoldMixedModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/unfoldtoolbox/UnfoldMixedModels.jl)
[![DOI](https://zenodo.org/badge/DOI/FIXME)](https://doi.org/FIXME)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
[![All Contributors](https://img.shields.io/github/all-contributors/unfoldtoolbox/UnfoldMixedModels.jl?labelColor=5e1ec7&color=c0ffee&style=flat-square)](#contributors)

|Estimation|Visualisation|Simulation|BIDS pipeline|Decoding|Statistics|MixedModelling|
|---|---|---|---|---|---|---|
| <a href="https://github.com/unfoldtoolbox/Unfold.jl/tree/main"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623787-757575d0-aeb9-4d94-a5f8-832f13dcd2dd.png" alt="Unfold.jl Logo"></a> | <a href="https://github.com/unfoldtoolbox/UnfoldMakie.jl"><img  src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623793-37af35a0-c99c-4374-827b-40fc37de7c2b.png" alt="UnfoldMakie.jl Logo"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldSim.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623795-328a4ccd-8860-4b13-9fb6-64d3df9e2091.png" alt="UnfoldSim.jl Logo"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldBIDS.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277622460-2956ca20-9c48-4066-9e50-c5d25c50f0d1.png" alt="UnfoldBIDS.jl Logo"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldDecode.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277622487-802002c0-a1f2-4236-9123-562684d39dcf.png" alt="UnfoldDecode.jl Logo"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldStats.jl"><img  src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623799-4c8f2b5a-ea84-4ee3-82f9-01ef05b4f4c6.png" alt="UnfoldStats.jl Logo"></a>|<a href=""><img src="https://github.com/user-attachments/assets/ffb2bba6-3a30-48b7-9849-7d4e7195b297" alt="UnfoldMixedModels.jl logo"></a>|

UnfoldMixedModels.jl is a package to perform hierarchical regression / **linear mixed models** on biological signals. As an experimental feature, it further allows to perform simultaneous overlap-correction / deconvolution.

This kind of modelling is also known as encoding modeling, linear deconvolution, Temporal Response Functions (TRFs), linear system identification, and probably under other names. fMRI models with HRF-basis functions and pupil-dilation bases are also supported.

## Getting started

### 🐍Python User?

We clearly recommend Julia 😉 - but [Python users can use juliacall/Unfold directly from python!](https://unfoldtoolbox.github.io/Unfold.jl/dev/generated/HowTo/juliacall_unfold/)

### Julia installation

<details>
<summary>Click to expand</summary>

The recommended way to install julia is [juliaup](https://github.com/JuliaLang/juliaup).
It allows you to, e.g., easily update Julia at a later point, but also test out alpha/beta versions etc.

TL:DR; If you dont want to read the explicit instructions, just copy the following command

#### Windows

AppStore -> JuliaUp,  or `winget install julia -s msstore` in CMD

#### Mac & Linux

`curl -fsSL https://install.julialang.org | sh` in any shell
</details>

### UnfoldMixedModels.jl installation

```julia
using Pkg
Pkg.add("UnfoldMixedModels")
```

## Usage

Please check out [the documentation](https://unfoldtoolbox.github.io/UnfoldMixedModels.jl) for extensive tutorials, explanations and more!

### Tipp on Docs

You can read the docs online: [![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://unfoldtoolbox.github.io/UnfoldMixedModels.jl/stable)  - or use the `?fit`, `?effects` julia-REPL feature. To filter docs, use e.g. `?fit(::UnfoldMixedModel)`

Here is a quick overview on what to expect.

### What you need

```julia
using UnfoldMixedModels

events::DataFrame

# formula with or without random effects

fLMM = @formula 0~1+condA+(1|subject) + (1|item)

# in case of [overlap-correction] we need continuous data plus per-eventtype one basisfunction (typically firbasis)
data::Array{Float64,2}
basis = firbasis(τ=(-0.3,0.5),srate=250) # for "timeexpansion" / deconvolution

# in case of [mass univariate] we need to epoch the data into trials, and a accompanying time vector
epochs::Array{Float64,3} # channel x time x epochs (n-epochs == nrows(events))
times = range(0,length=size(epochs,3),step=1/sampling_rate)
```

To fit any of the models, Unfold.jl offers a unified syntax:

| Overlap-Correction | Mixed Modelling | julia syntax |
|:---:|:---:|---|
|  | x | `fit(UnfoldModel,[Any=>(fLMM,times)),evts,data_epoch]` |
| x | x | `fit(UnfoldModel,[Any=>(fLMM,basis)),evts,data]` |

</details>

## Contributions

Contributions are very welcome. These could be typos, bugreports, feature-requests, speed-optimization, new solvers, better code, better documentation.

### How-to Contribute

You are very welcome to raise issues and start pull requests!

### Adding Documentation

1. We recommend to write a Literate.jl document and place it in `docs/literate/FOLDER/FILENAME.jl` with `FOLDER` being `HowTo`, `Explanation`, `Tutorial` or `Reference` ([recommended reading on the 4 categories](https://documentation.divio.com/)).
2. Literate.jl converts the `.jl` file to a `.md` automatically and places it in `docs/src/generated/FOLDER/FILENAME.md`.
3. Edit [make.jl](https://github.com/unfoldtoolbox/UnfoldMixedModels.jl/blob/main/docs/make.jl) with a reference to `docs/src/generated/FOLDER/FILENAME.md`.

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://allcontributors.org/docs/en/specification) specification.

Contributions of any kind welcome!

## Citation

For now, please cite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5759066.svg)](https://doi.org/10.5281/zenodo.5759066) and/or [Ehinger & Dimigen](https://peerj.com/articles/7838/)

## Acknowledgements

This work was initially supported by the Center for Interdisciplinary Research, Bielefeld (ZiF) Cooperation Group "Statistical models for psychological and linguistic data".

Funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany´s Excellence Strategy – EXC 2075 – 390740016
