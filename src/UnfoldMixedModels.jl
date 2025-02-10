module UnfoldMixedModels

using Unfold: BasisFunction
using Unfold

import Unfold: make_estimate, modelmatrices, typify
using MixedModels
#import MixedModels.FeMat # extended for sparse femats, type piracy => issue on MixedModels.jl github
using StaticArrays # for MixedModels extraction of parametrs (inherited from MixedModels.jl, not strictly needed )
import MixedModels: likelihoodratiotest, ranef
using StatsModels
import StatsModels: fit!, coef, coefnames, modelcols, modelmatrix
import StatsAPI: pvalue
using SparseArrays
using DocStringExtensions
using LinearAlgebra # LowerTriangular
using DataFrames
using ProgressMeter
using SimpleTraits
using Reexport
using StatsAPI

include("typedefinitions.jl")
include("basisfunctions.jl")
include("condense.jl")
include("designmatrix.jl")
include("fit.jl")
include("statistics.jl")
include("timeexpandedterm.jl")
include("effects.jl")

@reexport using Unfold
@reexport using MixedModels
export DesignMatrixLinearMixedModel, DesignMatrixLinearMixedModelContinuousTime
export LinearMixedModelFitCollection
export UnfoldLinearMixedModel, UnfoldLinearMixedModelContinuousTime
export likelihoodratiotest, pvalue

end
