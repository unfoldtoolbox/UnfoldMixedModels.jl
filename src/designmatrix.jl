
function StatsModels.coefnames(term::RandomEffectsTerm)
    coefnames(term.lhs)
end

check_groupsorting(r::MatrixTerm) = check_groupsorting(r.terms)
function check_groupsorting(r::Tuple)
    @debug "checking group sorting"
    ix = findall(isa.(r, MixedModels.AbstractReTerm))

    rhs(x::RandomEffectsTerm) = x.rhs
    rhs(x::MixedModels.ZeroCorr) = rhs(x.term)
    groupvars = [map(x -> rhs(x).sym, r[ix])...]


    @assert groupvars == sort(groupvars) "random effects have to be alphabetically ordered. e.g. (1+a|X) + (1+a|A) is not allowed. Please reorder"
end
function Unfold.unfold_apply_schema(
    type::Type{<:Union{<:UnfoldLinearMixedModel,<:UnfoldLinearMixedModelContinuousTime}},
    f,
    schema,
)
    @debug "LMM apply schema"
    f_new = apply_schema(f, schema, MixedModels.LinearMixedModel)
    check_groupsorting(f_new.rhs)
    return f_new

end



function StatsModels.coefnames(term::MixedModels.ZeroCorr)
    coefnames(term.term)
end

function lmm_combine_modelmatrices!(Xcomb, X1, X2)
    # we have random effects
    # combine REMats in single-eventtpe formulas ala y ~ (1|x) + (a|x)
    modelmatrix1 = MixedModels._amalgamate([X1.modelmatrix[2:end]...], Float64)
    modelmatrix2 = MixedModels._amalgamate([X2.modelmatrix[2:end]...], Float64)

    Xcomb = (Xcomb, modelmatrix1..., modelmatrix2...)

    # Next we make the ranefs all equal size
    equalize_ReMat_lengths!(Xcomb[2:end])

    # check if ranefs can be amalgamated. If this fails, then MixedModels tried to amalgamate over different eventtypes and we should throw the warning
    # if it success, we have to check if the size before and after is identical. If it is not, it tried to amalgamize over different eventtypes which were of the same length

    try
        reterms = MixedModels._amalgamate([Xcomb[2:end]...], Float64)
        if length(reterms) != length(Xcomb[2:end])
            throw("IncompatibleRandomGroupings")
        end
    catch e
        @error "Error, you seem to have two different eventtypes with the same random-effect grouping variable. \n
        This is not allowed, you have to rename one. Example:\n
        eventA: y~1+(1|item) \n
        eventB: y~1+(1|item)  \n
        This leads to this error. Rename the later one\n
        eventB: y~1+(1|itemB) "
        throw("IncompatibleRandomGroupings")
    end

    return Xcomb
end


function change_ReMat_size!(remat::MixedModels.AbstractReMat, m::Integer)

    n = m - length(remat.refs) # missing elements
    if n < 0
        deleteat!(remat.refs, range(m + 1, stop = length(remat.refs)))
        # MixedModels requires length(refs) == size(z,2) (and wtz, adjA respectively)
        remat.adjA = remat.adjA[:, 1:m]
        remat.wtz = remat.wtz[:, 1:m]
        remat.z = remat.z[:, 1:m]
    elseif n > 0
        append!(remat.refs, repeat([remat.refs[end]], n))
        remat.adjA = [remat.adjA sparse(repeat([0], size(remat.adjA)[1], n))]
        remat.wtz = [remat.wtz zeros((size(remat.wtz)[1], n))]
        remat.z = [remat.z zeros((size(remat.z)[1], n))]
    end

    #ReMat{T,S}(remat.trm, refs, levels, cnames, z, wtz, λ, inds, adjA, scratch)
end


function equalize_ReMat_lengths!(remats::NTuple{A,MixedModels.AbstractReMat}) where {A}
    # find max length
    m = maximum([x[1] for x in size.(remats)])
    @debug "combining lengths: $m"
    # for each reMat
    for k in range(1, length = length(remats))
        remat = remats[k]
        if size(remat)[1] == m
            continue
        end
        # prolong if necessary
        change_ReMat_size!(remat, m)

    end
end


mutable struct SparseReMat{T,S} <: MixedModels.AbstractReMat{T}
    trm::Any
    refs::Vector{Int32}
    levels::Any
    cnames::Vector{String}
    z::SparseMatrixCSC{T}
    wtz::SparseMatrixCSC{T}
    λ::LowerTriangular{T,Matrix{T}}
    inds::Vector{Int}
    adjA::SparseMatrixCSC{T,Int32}
    scratch::Matrix{T}
end

"""
$(SIGNATURES)
This function timeexpands the random effects and generates a ReMat object
"""
function StatsModels.modelcols(
    term::Unfold.TimeExpandedTerm{
        <:Union{<:RandomEffectsTerm,<:MixedModels.AbstractReTerm},
    },
    tbl,
)
    # exchange this to get ZeroCorr to work


    tbl = DataFrame(tbl)
    # get the non-timeexpanded reMat
    reMat = modelcols(term.term, tbl)

    # Timeexpand the designmatrix
    z = transpose(Unfold.time_expand(transpose(reMat.z), term, tbl))

    z = disallowmissing(z) # can't have missing in here


    if hasfield(typeof(term.term), :rhs)
        rhs = term.term.rhs


    elseif hasfield(typeof(term.term.term), :rhs)
        # we probably have something like zerocorr, which does not need to show a .rhs necessarily
        rhs = term.term.term.rhs
    else

        println("term.term: $(dump(term.term))")
        error("unknown RE structure, has no field .rhs:$(typeof(term.term))")

    end

    time = tbl[!, term.eventfields[1]]

    # Determine which grouping level covers each time sample.
    # Every event contributes to the time samples spanned by its basisfunction, and all
    # contributions to a sample have to belong to the same grouping level, otherwise the
    # random effects of different levels would leak into each other.
    # reMat.refs contains the level index (into reMat.levels) of each event.
    levels = reMat.levels
    refs = zeros(Int32, size(z, 2)) # grouping level per time sample, 0 = not assigned yet
    shift = Unfold.shift_onset(term.basisfunction)
    w = Unfold.width(term.basisfunction)
    for (onset, lev) in zip(time, reMat.refs)
        s = floor(Int, onset) + shift
        for t in max(s, 1):min(s + w - 1, length(refs))
            if refs[t] == 0
                refs[t] = lev
            elseif refs[t] != lev
                error(
                    "overlap in random effects structure detected: time sample $t is covered by " *
                    "the basisfunctions of more than one level of grouping variable $(rhs.sym) " *
                    "($(levels[refs[t]]) and $(levels[lev])). The random effects of different " *
                    "levels would leak into each other. Increase the distance between events of " *
                    "different levels or shorten the basisfunction.",
                )
            end
        end
    end

    # Time samples not covered by any event have z == 0, so the surrounding level can be extended into the gaps
    prev = one(Int32)
    for t in eachindex(refs)
        if refs[t] == 0
            refs[t] = prev
        else
            prev = refs[t]
        end
    end

    # Other variables with implementaions taken from the LinerMixedModel function
    wtz = z
    trm = term

    S = size(z, 1)
    T = eltype(z)
    λ = LowerTriangular(Matrix{T}(I, S, S))

    inds = MixedModels.sizehint!(Int[], (S * (S + 1)) >> 1) # double trouble? can this line go?
    m = reshape(1:abs2(S), (S, S))
    inds = sizehint!(Int[], (S * (S + 1)) >> 1)
    for j = 1:S
        for i = j:S
            # We currently restrict to diagonal entries
            # Once mixedmodels#293 is pushed, we can relax this and use zerocorr()
            if !(typeof(term.term) <: MixedModels.ZeroCorr) || (i == j) # for diagonal
                push!(inds, m[i, j])
            end
        end
    end

    cnames = coefnames(term)
    #print(refs)
    adjA = MixedModels.adjA(refs, z)
    scratch = Matrix{T}(undef, (S, length(levels)))

    # Once MixedModels.jl supports it, can be replaced with: SparseReMat
    ReMat{T,S}(rhs, refs, levels, cnames, z, wtz, λ, inds, adjA, scratch)
end
function change_modelmatrix_size!(m, fe::AbstractSparseMatrix, remats)
    change_ReMat_size!.(remats, Ref(m))
    @debug "changemodelmatrix" typeof(fe)
    fe = SparseMatrixCSC(m, fe.n, fe.colptr, fe.rowval, fe.nzval)
    return (fe, remats...)
end
function change_modelmatrix_size!(m, fe::AbstractMatrix, remats)
    change_ReMat_size!.(remats, Ref(m))
    fe = fe[1:m, :]
    return (fe, remats...)
end

"""
    modelmatrices(modelmatrix::Tuple)
in the case of a Tuple (MixedModels - FeMat/ReMat Tuple), returns only the FeMat part
"""
Unfold.modelmatrices(modelmatrix::Tuple) = modelmatrix[1]

#modelcols(rhs::MatrixTerm, tbl) = modelcols.(rhs, Ref(tbl))
