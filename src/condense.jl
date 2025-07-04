
function MixedModels.tidyσs(
    m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)
    #    using MixedModels: AbstractReTerm
    t = MixedModels.tidyσs(modelfit(m))
    reorder_tidyσs(t, Unfold.formulas(m))
end

MixedModels.tidyβ(m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime}) =
    MixedModels.tidyβ(modelfit(m))


"""
helper function because `coefnames` returns an array only if number of coefs is larger than 1
"""
function _coefnames(t)
    c = coefnames(t)
    return isa(c, Vector) ? c : [c]
end



"""
    random_effect_groupings(t::MixedModels.AbstractReTerm)
Returns the random effect grouping term (rhs), similar to coefnames, which returns the left hand sides
"""
random_effect_groupings(t::AbstractTerm) = repeat([nothing], length(_coefnames(t.terms)))
random_effect_groupings(t::Unfold.TimeExpandedTerm) =
    repeat(random_effect_groupings(t.term), length(Unfold.colnames(t.basisfunction)))
random_effect_groupings(t::MixedModels.AbstractReTerm) =
    repeat([t.rhs.sym], length(_coefnames(t.lhs)))

random_effect_groupings(f::FormulaTerm) = vcat(random_effect_groupings.(f.rhs)...)
random_effect_groupings(t::Vector) = vcat(random_effect_groupings.(t)...)

"""
    reorder_tidyσs(t, f)
This function reorders a MixedModels.tidyσs output, according to the formula and not according to the largest RandomGrouping.

"""
function reorder_tidyσs(t, f)
    #@debug typeof(f)
    # get the order from the formula, this is the target
    f_order = random_effect_groupings(f) # formula order
    @debug f_order
    f_order = vcat(f_order...)
    @debug f_order

    # find the fixefs
    fixef_ix = isnothing.(f_order)



    f_order = string.(f_order[.!fixef_ix])
    @debug fixef_ix
    @debug coefnames(f)

    f_name = vcat(coefnames(f)...)[.!fixef_ix]

    # get order from tidy object
    t_order = [string(i.group) for i in t if i.iter == 1]
    t_name = [string(i.column) for i in t if i.iter == 1]

    # combine for formula and tidy output the group + the coefname
    @debug "f" f_order f_name
    @debug "t" t_order t_name
    f_comb = f_order .* f_name
    t_comb = t_order .* t_name

    # find for each formula output, the fitting tidy permutation
    reorder_ix = Int[]
    for f_coef in f_comb
        ix = findall(t_comb .== f_coef)
        # @debug t_comb, f_coef
        @assert length(ix) == 1 "error in reordering of MixedModels - please file a bugreport!"
        push!(reorder_ix, ix[1])
    end
    @assert length(reorder_ix) == length(t_comb)
    #@debug reorder_ix
    # repeat and build the index for all timepoints
    reorder_ix_all = repeat(reorder_ix, length(t) ÷ length(reorder_ix))
    for k = 1:length(reorder_ix):length(t)
        reorder_ix_all[k:(k+length(reorder_ix)-1)] .+= (k - 1)
    end

    return t[reorder_ix_all]


end

"""
    Unfold.make_estimate(m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)
extracts betas (and sigma's for mixed models) with string grouping indicator

returns as a ch x beta, or ch x time x beta (for mass univariate)
"""
function Unfold.make_estimate(
    m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)

    coefs = coef(m)
    estimate = cat(coefs, ranef(m), dims = ndims(coefs))
    ranef_group = [x.group for x in MixedModels.tidyσs(m)]

    if ndims(coefs) == 3
        group_f =
            repeat([nothing], size(coefs, 1), size(coefs, 2), size(coefs, ndims(coefs)))


        # reshape to pred x time x chan and then invert to chan x time x pred
        ranef_group =
            permutedims(reshape(ranef_group, :, size(coefs, 2), size(coefs, 1)), [3 2 1])



        stderror_fixef = Unfold.stderror(m)
        stderror_ranef = fill(nothing, size(ranef(m)))
        stderror = cat(stderror_fixef, stderror_ranef, dims = 3)
    else
        group_f = repeat([nothing], size(coefs, 1), size(coefs, 2))

        # reshape to time x channel
        ranef_group = reshape(ranef_group, :, size(coefs, 1))
        # permute to channel x time
        ranef_group = permutedims(ranef_group, [2, 1])

        #@debug size(ranef_group)
        #
        #ranef_group = repeat(["ranef"], size(coefs, 1), size(ranef(m), 2))
        #@debug size(ranef_group)
        stderror = fill(nothing, size(estimate))
    end
    group = cat(group_f, ranef_group, dims = ndims(coefs)) |> Unfold.poolArray
    return Float64.(estimate), stderror, group
end

function Unfold.stderror(
    m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)
    return permutedims(
        reshape(vcat([[b.se...] for b in modelfit(m).fits]...), reverse(size(coef(m)))),
        [3, 2, 1],
    )
end


function Unfold.get_coefnames(uf::UnfoldLinearMixedModelContinuousTime)
    # special case here, because we have to reorder the random effects to the end, else labels get messed up as we concat (coefs,ranefs)
    #   coefnames = Unfold.coefnames(formula(uf))
    #    coefnames(formula(uf)[1].rhs[1])
    formulas = Unfold.formulas(uf)
    if !isa(formulas, AbstractArray) # in case we have only a single basisfunction
        formulas = [formulas]
    end
    fe_coefnames = vcat([coefnames(f.rhs[1]) for f in formulas]...)
    re_coefnames = vcat([coefnames(f.rhs[2:end]) for f in formulas]...)
    return vcat(fe_coefnames, re_coefnames)
end


Unfold.modelfit(uf::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime}) =
    uf.modelfit.collection
