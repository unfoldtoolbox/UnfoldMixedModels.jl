
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
    UnfoldMixedModels.tidyρs(m)

Correlations between the random effects, similar to `MixedModels.tidyσs`: one entry per fit & correlated random-effect pair, with fields `(iter, group, column1, column2, ρ)`.

The correlations are extracted by looping over `modelfit(m).fits`, rebuilding each λ block from the stored `θ`. Correlations that are structurally zero (e.g. `zerocorr`) are not included.
"""
function tidyρs(m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime})
    coll = modelfit(m)
    out = Vector{NamedTuple}()
    for (iter, f) in enumerate(coll.fits)
        offset = 0
        for (r, λtemplate) in enumerate(coll.λ)
            inds = coll.inds[r]
            k = size(λtemplate, 1)
            # rebuild the λ of this fit from θ, as MixedModels.setθ! does
            λdata = zeros(size(λtemplate))
            λdata[inds] .= f.θ[offset .+ (1:length(inds))]
            offset += length(inds)
            # free (non-zerocorr) λ entries, as MixedModels.indmat
            im = reshape([ix ∈ inds for ix = 1:abs2(k)], k, k)
            cnames = coll.fcnames[r]
            σs = [norm(view(λdata, i, 1:i)) for i = 1:k]
            for i = 2:k, j = 1:(i-1)
                iszero(dot(view(im, i, :), view(im, j, :))) && continue # zerocorr - no free shared entry
                ρ =
                    σs[i] * σs[j] == 0 ? NaN :
                    dot(view(λdata, i, :), view(λdata, j, :)) / (σs[i] * σs[j])
                push!(
                    out,
                    (
                        iter = iter,
                        group = keys(coll.fcnames)[r],
                        column1 = string(cnames[i]),
                        column2 = string(cnames[j]),
                        ρ = ρ,
                    ),
                )
            end
        end
    end
    return out
end


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
    corrs = tidyρs(m)
    ρs = [c.ρ for c in corrs]
    if ndims(coefs) == 3
        corr_est = permutedims(reshape(ρs, :, size(coefs, 2), size(coefs, 1)), [3, 2, 1])
        corr_group = permutedims(
            reshape([c.group for c in corrs], :, size(coefs, 2), size(coefs, 1)),
            [3, 2, 1],
        )
    else
        corr_est = permutedims(reshape(ρs, :, size(coefs, 1)), [2, 1])
        corr_group =
            permutedims(reshape([c.group for c in corrs], :, size(coefs, 1)), [2, 1])
    end

    estimate = cat(coefs, ranef(m), corr_est, dims = ndims(coefs))
    ranef_group = [x.group for x in MixedModels.tidyσs(m)]

    if ndims(coefs) == 3
        group_f =
            repeat([nothing], size(coefs, 1), size(coefs, 2), size(coefs, ndims(coefs)))


        # reshape to pred x time x chan and then invert to chan x time x pred
        ranef_group =
            permutedims(reshape(ranef_group, :, size(coefs, 2), size(coefs, 1)), [3 2 1])



        stderror_fixef = Unfold.stderror(m)
        stderror_ranef = fill(nothing, size(ranef(m)))
        stderror_corr = fill(nothing, size(corr_est))
        stderror = cat(stderror_fixef, stderror_ranef, stderror_corr, dims = 3)
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
    group = cat(group_f, ranef_group, corr_group, dims = ndims(coefs)) |> Unfold.poolArray
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


"""
correlation coefficient names for the coeftable, one per correlated random-effect pair, in the same order as `tidyρs`.
"""
function _corr_coefnames(uf)
    out = String[]
    for c in tidyρs(uf)
        c.iter == 1 || continue
        if uf isa UnfoldLinearMixedModelContinuousTime
            # coeftable splits names on ' : ', so eventname & colname components are needed
            event = split(c.column1, " : ")[1]
            push!(
                out,
                "$event : ρ: $(_compact_cname(c.column1)) - $(_compact_cname(c.column2)) : NaN",
            )
        else
            push!(out, "ρ: $(c.column1) - $(c.column2)")
        end
    end
    return out
end

# "Any : (Intercept) : -0.2" => "(Intercept)@-0.2"
_compact_cname(cname) = (p = split(cname, " : "); "$(p[2])@$(p[3])")

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
    return vcat(fe_coefnames, re_coefnames, _corr_coefnames(uf))
end

function Unfold.get_coefnames(uf::UnfoldLinearMixedModel)
    return vcat(vcat(Unfold.coefnames(formulas(uf))...), _corr_coefnames(uf))
end


Unfold.modelfit(uf::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime}) =
    uf.modelfit.collection
