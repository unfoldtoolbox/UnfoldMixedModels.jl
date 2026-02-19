
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
    tidyρs(m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime})

Extract correlations from random effects in the mixed model.
Returns a vector of NamedTuples with fields :iter, :group, :column1, :column2, :ρ
"""
function tidyρs(
    m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)
    t = tidyρs(modelfit(m))
    reorder_tidyρs(t, Unfold.formulas(m))
end

"""
    tidyρs(bsamp::LinearMixedModelFitCollection)

Extract correlations from a LinearMixedModelFitCollection.
Returns a vector of NamedTuples with fields :iter, :group, :column1, :column2, :ρ
"""
function tidyρs(bsamp::LinearMixedModelFitCollection{T}) where {T}
    fits = bsamp.fits
    fcnames = bsamp.fcnames
    λ = bsamp.λ
    colnms = (:iter, :group, :column1, :column2, :ρ)
    result = sizehint!(
        NamedTuple{colnms,Tuple{Int,Symbol,Symbol,Symbol,T}}[], 
        length(fits) * sum(λ_mat -> (size(λ_mat, 1) * (size(λ_mat, 1) - 1)) ÷ 2, λ)
    )
    
    for (iter, r) in enumerate(fits)
        MixedModels.setθ!(bsamp, iter)    # install r.θ in λ
        σ = coalesce(r.σ, one(T))
        for (grp, ll) in zip(keys(fcnames), λ)
            cnames = getproperty(fcnames, grp)
            # Only process if we have multiple random effects (and thus correlations)
            if length(cnames) > 1
                # Extract correlations from the lower triangular matrix
                # Following the logic from MixedModels.jl σρ! function
                dat = ll.data
                k = size(dat, 1)
                # Normalize rows and compute correlations
                normalized = copy(dat)
                for i in 1:k
                    len = sqrt(sum(abs2(dat[i, j]) for j in 1:i))
                    if len > 0
                        for j in 1:i
                            normalized[i, j] = dat[i, j] / len
                        end
                    end
                end
                # Compute correlations as dot products of normalized rows
                for i in 2:k
                    for j in 1:(i-1)
                        corr = sum(normalized[i, k] * normalized[j, k] for k in 1:i)
                        push!(result, NamedTuple{colnms}((iter, grp, Symbol(cnames[i]), Symbol(cnames[j]), corr)))
                    end
                end
            end
        end
    end
    return result
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
    reorder_tidyρs(t, f)
This function reorders a tidyρs output, according to the formula and not according to the largest RandomGrouping.

"""
function reorder_tidyρs(t, f)
    # If there are no correlations, return empty
    if isempty(t)
        return t
    end
    
    # get the order from the formula, this is the target
    f_order = random_effect_groupings(f) # formula order
    f_order = vcat(f_order...)
    
    # find the fixefs
    fixef_ix = isnothing.(f_order)
    f_order = string.(f_order[.!fixef_ix])
    f_name = vcat(coefnames(f)...)[.!fixef_ix]
    
    # get order from tidy object
    t_order = [string(i.group) for i in t if i.iter == 1]
    t_name1 = [string(i.column1) for i in t if i.iter == 1]
    t_name2 = [string(i.column2) for i in t if i.iter == 1]
    
    # combine for formula and tidy output the group + the coefnames
    f_comb = f_order .* f_name
    t_comb = t_order .* t_name1 .* t_name2
    
    # For correlations, we need to match pairs
    # Build a mapping based on group and column names
    reorder_ix = Int[]
    for (i, (grp_f, name_f)) in enumerate(zip(f_order, f_name))
        # Find all correlations involving this group and name
        for (j, row) in enumerate(t)
            if row.iter == 1 && string(row.group) == grp_f && 
               (string(row.column1) == name_f || string(row.column2) == name_f)
                # This correlation involves our target
                # Check if we haven't already added it
                if !(j in reorder_ix)
                    push!(reorder_ix, j)
                end
            end
        end
    end
    
    # If we found a different number of correlations, just return original order
    # This handles the complex case where formula order doesn't match
    if length(reorder_ix) != length(t_comb)
        @debug "Warning: Could not fully reorder correlations, using original order"
        return t
    end
    
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
extracts betas (and sigma's and correlations for mixed models) with string grouping indicator

returns as a ch x beta, or ch x time x beta (for mass univariate)
"""
function Unfold.make_estimate(
    m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)

    coefs = coef(m)
    ranef_sigma = ranef(m)
    corrs = ranefcorr(m)
    
    ranef_group = [x.group for x in MixedModels.tidyσs(m)]
    corr_group_tidy = tidyρs(m)
    corr_group = [x.group for x in corr_group_tidy]

    if ndims(coefs) == 3
        group_f =
            repeat([nothing], size(coefs, 1), size(coefs, 2), size(coefs, ndims(coefs)))


        # reshape to pred x time x chan and then invert to chan x time x pred
        ranef_group =
            permutedims(reshape(ranef_group, :, size(coefs, 2), size(coefs, 1)), [3 2 1])

        # reshape correlations similarly
        if !isempty(corr_group) && size(corrs, 3) > 0
            corr_group =
                permutedims(reshape(corr_group, :, size(coefs, 2), size(coefs, 1)), [3 2 1])
        else
            corr_group = Array{Union{Nothing,Symbol}}(nothing, size(coefs, 1), size(coefs, 2), 0)
        end

        stderror_fixef = Unfold.stderror(m)
        stderror_ranef = fill(nothing, size(ranef_sigma))
        
        # Only concatenate if we have correlations
        if size(corrs, 3) > 0
            estimate = cat(coefs, ranef_sigma, corrs, dims = ndims(coefs))
            stderror_corr = fill(nothing, size(corrs))
            stderror = cat(stderror_fixef, stderror_ranef, stderror_corr, dims = 3)
        else
            estimate = cat(coefs, ranef_sigma, dims = ndims(coefs))
            stderror = cat(stderror_fixef, stderror_ranef, dims = 3)
        end
    else
        group_f = repeat([nothing], size(coefs, 1), size(coefs, 2))

        # reshape to time x channel
        ranef_group = reshape(ranef_group, :, size(coefs, 1))
        # permute to channel x time
        ranef_group = permutedims(ranef_group, [2, 1])

        # reshape correlations similarly
        if !isempty(corr_group) && size(corrs, 2) > 0
            corr_group = reshape(corr_group, :, size(coefs, 1))
            corr_group = permutedims(corr_group, [2, 1])
        else
            corr_group = Array{Union{Nothing,Symbol}}(nothing, size(coefs, 1), 0)
        end

        # Only concatenate if we have correlations
        if size(corrs, 2) > 0
            estimate = cat(coefs, ranef_sigma, corrs, dims = ndims(coefs))
        else
            estimate = cat(coefs, ranef_sigma, dims = ndims(coefs))
        end
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
    corr_names = get_corrnames(uf)
    return vcat(fe_coefnames, re_coefnames, corr_names)
end

"""
    get_corrnames(uf::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime})

Get correlation parameter names from a mixed model.
Returns names in the format "ρ: col1, col2" for correlations between random effects.
"""
function get_corrnames(uf::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime})
    corrs = tidyρs(uf)
    if isempty(corrs)
        return String[]
    end
    # Get unique correlation names from the first iteration
    unique_corrs = filter(x -> x.iter == 1, corrs)
    return ["ρ: $(x.column1), $(x.column2)" for x in unique_corrs]
end


Unfold.modelfit(uf::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime}) =
    uf.modelfit.collection
