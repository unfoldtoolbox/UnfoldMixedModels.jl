
"""
    _allpars_to_tidy_β(allpars)

Convert MixedModels.allpars output to tidyβ format.
Returns a vector of NamedTuples with fields (:iter, :coefname, :β).
"""
function _allpars_to_tidy_β(allpars)
    # Filter for β parameters
    mask = allpars.type .== "β"
    iters = allpars.iter[mask]
    names = allpars.names[mask]
    values = allpars.value[mask]
    
    # Convert to NamedTuple format like tidyβ
    colnms = (:iter, :coefname, :β)
    T = eltype(values)
    return [NamedTuple{colnms,Tuple{Int,Symbol,T}}((i, Symbol(nm), v)) 
            for (i, nm, v) in zip(iters, names, values)]
end

"""
    _allpars_to_tidy_σ(allpars)

Convert MixedModels.allpars output to tidyσs format.
Returns a vector of NamedTuples with fields (:iter, :group, :column, :σ).
"""
function _allpars_to_tidy_σ(allpars)
    # Filter for σ parameters (excluding residual)
    mask = (allpars.type .== "σ") .& .!ismissing.(allpars.group) .& (allpars.group .!= "residual")
    iters = allpars.iter[mask]
    groups = allpars.group[mask]
    columns = allpars.names[mask]
    values = allpars.value[mask]
    
    # Convert to NamedTuple format like tidyσs
    colnms = (:iter, :group, :column, :σ)
    T = eltype(values)
    return [NamedTuple{colnms,Tuple{Int,Symbol,Symbol,T}}((i, Symbol(g), Symbol(c), v)) 
            for (i, g, c, v) in zip(iters, groups, columns, values)]
end

"""
    _allpars_to_tidy_ρ(allpars)

Convert MixedModels.allpars output to tidyρs format.
Returns a vector of NamedTuples with fields (:iter, :group, :column1, :column2, :ρ).

Note: This function parses correlation parameter names in the format "col1, col2" 
as produced by MixedModels.allpars(). If the format changes upstream, this will need updating.
"""
function _allpars_to_tidy_ρ(allpars)
    # Filter for ρ parameters
    mask = allpars.type .== "ρ"
    iters = allpars.iter[mask]
    groups = allpars.group[mask]
    names = allpars.names[mask]
    values = allpars.value[mask]
    
    # Convert to NamedTuple format like tidyρs
    # names are in format "col1, col2" as produced by MixedModels.allpars
    colnms = (:iter, :group, :column1, :column2, :ρ)
    T = eltype(values)
    result = NamedTuple{colnms,Tuple{Int,Symbol,Symbol,Symbol,T}}[]
    for (i, g, n, v) in zip(iters, groups, names, values)
        # Parse "col1, col2" format - note the comma-space delimiter
        cols = split(n, ", ")
        if length(cols) == 2
            push!(result, NamedTuple{colnms,Tuple{Int,Symbol,Symbol,Symbol,T}}((i, Symbol(g), Symbol(cols[1]), Symbol(cols[2]), v)))
        else
            # This should not happen with standard MixedModels.allpars output
            # If it does, it indicates a format change or data issue
            @warn "Unexpected correlation parameter name format: \"$n\" (expected \"col1, col2\")"
        end
    end
    return result
end

function MixedModels.tidyσs(
    m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)
    # Use allpars to get all parameters, then filter for σ
    allpars_result = MixedModels.allpars(modelfit(m))
    t = _allpars_to_tidy_σ(allpars_result)
    reorder_tidyσs(t, Unfold.formulas(m))
end

function MixedModels.tidyβ(m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime})
    # Use allpars to get all parameters, then filter for β
    allpars_result = MixedModels.allpars(modelfit(m))
    _allpars_to_tidy_β(allpars_result)
end

"""
    tidyρs(m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime})

Extract correlations from random effects in the mixed model.
Returns a vector of NamedTuples with fields :iter, :group, :column1, :column2, :ρ
"""
function tidyρs(
    m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)
    # Use allpars to get all parameters, then filter for ρ
    allpars_result = MixedModels.allpars(modelfit(m))
    t = _allpars_to_tidy_ρ(allpars_result)
    reorder_tidyρs(t, Unfold.formulas(m))
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
        @debug "Could not fully reorder correlations due to formula/tidy mismatch, using original order"
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
    Unfold.make_estimate(m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime})

Extract all parameters (betas, sigmas, and correlations) from mixed model using MixedModels.pbstrtbl.
Returns estimates as channel x time x parameter arrays.

This function uses MixedModels.pbstrtbl() which provides all parameters in a wide table format,
avoiding the need to separately iterate through fits multiple times.
"""
function Unfold.make_estimate(
    m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)
    # Get the wide-format table with all parameters using pbstrtbl
    tbl = MixedModels.pbstrtbl(modelfit(m))
    
    # Get column names and identify parameter types
    colnames_tbl = Tables.columnnames(tbl)
    
    # Separate fixed effects (β), random effects (σ), and correlations (ρ)
    # Exclude :obj, :σ (residual), and :θ columns
    β_cols = filter(s -> startswith(string(s), "β"), colnames_tbl)
    σ_cols = filter(s -> startswith(string(s), "σ") && string(s) != "σ", colnames_tbl)
    ρ_cols = filter(s -> startswith(string(s), "ρ"), colnames_tbl)
    
    # Extract values as matrices
    β_vals = hcat([Tables.getcolumn(tbl, col) for col in β_cols]...)
    σ_vals = isempty(σ_cols) ? zeros(Float64, length(tbl), 0) : hcat([Tables.getcolumn(tbl, col) for col in σ_cols]...)
    ρ_vals = isempty(ρ_cols) ? zeros(Float64, length(tbl), 0) : hcat([Tables.getcolumn(tbl, col) for col in ρ_cols]...)
    
    # Combine all parameters
    all_vals = hcat(β_vals, σ_vals, ρ_vals)
    
    # Get grouping information - we still need tidyσs and tidyρs for this
    # as pbstrtbl doesn't include grouping metadata
    n_fixef = length(β_cols)
    n_ranef = length(σ_cols)
    n_corr = length(ρ_cols)
    
    # Reshape based on model type
    if m isa UnfoldLinearMixedModel
        # Mass univariate: channel x time x parameter
        ntime = length(Unfold.times(m)[1])
        nchan = modelfit(m).fits[end].channel
        # Transform from flat table (iterations x parameters) to (parameters x time x channel)
        # then permute to (channel x time x parameters)
        estimate = permutedims(reshape(all_vals', size(all_vals, 2), ntime, nchan), [3, 2, 1])
        
        # Create group labels
        group_f = repeat([nothing], nchan, ntime, n_fixef)
        
        # Get group names from model structure
        if n_ranef > 0
            ranef_group = [x.group for x in MixedModels.tidyσs(m)]
            ranef_group = permutedims(reshape(ranef_group, :, ntime, nchan), [3, 2, 1])
        else
            ranef_group = Array{Union{Nothing,Symbol}}(nothing, nchan, ntime, 0)
        end
        
        if n_corr > 0
            corr_group_tidy = tidyρs(m)
            corr_group = [x.group for x in corr_group_tidy]
            corr_group = permutedims(reshape(corr_group, :, ntime, nchan), [3, 2, 1])
        else
            corr_group = Array{Union{Nothing,Symbol}}(nothing, nchan, ntime, 0)
        end
        
        # Standard errors - only for fixed effects
        stderror_fixef = Unfold.stderror(m)
        stderror_ranef = fill(nothing, nchan, ntime, n_ranef)
        stderror_corr = fill(nothing, nchan, ntime, n_corr)
        stderror = cat(stderror_fixef, stderror_ranef, stderror_corr, dims=3)
    else
        # Continuous time: channel x parameter
        nchan = modelfit(m).fits[end].channel
        # Transform from flat table to (parameter x channel), then transpose to (channel x parameter)
        estimate = reshape(all_vals', size(all_vals, 2), nchan)'
        
        # Create group labels
        group_f = repeat([nothing], nchan, n_fixef)
        
        if n_ranef > 0
            ranef_group = [x.group for x in MixedModels.tidyσs(m)]
            ranef_group = reshape(ranef_group, :, nchan)'
        else
            ranef_group = Array{Union{Nothing,Symbol}}(nothing, nchan, 0)
        end
        
        if n_corr > 0
            corr_group_tidy = tidyρs(m)
            corr_group = [x.group for x in corr_group_tidy]
            corr_group = reshape(corr_group, :, nchan)'
        else
            corr_group = Array{Union{Nothing,Symbol}}(nothing, nchan, 0)
        end
        
        stderror = fill(nothing, size(estimate))
    end
    
    group = cat(group_f, ranef_group, corr_group, dims=ndims(estimate)) |> Unfold.poolArray
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
    # Get coefficient names from pbstrtbl column names
    tbl = MixedModels.pbstrtbl(modelfit(uf))
    colnames_tbl = Tables.columnnames(tbl)
    
    # Extract parameter names, excluding :obj, :σ (residual), and :θ
    β_cols = filter(s -> startswith(string(s), "β"), colnames_tbl)
    σ_cols = filter(s -> startswith(string(s), "σ") && string(s) != "σ", colnames_tbl)
    ρ_cols = filter(s -> startswith(string(s), "ρ"), colnames_tbl)
    
    # Get actual parameter names for display
    # For β: use the formula coefficient names
    formulas = Unfold.formulas(uf)
    if !isa(formulas, AbstractArray)
        formulas = [formulas]
    end
    fe_coefnames = vcat([coefnames(f.rhs[1]) for f in formulas]...)
    
    # For σ: use tidyσs to get the names with groups
    re_coefnames = isempty(σ_cols) ? String[] : vcat([coefnames(f.rhs[2:end]) for f in formulas]...)
    
    # For ρ: use get_corrnames
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
