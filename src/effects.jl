function StatsModels.modelmatrix(model::UnfoldLinearMixedModel, bool)
    @assert bool == false "time continuous model matrix is not implemented for a `UnfoldLinearMixedModel`"
    return modelmatrix(model)
end


# mixedModels case - just use the FixEff, ignore the ranefs
Unfold.typify(reference_grid, form, m::Tuple; typical) =
    Unfold.typify(reference_grid, form, m[1]; typical)
