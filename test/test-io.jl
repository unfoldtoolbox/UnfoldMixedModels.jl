
save_path = mktempdir(; cleanup = true)#tempdir()
#----
## 2. Test data set:
# - Generate a 2x2 design with Hanning window for multiple subjects (using UnfoldSim)
# - Use a Mixed-effects Unfold model

data, evts = UnfoldSim.predef_2x2(; n_subjects = 5, return_epoched = true);
data = reshape(data, size(data, 1), :)
subj_idx = [parse(Int, split(string(s), 'S')[2]) for s in evts.subject]
evts.latency .+= size(data, 1) .* (subj_idx .- 1)


# Define a model formula with interaction term and random effects for subjects
f2 = @formula(0 ~ 1 + A * B + (1 | subject));
τ2 = [-0.1, 1];
sfreq2 = 100;
times2 = range(τ2[1], length = size(data, 1), step = 1 ./ sfreq2);

m2 = Unfold.fit(
    UnfoldModel,
    Dict(Any => (f2, times2)),
    evts,
    reshape(data, 1, size(data)...),
);

save(joinpath(save_path, "m2_compressed2.jld2"), m2; compress = true)
m2_loaded =
    load(joinpath(save_path, "m2_compressed2.jld2"), UnfoldModel, generate_Xs = true)

@testset "2x2 MultiSubjectDesign Mixed-effects model" begin
    # save the model to a compressed .jld2 file and load it again
    save(joinpath(save_path, "m2_compressed2.jld2"), m2; compress = true)
    m2_loaded =
        load(joinpath(save_path, "m2_compressed2.jld2"), UnfoldModel, generate_Xs = true)


    @test isempty(Unfold.modelmatrices(designmatrix(m2_loaded))[1]) == false

    @test typeof(m2) == typeof(m2_loaded)
    @test coeftable(m2) == coeftable(m2_loaded)
    @test modelfit(m2).fits == modelfit(m2_loaded).fits
    @test Unfold.events(m2) == Unfold.events(m2_loaded)
    @test modelmatrix(m2) == modelmatrix(m2_loaded)

    # Test whether the effects function works with the loaded models
    # and the results match the ones of the original model
    conditions = Dict(:A => levels(evts.A), :B => levels(evts.B))

    # The effects function is currently not defined for UnfoldLinearMixedModel
    eff2 = effects(conditions, m2)
    eff2_loaded = effects(conditions, m2_loaded)

    @test eff2 == eff2_loaded

    # load the model without reconstructing the designmatrix
    m2_loaded_without_dm =
        load(joinpath(save_path, "m2_compressed2.jld2"), UnfoldModel, generate_Xs = false)

    @test isempty(modelmatrix(designmatrix(m2_loaded_without_dm))[1]) == true
end
