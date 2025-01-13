
@testset "effects MixedModel" begin
    data, evts = UnfoldSim.predef_eeg(10; return_epoched = true)
    data = reshape(data, size(data, 1), :)

    m = fit(
        UnfoldModel,
        @formula(0 ~ 1 + condition + (1 + condition | subject)),
        evts,
        data,
        1:size(data, 1),
    )
    eff = effects(Dict(:condition => ["car", "face"]), m)
end


@testset "effects MixedModelContinuousTime" begin
    data, evts = UnfoldSim.predef_eeg(10; sfreq = 10, return_epoched = false)
    data = data[:]
    subj_idx = [parse(Int, split(string(s), 'S')[2]) for s in evts.subject]
    evts.latency .+= size(data, 1) .* (subj_idx .- 1)


    m = fit(
        UnfoldModel,
        [
            Any => (
                @formula(0 ~ 1 + condition + zerocorr(1 + condition | subject)),
                firbasis([0.0, 1], 10; interpolate = false),
            ),
        ],
        evts,
        data,
    )
    @test_broken eff = effects(Dict(:condition => ["car", "face"]), m)
end
