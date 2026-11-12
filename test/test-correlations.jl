
@testset "LMM correlations in coeftable" begin
    data, evts = UnfoldSim.predef_2x2(
        StableRNG(1);
        return_epoched = false,
        n_subjects = 5,
        noiselevel = 1,
        signalsize = 10,
        n_items = 16,
    )
    subj_idx = [parse(Int, split(string(s), 'S')[2]) for s in evts.subject]
    evts.latency .+= size(data, 1) .* (subj_idx .- 1)
    data = reshape(data, 1, :)
    data = vcat(data, data)
    transform!(evts, :subject => categorical => :subject)

    data_e, times = Unfold.epoch(data = data, tbl = evts, τ = (-1.0, 1.9), sfreq = 10)
    evts_e, data_e = Unfold.drop_missing_epochs(copy(evts), data_e)

    f = @formula 0 ~ 1 + A + B + (1 + A + B | subject)
    m = fit(
        UnfoldModel,
        f,
        evts_e,
        data_e,
        times,
        contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding()),
        show_progress = false,
    )

    # tidyρs: one entry per fit & correlated pair
    corrs = tidyρs(m)
    @test length(corrs) == 3 * 2 * length(times) # 3 pairs x 2 channels x ntime
    @test all(
        haskey(c, :iter) &&
            haskey(c, :group) &&
            haskey(c, :column1) &&
            haskey(c, :column2) &&
            haskey(c, :ρ) for c in corrs
    )
    @test all(c.group == :subject for c in corrs)
    @test all(-1 .<= c.ρ .<= 1 for c in corrs)

    # correlations are part of the coeftable
    df = coeftable(m)
    isρ(name) = occursin("ρ", string(name))
    corr_rows = df[isρ.(string.(df.coefname)), :]
    @test nrow(corr_rows) == 3 * 2 * length(times)
    @test all(corr_rows.group .== :subject)
    @test all(isnothing.(corr_rows.stderror))

    # values match MixedModels' own extraction for the first fit (channel 1, first timepoint)
    evts_y = copy(evts_e)
    evts_y.y = data_e[1, 1, :]
    m_mm = fit(
        MixedModel,
        @formula(y ~ 1 + A + B + (1 + A + B | subject)),
        evts_y;
        contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding()),
    )
    ρ_mm = collect(MixedModels.σρs(m_mm)[1].ρ)
    our = corr_rows[(corr_rows.channel .== 1) .& (corr_rows.time .== -1.0), :estimate]
    @test isapprox(ρ_mm, our, atol = 1e-6)

    # zerocorr: no correlations, table unchanged otherwise
    f0 = @formula 0 ~ 1 + A + B + zerocorr(1 + A + B | subject)
    m0 = fit(
        UnfoldModel,
        f0,
        evts_e,
        data_e,
        times,
        contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding()),
        show_progress = false,
    )
    df0 = coeftable(m0)
    @test nrow(df0[isρ.(string.(df0.coefname)), :]) == 0
    @test nrow(df0) == nrow(df) - nrow(corr_rows)

    # single random effect column: no correlations
    f1 = @formula 0 ~ 1 + A + B + (1 | subject)
    m1 = fit(
        UnfoldModel,
        f1,
        evts_e,
        data_e,
        times,
        contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding()),
        show_progress = false,
    )
    @test nrow(coeftable(m1)[isρ.(string.(coeftable(m1).coefname)), :]) == 0

    # continuous-time model: correlations across basisfunction columns, one row per pair & channel
    f2 = @formula 0 ~ 1 + A + B + (1 + A | subject)
    m2 = fit(
        UnfoldModel,
        f2,
        evts,
        data,
        firbasis(τ = (-0.2, 0.3), sfreq = 10),
        contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding()),
        show_progress = false,
    )
    df2 = coeftable(m2)
    corr_rows2 = df2[isρ.(string.(df2.coefname)), :]
    @test nrow(corr_rows2) == 66 * 2 # 12 timeexpanded re-columns: 66 pairs x 2 channels
    @test all(corr_rows2.group .== :subject)
    @test all(isnan.(corr_rows2.time))
end

@testset "LMM correlations, two random groupings" begin
    data, evts = UnfoldSim.predef_2x2(
        StableRNG(1);
        return_epoched = false,
        n_subjects = 5,
        noiselevel = 1,
        signalsize = 10,
        n_items = 16,
    )
    subj_idx = [parse(Int, split(string(s), 'S')[2]) for s in evts.subject]
    evts.latency .+= size(data, 1) .* (subj_idx .- 1)
    data = reshape(data, 1, :)
    data = vcat(data, data)
    transform!(evts, :subject => categorical => :subject)

    data_e, times = Unfold.epoch(data = data, tbl = evts, τ = (-1.0, 1.9), sfreq = 10)
    evts_e, data_e = Unfold.drop_missing_epochs(copy(evts), data_e)

    # absurd is a purely random grouping, unrelated to the signal
    evts2 = copy(evts_e)
    evts2.absurd = string.(rand(StableRNG(1), 1:10, size(evts2, 1)))
    transform!(evts2, :absurd => categorical => :absurd)

    f = @formula 0 ~ 1 + A + B + (1 + A | absurd) + (1 + B | item)
    m = fit(
        UnfoldModel,
        f,
        evts2,
        data_e,
        times,
        contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding()),
        show_progress = false,
    )

    isρ(name) = occursin("ρ", string(name))
    corr_rows = coeftable(m)[isρ.(string.(coeftable(m).coefname)), :]

    # one pair per group: 2 groups x 2 channels x ntime
    @test nrow(corr_rows) == 2 * 2 * length(times)
    @test sort(unique(corr_rows.group)) == [:absurd, :item]
    @test all(isnothing.(corr_rows.stderror))

    # values match MixedModels' own extraction for the first fit (channel 1, first timepoint);
    # this catches θ-offset mix-ups between the two random-effect blocks
    evts_y = copy(evts2)
    evts_y.y = data_e[1, 1, :]
    m_mm = fit(
        MixedModel,
        @formula(y ~ 1 + A + B + (1 + A | absurd) + (1 + B | item)),
        evts_y;
        contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding()),
    )
    σρ = MixedModels.σρs(m_mm)
    for (g, ρ_mm) in zip(keys(σρ), σρ)
        our = corr_rows[
            (corr_rows.group .== g) .& (corr_rows.channel .== 1) .& (corr_rows.time .== -1.0),
            :estimate,
        ]
        @test length(our) == 1 # one correlation per group
        @test all(
            isapprox.(collect(ρ_mm.ρ), our, atol = 1e-6) .|
            (isnan.(collect(ρ_mm.ρ)) .& isnan.(our)),
        )
    end
end
