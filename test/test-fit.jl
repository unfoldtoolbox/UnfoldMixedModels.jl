
@testset "autodetect type" begin
    @test Unfold.design_to_modeltype([Any => (@formula(0 ~ 1 + (1 | test)), 0:10)]) ==
          UnfoldLinearMixedModel
    @test Unfold.design_to_modeltype([
        Any => (@formula(0 ~ 1 + (1 | test)), firbasis(τ = (-1, 1), sfreq = 20)),
    ],) == UnfoldLinearMixedModelContinuousTime
end


@testset "lmm tests" begin
    ###############################
    ##  Mixed Model tests
    ###############################

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

    evts.C = rand(StableRNG(1), ["a", "b", "c"], size(evts, 1))
    data = reshape(data, 1, :)
    #append!(data, zeros(1000))
    data = vcat(data, data)
    #data = data .+ 1 * randn(size(data)) # we have to add minimal noise, else mixed models crashes.
    data_missing = Array{Union{Missing,Number}}(undef, size(data))
    data_missing .= deepcopy(data)

    data_missing[450:460] .= missing

    transform!(evts, :subject => categorical => :subject)

    f = @formula 0 ~ 1 + A + B + (1 + A + B | subject)
    #f  = @formula 0~1 + (1|subject)



    # cut the data into epochs
    # TODO This ignores subject bounds
    data_e, times = Unfold.epoch(data = data, tbl = evts, τ = (-1.0, 1.9), sfreq = 10)
    data_missing_e, times =
        Unfold.epoch(data = data_missing, tbl = evts, τ = (-1.0, 1.9), sfreq = 10)
    evts_e, data_e = Unfold.drop_missing_epochs(copy(evts), data_e)
    evts_missing_e, data_missing_e = Unfold.drop_missing_epochs(copy(evts), data_missing_e)

    ######################
    ##  Mass Univariate Mixed
    @time m_mum = fit(
        UnfoldModel,
        f,
        evts_e,
        data_e,
        times,
        contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding()),
        show_progress = false,
    )
    df = Unfold.coeftable(m_mum)
    @test isapprox(
        df[
            (df.channel .== 1) .& (df.coefname .== "A: a_small") .& (df.time .== 0.0),
            :estimate,
        ],
        [-0.02, 0.054],
        atol = 0.01,
    )



    # with missing
    @time m_mum = fit(
        UnfoldModel,
        f,
        evts_missing_e,
        data_missing_e,
        times,
        contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding()),
        show_progress = false,
    )
    df = coeftable(m_mum)
    @test isapprox(
        df[
            (df.channel .== 1) .& (df.coefname .== "A: a_small") .& (df.time .== 0.0),
            :estimate,
        ],
        [0.031, 0.05],
        atol = 0.1,
    )


    # Timexpanded Univariate Mixed
    f = @formula 0 ~ 1 + A + B + (1 + A | subject)
    basisfunction = firbasis(τ = (-0.2, 0.3), sfreq = 10)
    @time m_tum = fit(
        UnfoldModel,
        f,
        evts,
        data,
        basisfunction,
        contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding()),
        show_progress = false,
    )
    df = coeftable(m_tum)
    @test isapprox(
        df[
            (df.channel .== 1) .& (df.coefname .== "A: a_small") .& (df.time .== 0.0),
            :estimate,
        ],
        [-0.03, 0.064],
        atol = 0.1,
    )


    # missing data in LMMs
    # not yet implemented
    Test.@test_broken m_tum = fit(
        UnfoldModel,
        f,
        evts,
        data_missing,
        basisfunction,
        contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding()),
    )


    evts.subjectB = evts.subject
    evts1 = evts[evts.A .== "a_small", :]
    evts2 = evts[evts.A .== "a_big", :]

    f0_lmm = @formula 0 ~ 1 + B + (1 | subject) + (1 + C | subjectB)
    @time m = fit(UnfoldModel, f0_lmm, evts, data, basisfunction; show_progress = false)
    @time m_tum = coeftable(m)


    f1_lmm = @formula 0 ~ 1 + B + (1 | subject)
    f2_lmm = @formula 0 ~ 1 + B + (1 | subjectB)

    b1 = firbasis(τ = (-0.2, 0.3), sfreq = 10, name = "a_small")
    b2 = firbasis(τ = (-0.1, 0.3), sfreq = 10, name = "a_big")

    X1_lmm = designmatrix(UnfoldLinearMixedModelContinuousTime, f1_lmm, evts1, b1)
    X2_lmm = designmatrix(UnfoldLinearMixedModelContinuousTime, f2_lmm, evts2, b2)

    r = fit(
        UnfoldLinearMixedModelContinuousTime,
        X1_lmm + X2_lmm,
        data;
        show_progress = false,
    )
    df = coeftable(r)

    @test isapprox(
        df[
            (df.channel .== 1) .& (df.coefname .== "B: b_tiny") .& (df.time .== 0.0),
            :estimate,
        ],
        [0.65, 0.69],
        rtol = 0.1,
    )

    # Fast-lane new implementation
    m = coeftable(
        fit(
            UnfoldModel,
            ["a_small" => (f1_lmm, b1), "a_big" => (f2_lmm, b2)],
            evts,
            data,
            eventcolumn = "A",
        ),
    )


    #----
    # # test #13, 2x3 design

    f = @formula 0 ~ 1 + A + C + (1 + A + C | subject)
    #f  = @formula 0~1 + (1|subject)



    # cut the data into epochs
    # TODO This ignores subject bounds
    data_e, times = Unfold.epoch(data = data, tbl = evts, τ = (-1.0, 1.9), sfreq = 10)
    data_missing_e, times =
        Unfold.epoch(data = data_missing, tbl = evts, τ = (-1.0, 1.9), sfreq = 10)
    evts_e, data_e = Unfold.drop_missing_epochs(copy(evts), data_e)
    evts_missing_e, data_missing_e = Unfold.drop_missing_epochs(copy(evts), data_missing_e)

    ######################
    ##  Mass Univariate Mixed
    @time m_mum = fit(
        UnfoldModel,
        f,
        evts_e,
        data_e,
        times,
        contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding()),
        show_progress = false,
    )
    df = Unfold.coeftable(m_mum)
end
## Condense check for multi channel, multi
@testset "LMM multi channel, multi basisfunction" begin
    data, evts = UnfoldSim.predef_2x2(;
        return_epoched = true,
        n_subjects = 5,
        noiselevel = 1,
        signalsize = 10,
        n_items = 16,
    )
    subj_idx = [parse(Int, split(string(s), 'S')[2]) for s in evts.subject]
    evts.latency .+= size(data, 1) .* (subj_idx .- 1)

    data = reshape(data, 1, :)
    transform!(evts, :subject => categorical => :subject)
    data = vcat(data, data)

    bA0 = firbasis(τ = (-0.0, 0.1), sfreq = 10)
    bA1 = firbasis(τ = (0.1, 0.2), sfreq = 10)
    evts.subject2 = evts.subject
    fA0 = @formula 0 ~ 1 + B + zerocorr(1 | subject)
    fA1 = @formula 0 ~ 1 + B + zerocorr(1 | subject2)
    m = fit(
        UnfoldModel,
        ["a_small" => (fA0, bA0), "a_big" => (fA1, bA1)],
        evts,
        data;
        eventcolumn = "A",
        show_progress = false,
    )

    res = coeftable(m)

    @test all(last(.!isnothing.(res.group), 8))
    @test all(last(res.coefname, 8) .== "(Intercept)")

    # test more complex formulas
    fA0 = @formula 0 ~ 1 + zerocorr(1 + C | subject)
    fA1 = @formula 0 ~ 1 + B +C + zerocorr(1 + C | subject2)
    evts.C = rand(StableRNG(1), ["a", "b", "c"], size(evts, 1))
    m = fit(
        UnfoldModel,
        ["a_small" => (fA0, bA0), "a_big" => (fA1, bA1)],
        evts,
        data;
        eventcolumn = "A",
        show_progress = false,
    )

    res = coeftable(m)

end


@testset "LMM bug reorder #115" begin

    data, evts = UnfoldSim.predef_2x2(;
        return_epoched = true,
        n_subjects = 10,
        noiselevel = 1,
        onset = NoOnset(),
    )

    data = reshape(data, size(data, 1), :)

    designList = [
        [
            Any => (
                @formula(
                    0 ~ 1 + A + B + zerocorr(1 + B + A | subject) + zerocorr(1 + B | item)
                ),
                range(0, 1, length = size(data, 1)),
            ),
        ],
        [
            Any => (
                @formula(
                    0 ~ 1 + A + B + zerocorr(1 + A + B | subject) + zerocorr(1 + B | item)
                ),
                range(0, 1, length = size(data, 1)),
            ),
        ],
        [
            Any => (
                @formula(0 ~ 1 + zerocorr(1 + A + B | subject) + zerocorr(1 | item)),
                range(0, 1, length = size(data, 1)),
            ),
        ],
    ]
    #des = designList[1]
    #des = designList[2]
    for des in designList
        @test_throws AssertionError fit(UnfoldModel, des, evts, data)
        #
    end

    #counter check

    des = [
        Any => (
            @formula(0 ~ 1 + zerocorr(1 | item) + zerocorr(1 + A + B | subject)),
            range(0, 1, length = size(data, 1)),
        ),
    ]

    #= fails but not in the repl...?
    uf = fit(UnfoldModel, des, evts, data; show_progress = false)
    @test 3 ==
          unique(
        @subset(
            coeftable(uf),
            @byrow(:group == Symbol("subject")),
            @byrow :time == 0.0
        ).coefname,
    ) |> length
    =#
end


@testset "LMM bug reshape #110" begin
    data, evts =
        UnfoldSim.predef_2x2(; return_epoched = true, n_subjects = 10, noiselevel = 1)
    data = reshape(data, size(data, 1), :)
    subj_idx = [parse(Int, split(string(s), 'S')[2]) for s in evts.subject]
    evts.latency .+= size(data, 1) .* (subj_idx .- 1)

    des = [
        Any => (
            @formula(
                0 ~ 1 + A + B + zerocorr(1 + B + A | item) + zerocorr(1 + B | subject)
            ),
            range(0, 1, length = size(data, 1)),
        ),
    ]
    uf = fit(UnfoldModel, des, evts, data; show_progress = false)
    @test size(coef(uf)) == (1, 100, 3)
end
