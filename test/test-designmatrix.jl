@testset "combining MixedModel Designmatrices" begin

    basisfunction1 = firbasis(τ = (0, 1), sfreq = 10, name = "basis1")
    basisfunction2 = firbasis(τ = (0, 0.5), sfreq = 10, name = "basis2")

    tbl = DataFrame(
        [1 12 23 34 45 56 67 78; 1 1 1 2 2 2 3 3; 1 2 3 1 2 3 1 2]',
        [:latency, :subject, :item],
    )
    tbl2 = DataFrame(
        [2 9 16 23 30 37 44 51; 1 1 1 2 2 2 3 3; 1 2 3 1 2 3 1 2]',
        [:latency, :subject, :itemB],
    )
    y = Float64.([collect(range(1, stop = 100))...])'
    transform!(tbl, :subject => categorical => :subject)
    transform!(tbl2, :itemB => categorical => :itemB)
    transform!(tbl, :item => categorical => :item)
    #tbl.itemB = tbl.item
    f3 = @formula 0 ~ 1 + (1 | item) + (1 | subject)
    f4 = @formula 0 ~ 1 + (1 | itemB)
    f4_wrong = @formula 0 ~ 1 + (1 | item)

    Xdc3 = designmatrix(UnfoldLinearMixedModel, f3, tbl, basisfunction1)
    Xdc4 = designmatrix(UnfoldLinearMixedModel, f4, tbl2, basisfunction2)
    Xdc4_wrong = designmatrix(UnfoldLinearMixedModel, f4_wrong, tbl, basisfunction2)

    Xdc = Xdc3 + Xdc4
    @test typeof(modelmatrix(Xdc)[1]) <: SparseArrays.SparseMatrixCSC
    @test length(modelmatrix(Xdc)) == 4 # one FeMat  + 3 ReMat
    @test_throws String modelmatrix(Xdc3 + Xdc4_wrong)



end

@testset "equalizeReMatLengths" begin
    bf1 = firbasis(τ = (0, 1), sfreq = 10, name = "basis1")
    bf2 = firbasis(τ = (0, 0.5), sfreq = 10, name = "basis2")

    tbl1 = DataFrame(
        [1 12 23 34 45 56 67 78; 1 1 1 2 2 2 3 3; 1 2 3 1 2 3 1 2]',
        [:latency, :subject, :item],
    )
    tbl2 = DataFrame(
        [2 9 16 23 30 37 44 51; 1 1 1 2 2 2 3 3; 1 2 3 1 2 3 1 2]',
        [:latency, :subject, :itemB],
    )

    transform!(tbl1, :subject => categorical => :subject)
    transform!(tbl1, :item => categorical => :item)
    transform!(tbl2, :itemB => categorical => :itemB)
    #tbl.itemB = tbl.item
    f1 = @formula 0 ~ 1 + (1 | item) + (1 | subject)
    f2 = @formula 0 ~ 1 + (1 | itemB)

    form = Unfold.apply_schema(f1, Unfold.schema(f1, tbl1), MixedModels.LinearMixedModel)
    form = Unfold.apply_basisfunction(form, bf1, nothing, Any)
    X1 = Unfold.modelcols.(form.rhs, Ref(tbl1))

    form = Unfold.apply_schema(f2, Unfold.schema(f2, tbl2), MixedModels.LinearMixedModel)
    form = Unfold.apply_basisfunction(form, bf2, nothing, Any)
    X2 = Unfold.modelcols.(form.rhs, Ref(tbl2))

    # no missmatch, shouldnt change anything then
    X = deepcopy(X1[2:end])

    UnfoldMixedModels.equalize_ReMat_lengths!(X)
    @test all([x[1] for x in size.(X)] .== 88)

    X = (deepcopy(X1[2:end])..., deepcopy(X2[2:end])...)
    @test !all([x[1] for x in size.(X)] .== 56) # not alllenghts the same
    UnfoldMixedModels.equalize_ReMat_lengths!(X)
    @test all([x[1] for x in size.(X)] .== 88) # now all lengths the same :-)


    X = deepcopy(X2[2])

    @test size(X)[1] == 56
    UnfoldMixedModels.change_ReMat_size!(X, 60)
    @test size(X)[1] == 60

    X = deepcopy(X2[2])
    @test size(X)[1] == 56
    UnfoldMixedModels.change_ReMat_size!(X, 40)
    @test size(X)[1] == 40


    X = (deepcopy(X1)..., deepcopy(X2[2:end])...)
    @test size(X[1])[1] == 88
    @test size(X[2])[1] == 88
    @test size(X[3])[1] == 88
    @test size(X[4])[1] == 56
    XA, XB = UnfoldMixedModels.change_modelmatrix_size!(95, X[1], X[2:end])
    @test size(XA)[1] == 95
    @test size(XB)[1] == 95

    XA, XB = UnfoldMixedModels.change_modelmatrix_size!(40, X[1], X[2:end])
    @test size(XA)[1] == 40
    @test size(XB)[1] == 40

    XA, XB = UnfoldMixedModels.change_modelmatrix_size!(30, Matrix(X[1]), X[2:end])
    @test size(XA)[1] == 30
    @test size(XB)[1] == 30
end

@testset "designmatrix zerocorr, non-sequenial" begin

    data, evts =
        UnfoldSim.predef_2x2(; return_epoched = true, n_subjects = 5, noiselevel = 1)
    evts.subject = categorical(evts.subject)

    # shift each subject's events into its own non-overlapping data segment: the
    # time-expanded random effects require each time sample to be covered by at most one subject
    subj_idx = [parse(Int, split(string(s), 'S')[2]) for s in evts.subject]
    evts.latency .+= maximum(evts.latency) .* (subj_idx .- 1)
    sort!(evts, :latency)


    f_zc = @formula 0 ~ 1 + A + B + zerocorr(1 + A + B | subject)
    basisfunction = firbasis(τ = (-0.1, 0.1), sfreq = 10, name = "ABC")
    Xdc_zc = designmatrix(UnfoldLinearMixedModel, f_zc, evts, basisfunction)

    @test length(Xdc_zc.modelmatrix[2].inds) == 9
    f = @formula 0 ~ 1 + A + B + (1 + A + B | subject)
    Xdc = designmatrix(UnfoldLinearMixedModel, f, evts, basisfunction)
    @test length(Xdc.modelmatrix[2].inds) == (9 * 9 + 9) / 2

    # test bug with not sequential subjects
    evts_nonseq = copy(evts)
    evts_nonseq = evts_nonseq[.!(evts_nonseq.subject .== 2), :]
    Xdc_nonseq = designmatrix(UnfoldLinearMixedModel, f_zc, evts_nonseq, basisfunction)


end

@testset "random effects: overlap & refs/levels correspondence" begin
    bf = firbasis(τ = (0, 0.5), sfreq = 10) # 6 taps, span [latency, latency + 5]
    f = @formula 0 ~ 1 + (1 | subject)

    function get_remat(f, bf, tbl)
        form = Unfold.apply_schema(f, Unfold.schema(f, tbl), MixedModels.LinearMixedModel)
        form = Unfold.apply_basisfunction(form, bf, nothing, Any)
        return Unfold.modelcols.(form.rhs, Ref(tbl))[2]
    end

    # events of different subjects may not overlap in time (would leak random effects)
    tbl = DataFrame(latency = [20, 24], subject = categorical(["S1", "S2"]))
    @test_throws ErrorException get_remat(f, bf, tbl)

    # non-overlapping subjects: each time sample is owned by the level covering it
    tbl = DataFrame(latency = [10, 30], subject = categorical(["S1", "S2"]))
    re = get_remat(f, bf, tbl)
    @test re.levels == ["S1", "S2"]
    @test re.refs[10:15] == fill(Int32(1), 6) # design length = 30 + 6 - 1 = 35
    @test re.refs[30:35] == fill(Int32(2), 6)

    # events of the same subject may overlap
    tbl = DataFrame(latency = [20, 24, 40], subject = categorical(["S1", "S1", "S2"]))
    re = get_remat(f, bf, tbl)
    @test re.refs[20:29] == fill(Int32(1), 10) # overlapping S1 events at 20, 24
    @test re.refs[40:45] == fill(Int32(2), 6)

    # interleaved subjects in the event table are assigned correctly
    tbl = DataFrame(
        latency = [10, 20, 12, 30],
        subject = categorical(["S2", "S1", "S2", "S1"]),
    )
    re = get_remat(f, bf, tbl)
    @test re.refs[10:17] == fill(Int32(2), 8) # S2 events at 10, 12
    @test re.refs[20:25] == fill(Int32(1), 6) # S1 events at 20, 30
    @test re.refs[30:35] == fill(Int32(1), 6)

    # refs correspond to levels, not to the order of first appearance in the table
    tbl = DataFrame(latency = [10, 30], subject = categorical(["S2", "S1"]))
    levels!(tbl.subject, ["S1", "S2"])
    re = get_remat(f, bf, tbl)
    @test re.levels == ["S1", "S2"]
    @test re.refs[10:15] == fill(Int32(2), 6) # S2 = level 2
    @test re.refs[30:35] == fill(Int32(1), 6) # S1 = level 1

    # non-integer latencies are floored
    tbl = DataFrame(latency = [10.5, 30.0], subject = categorical(["S1", "S2"]))
    re = get_remat(f, bf, tbl)
    @test re.refs[10] == Int32(1) # floor(10.5) = 10
end
