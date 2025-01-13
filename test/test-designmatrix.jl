@testset "combining MixedModel Designmatrices" begin

    basisfunction1 = firbasis(τ = (0, 1), sfreq = 10, name = "basis1")
    basisfunction2 = firbasis(τ = (0, 0.5), sfreq = 10, name = "basis2")

    tbl = DataFrame(
        [1 4 10 15 20 22 31 37; 1 1 1 2 2 2 3 3; 1 2 3 1 2 3 1 2]',
        [:latency, :subject, :item],
    )
    tbl2 = DataFrame(
        [2 3 12 18 19 25 40 43; 1 1 1 2 2 2 3 3; 1 2 3 1 2 3 1 2]',
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
        [1 4 10 15 20 22 31 37; 1 1 1 2 2 2 3 3; 1 2 3 1 2 3 1 2]',
        [:latency, :subject, :item],
    )
    tbl2 = DataFrame(
        [2 3 12 18 19 25 40 43; 1 1 1 2 2 2 3 3; 1 2 3 1 2 3 1 2]',
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
    @test all([x[1] for x in size.(X)] .== 47)

    X = (deepcopy(X1[2:end])..., deepcopy(X2[2:end])...)
    @test !all([x[1] for x in size.(X)] .== 48) # not alllenghts the same
    UnfoldMixedModels.equalize_ReMat_lengths!(X)
    @test all([x[1] for x in size.(X)] .== 48) # now all lengths the same :-)


    X = deepcopy(X2[2])

    @test size(X)[1] == 48
    UnfoldMixedModels.change_ReMat_size!(X, 52)
    @test size(X)[1] == 52

    X = deepcopy(X2[2])
    @test size(X)[1] == 48
    UnfoldMixedModels.change_ReMat_size!(X, 40)
    @test size(X)[1] == 40


    X = (deepcopy(X1)..., deepcopy(X2[2:end])...)
    @test size(X[1])[1] == 47
    @test size(X[2])[1] == 47
    @test size(X[3])[1] == 47
    @test size(X[4])[1] == 48
    XA, XB = UnfoldMixedModels.change_modelmatrix_size!(52, X[1], X[2:end])
    @test size(XA)[1] == 52
    @test size(XB)[1] == 52

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
