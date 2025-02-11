var documenterSearchIndex = {"docs":
[{"location":"howto/lmm_pvalues/#lmm_pvalues","page":"P-values for mixedModels","title":"How To get P-Values for Mass-Univariate LMM","text":"","category":"section"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"There are currently two ways to obtain p-values for LMMs: Wald's t-test and likelihood ratio tests (mass univariate only).","category":"page"},{"location":"howto/lmm_pvalues/#Setup","page":"P-values for mixedModels","title":"Setup","text":"","category":"section"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"using UnfoldMixedModels # we require to load MixedModels to load the PackageExtension\nusing DataFrames\nusing UnfoldSim\nusing CairoMakie\nusing DisplayAs # hide\ndata_epoch, evts =\n    UnfoldSim.predef_2x2(; n_items = 52, n_subjects = 40, return_epoched = true)\ndata_epoch = reshape(data_epoch, size(data_epoch, 1), :) #\ntimes = range(0, 1, length = size(data_epoch, 1))","category":"page"},{"location":"howto/lmm_pvalues/#Define-f0-and-f1-and-fit","page":"P-values for mixedModels","title":"Define f0 & f1 and fit","text":"","category":"section"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"\nf0 = @formula 0 ~ 1 + A + (1 + A | subject);\nf1 = @formula 0 ~ 1 + A + B + (1 + A | subject); # could also differ in random effects\n\nm0 = fit(UnfoldModel,[Any=>(f0,times)],evts,data_epoch);\nm1 = fit(UnfoldModel,[Any=>(f1,times)],evts,data_epoch);\n\nm1|> DisplayAs.withcontext(:is_pluto=>true) # hide","category":"page"},{"location":"howto/lmm_pvalues/#Likelihood-ratio","page":"P-values for mixedModels","title":"Likelihood ratio","text":"","category":"section"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"uf_lrt = likelihoodratiotest(data_epoch, m0, m1)\nuf_lrt[1]","category":"page"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"As you can see, we have some likelihood ratio outcomes, exciting!","category":"page"},{"location":"howto/lmm_pvalues/#Extract-p-values","page":"P-values for mixedModels","title":"Extract p-values","text":"","category":"section"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"pvalue(uf_lrt)","category":"page"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"We have extracted the p-values and now need to make them usable.     The solution can be found in the documentation under ?pvalue.","category":"page"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"pvals_lrt = vcat(pvalue(uf_lrt)...)\nnchan = 1\nntime = length(times)\nreshape(pvals_lrt, ntime, nchan)' # note the last transpose via ' !","category":"page"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"Perfecto, these are the LRT p-values of a model condA vs. condA+condB with same random effect structure.","category":"page"},{"location":"howto/lmm_pvalues/#Walds-T-Test","page":"P-values for mixedModels","title":"Walds T-Test","text":"","category":"section"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"This method is easier to calculate but has limitations in accuracy and scope. It may also be less accurate due to the liberal estimation of degrees of freedom. Testing is limited in this case, as random effects cannot be tested and only single predictors can be used, which may not be appropriate for spline effects. It is important to note that this discussion is beyond the scope of this LMM package.","category":"page"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"res = coeftable(m1)\n# only fixed effects: what is not in a ranef group is a fixef.\nres = res[isnothing.(res.group), :]\n# calculate z-value\nres[:, :zvalue] = res.estimate ./ res.stderror","category":"page"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"We obtained Walds z, but how to translate them to a p-value?","category":"page"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"Determining the necessary degrees of freedom for the z/t-distribution is a complex issue with much debate surrounding it. One approach is to use the number of subjects as an upper bound for the p-value (your df will be between n_subject and sumn_trials).","category":"page"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"df = length(unique(evts.subject))","category":"page"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"Plug it into the t-distribution.","category":"page"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"using Distributions\nres.pvalue = pdf.(TDist(df),res.zvalue)","category":"page"},{"location":"howto/lmm_pvalues/#Comparison-of-methods","page":"P-values for mixedModels","title":"Comparison of methods","text":"","category":"section"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"Cool! Let's compare both methods of p-value calculation!","category":"page"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"df = DataFrame(:walds => res[res.coefname.==\"B: b_tiny\", :pvalue], :lrt => pvals_lrt)\nf = Figure()\n\nscatter(f[1,1],times,res[res.coefname .== \"B: b_tiny\",:estimate],axis=(;xlabel=\"time\",title=\"coef: B:b_tiny\"))\nscatter(f[1,2],df.walds,df.lrt,axis=(;xlabel=\"walds-z pvalue\",ylabel=\"LRT pvalue\"))\nscatter(f[2,1],times,df.walds,axis=(;title=\"walds-z pvalue\",xlabel=\"time\"))\nscatter(f[2,2],times,df.lrt,axis=(;title=\"lrt pvalue\",xlabel=\"time\"))\n\nf","category":"page"},{"location":"howto/lmm_pvalues/","page":"P-values for mixedModels","title":"P-values for mixedModels","text":"Note that the Walds-z is typically too liberal (LRT also, but to a lesser exted). Best is to use the forthcoming MixedModelsPermutations.jl or go the route via R and use KenwardRogers (data not yet published)","category":"page"},{"location":"91-developer/#dev_docs","page":"Developer documentation","title":"Developer documentation","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"note: Contributing guidelines\nIf you haven't, please read the Contributing guidelines first.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If you want to make contributions to this package that involves code, then this guide is for you.","category":"page"},{"location":"91-developer/#First-time-clone","page":"Developer documentation","title":"First time clone","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"tip: If you have writing rights\nIf you have writing rights, you don't have to fork. Instead, simply clone and skip ahead. Whenever upstream is mentioned, use origin instead.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If this is the first time you work with this repository, follow the instructions below to clone the repository.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Fork this repo\nClone your repo (this will create a git remote called origin)\nAdd this repo as a remote:\ngit remote add upstream https://github.com/unfoldtoolbox/UnfoldMixedModels.jl","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"This will ensure that you have two remotes in your git: origin and upstream. You will create branches and push to origin, and you will fetch and update your local main branch from upstream.","category":"page"},{"location":"91-developer/#Linting-and-formatting","page":"Developer documentation","title":"Linting and formatting","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Install a plugin on your editor to use EditorConfig. This will ensure that your editor is configured with important formatting settings.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"We use https://pre-commit.com to run the linters and formatters. In particular, the Julia code is formatted using JuliaFormatter.jl, so please install it globally first:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"julia> # Press ]\npkg> activate\npkg> add JuliaFormatter","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To install pre-commit, we recommend using pipx as follows:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"# Install pipx following the link\npipx install pre-commit","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"With pre-commit installed, activate it as a pre-commit hook:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"pre-commit install","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To run the linting and formatting manually, enter the command below:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"pre-commit run -a","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Now, you can only commit if all the pre-commit tests pass.","category":"page"},{"location":"91-developer/#Testing","page":"Developer documentation","title":"Testing","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"As with most Julia packages, you can just open Julia in the repository folder, activate the environment, and run test:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"julia> # press ]\npkg> activate .\npkg> test","category":"page"},{"location":"91-developer/#Working-on-a-new-issue","page":"Developer documentation","title":"Working on a new issue","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"We try to keep a linear history in this repo, so it is important to keep your branches up-to-date.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Fetch from the remote and fast-forward your local main\ngit fetch upstream\ngit switch main\ngit merge --ff-only upstream/main\nBranch from main to address the issue (see below for naming)\ngit switch -c 42-add-answer-universe\nPush the new local branch to your personal remote repository\ngit push -u origin 42-add-answer-universe\nCreate a pull request to merge your remote branch into the org main.","category":"page"},{"location":"91-developer/#Branch-naming","page":"Developer documentation","title":"Branch naming","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If there is an associated issue, add the issue number.\nIf there is no associated issue, and the changes are small, add a prefix such as \"typo\", \"hotfix\", \"small-refactor\", according to the type of update.\nIf the changes are not small and there is no associated issue, then create the issue first, so we can properly discuss the changes.\nUse dash separated imperative wording related to the issue (e.g., 14-add-tests, 15-fix-model, 16-remove-obsolete-files).","category":"page"},{"location":"91-developer/#Commit-message","page":"Developer documentation","title":"Commit message","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Use imperative or present tense, for instance: Add feature or Fix bug.\nHave informative titles.\nWhen necessary, add a body with details.\nIf there are breaking changes, add the information to the commit message.","category":"page"},{"location":"91-developer/#Before-creating-a-pull-request","page":"Developer documentation","title":"Before creating a pull request","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"tip: Atomic git commits\nTry to create \"atomic git commits\" (recommended reading: The Utopic Git History).","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Make sure the tests pass.\nMake sure the pre-commit tests pass.\nFetch any main updates from upstream and rebase your branch, if necessary:\ngit fetch upstream\ngit rebase upstream/main BRANCH_NAME\nThen you can open a pull request and work with the reviewer to address any issues.","category":"page"},{"location":"91-developer/#Building-and-viewing-the-documentation-locally","page":"Developer documentation","title":"Building and viewing the documentation locally","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Following the latest suggestions, we recommend using LiveServer to build the documentation. Here is how you do it:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Run julia --project=docs to open Julia in the environment of the docs.\nIf this is the first time building the docs\nPress ] to enter pkg mode\nRun pkg> dev . to use the development version of your package\nPress backspace to leave pkg mode\nRun julia> using LiveServer\nRun julia> servedocs()","category":"page"},{"location":"references/functions/","page":"API: Functions","title":"API: Functions","text":"Modules = [UnfoldMixedModels]\nOrder   = [:function]","category":"page"},{"location":"references/functions/#MixedModels.likelihoodratiotest-Tuple{AbstractArray, Vararg{UnfoldLinearMixedModel}}","page":"API: Functions","title":"MixedModels.likelihoodratiotest","text":"likelihoodratiotest(data, m)\n\n\nCalculate likelihoodratiotest\n\n\n\n\n\n","category":"method"},{"location":"references/functions/#StatsAPI.fit!-Union{Tuple{T}, Tuple{Union{UnfoldLinearMixedModel, UnfoldLinearMixedModelContinuousTime}, AbstractArray{T}}} where T","page":"API: Functions","title":"StatsAPI.fit!","text":"fit!(uf::UnfoldModel,data::Union{<:AbstractArray{T,2},<:AbstractArray{T,3}}) where {T<:Union{Missing, <:Number}}\n\nFit a DesignMatrix against a 2D/3D Array data along its last dimension Data is typically interpreted as channel x time (with basisfunctions) or channel x time x epoch (for mass univariate)\n\nshow_progress (default:true), deactivate the progressmeter\n\nReturns an UnfoldModel object\n\nExamples\n\n\n\n\n\n\n\n","category":"method"},{"location":"references/functions/#StatsAPI.pvalue-Tuple{Vector{MixedModels.LikelihoodRatioTest}}","page":"API: Functions","title":"StatsAPI.pvalue","text":"pvalue(lrtvec)\n\n\nUnfold-Method: return pvalues of likelihoodratiotests, typically calculated:\n\nExamples\n\njulia> pvalue(likelihoodratiotest(m1,m2))\n\nwhere m1/m2 are UnfoldLinearMixedModel's\n\nTipp: if you only compare two models you can easily get a vector of p-values:\n\njulia> vcat(pvalues(likelihoodratiotest(m1,m2))...)\n\nMultiple channels are returned linearized at the moment, as we do not have access to the amount of channels after the LRT, you can do:\n\njulia> reshape(vcat(pvalues(likelihoodratiotest(m1,m2))...),ntimes,nchan)'\n\n\n\n\n\n","category":"method"},{"location":"references/functions/#StatsModels.modelcols-Tuple{Unfold.TimeExpandedTerm{<:Union{var\"#s1\", var\"#s118\"} where {var\"#s1\"<:RandomEffectsTerm, var\"#s118\"<:MixedModels.AbstractReTerm}}, Any}","page":"API: Functions","title":"StatsModels.modelcols","text":"modelcols(term, tbl)\n\n\nThis function timeexpands the random effects and generates a ReMat object\n\n\n\n\n\n","category":"method"},{"location":"references/functions/#Unfold.make_estimate-Tuple{Union{UnfoldLinearMixedModel, UnfoldLinearMixedModelContinuousTime}}","page":"API: Functions","title":"Unfold.make_estimate","text":"Unfold.make_estimate(m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},\n\n) extracts betas (and sigma's for mixed models) with string grouping indicator\n\nreturns as a ch x beta, or ch x time x beta (for mass univariate)\n\n\n\n\n\n","category":"method"},{"location":"references/functions/#Unfold.modelmatrices-Tuple{Tuple}","page":"API: Functions","title":"Unfold.modelmatrices","text":"modelmatrices(modelmatrix::Tuple)\n\nin the case of a Tuple (MixedModels - FeMat/ReMat Tuple), returns only the FeMat part\n\n\n\n\n\n","category":"method"},{"location":"references/functions/#UnfoldMixedModels.LinearMixedModel_wrapper-Union{Tuple{TData}, Tuple{Any, AbstractVector{<:TData}, Any}} where TData<:Number","page":"API: Functions","title":"UnfoldMixedModels.LinearMixedModel_wrapper","text":"LinearMixedModel_wrapper(form, data, Xs; wts)\n\n\nWrapper to generate a LinearMixedModel. Code taken from MixedModels.jl and slightly adapted.\n\n\n\n\n\n","category":"method"},{"location":"references/functions/#UnfoldMixedModels.fake_lmm-Union{Tuple{N}, Tuple{AbstractArray{<:Number, N}, UnfoldLinearMixedModel, Int64}} where N","page":"API: Functions","title":"UnfoldMixedModels.fake_lmm","text":"fake_lmm(data, m, k)\n\n\nReturns a partial LMM model (non-functional due to lacking data) to be used in likelihoodratiotests. k to selcet which of the modelfit's to fake\n\n\n\n\n\n","category":"method"},{"location":"references/functions/#UnfoldMixedModels.get_timeexpanded_random_grouping-Tuple{Any, Any, Any}","page":"API: Functions","title":"UnfoldMixedModels.get_timeexpanded_random_grouping","text":"get_timeexpanded_random_grouping(\n    tbl_group,\n    tbl_latencies,\n    basisfunction\n)\n\n\nGet the timeranges where the random grouping variable was applied\n\n\n\n\n\n","category":"method"},{"location":"references/functions/#UnfoldMixedModels.isa_lmm_formula-Tuple{typeof(zerocorr)}","page":"API: Functions","title":"UnfoldMixedModels.isa_lmm_formula","text":"isa_lmm_formula\n\niterates over all parts of a formula until either a MixedModels.zerocorr, or a | was found. Then returns true, else returns false.\n\n\n\n\n\n","category":"method"},{"location":"references/functions/#UnfoldMixedModels.random_effect_groupings-Tuple{StatsModels.AbstractTerm}","page":"API: Functions","title":"UnfoldMixedModels.random_effect_groupings","text":"random_effect_groupings(t::MixedModels.AbstractReTerm)\n\nReturns the random effect grouping term (rhs), similar to coefnames, which returns the left hand sides\n\n\n\n\n\n","category":"method"},{"location":"references/functions/#UnfoldMixedModels.reorder_tidyσs-Tuple{Any, Any}","page":"API: Functions","title":"UnfoldMixedModels.reorder_tidyσs","text":"reorder_tidyσs(t, f)\n\nThis function reorders a MixedModels.tidyσs output, according to the formula and not according to the largest RandomGrouping.\n\n\n\n\n\n","category":"method"},{"location":"references/types/","page":"API: Types","title":"API: Types","text":"Modules = [UnfoldMixedModels]\nOrder   = [:type]","category":"page"},{"location":"references/types/#UnfoldMixedModels.UnfoldLinearMixedModel","page":"API: Types","title":"UnfoldMixedModels.UnfoldLinearMixedModel","text":"Concrete type to implement an Mass-Univariate LinearMixedModel. .design contains the formula + times dict .designmatrix contains a DesignMatrix modelfit is a Any container for the model results\n\n\n\n\n\n","category":"type"},{"location":"references/types/#UnfoldMixedModels.UnfoldLinearMixedModelContinuousTime","page":"API: Types","title":"UnfoldMixedModels.UnfoldLinearMixedModelContinuousTime","text":"Concrete type to implement an deconvolution LinearMixedModel.\n\nWarning This is to be treated with care, not much testing went into it.\n\n.design contains the formula + times dict .designmatrix contains a DesignMatrix .modelfit is a Any container for the model results\n\n\n\n\n\n","category":"type"},{"location":"90-contributing/#contributing","page":"Contributing guidelines","title":"Contributing guidelines","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"First of all, thanks for the interest!","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"We welcome all kinds of contribution, including, but not limited to code, documentation, examples, configuration, issue creating, etc.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"Be polite and respectful, and follow the code of conduct.","category":"page"},{"location":"90-contributing/#Bug-reports-and-discussions","page":"Contributing guidelines","title":"Bug reports and discussions","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If you think you found a bug, feel free to open an issue. Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.","category":"page"},{"location":"90-contributing/#Working-on-an-issue","page":"Contributing guidelines","title":"Working on an issue","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If you found an issue that interests you, comment on that issue what your plans are. If the solution to the issue is clear, you can immediately create a pull request (see below). Otherwise, say what your proposed solution is and wait for a discussion around it.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"tip: Tip\nFeel free to ping us after a few days if there are no responses.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If your solution involves code (or something that requires running the package locally), check the developer documentation. Otherwise, you can use the GitHub interface directly to create your pull request.","category":"page"},{"location":"tutorials/lmm_overlap/#lmm_overlap","page":"lmmERP (overlap correction)","title":"Overlap Correction with Linear Mixed Models","text":"","category":"section"},{"location":"tutorials/lmm_overlap/","page":"lmmERP (overlap correction)","title":"lmmERP (overlap correction)","text":"using UnfoldMixedModels\n\nusing UnfoldSim\n\nusing CategoricalArrays\nusing UnfoldMakie, CairoMakie\nusing DataFrames\n\nnothing;#hide","category":"page"},{"location":"tutorials/lmm_overlap/","page":"lmmERP (overlap correction)","title":"lmmERP (overlap correction)","text":"This notebook is similar to the Linear Model with Overlap Correction tutorial, but fits mixed models with overlap correction","category":"page"},{"location":"tutorials/lmm_overlap/","page":"lmmERP (overlap correction)","title":"lmmERP (overlap correction)","text":"warning: Warning\nLimitation: This functionality is not ready for general use. There are still a lot of things to find out and tinker with. Don't use this if you haven't looked under the hood of the toolbox! Be aware of crashes / timeouts for non-trivial problems","category":"page"},{"location":"tutorials/lmm_overlap/#Get-some-data","page":"lmmERP (overlap correction)","title":"Get some data","text":"","category":"section"},{"location":"tutorials/lmm_overlap/","page":"lmmERP (overlap correction)","title":"lmmERP (overlap correction)","text":"dat, evts = UnfoldSim.predef_2x2(; signalsize=20, n_items=16, n_subjects=16)\n\n# We also need to fix the latencies, they are now relative to 1:size(data, 1), but we want a continuous long EEG.\nsubj_idx = [parse(Int, split(string(s), 'S')[2]) for s in evts.subject]\nevts.latency .+= size(dat, 1) .* (subj_idx .- 1)\n\ndat = dat[:] # we need all data concatenated over subjects\nevts.subject  = categorical(Array(evts.subject))\nnothing #hide","category":"page"},{"location":"tutorials/lmm_overlap/#Linear-**Mixed**-Model-Continuous-Time","page":"lmmERP (overlap correction)","title":"Linear Mixed Model Continuous Time","text":"","category":"section"},{"location":"tutorials/lmm_overlap/","page":"lmmERP (overlap correction)","title":"lmmERP (overlap correction)","text":"Again we have 4 steps:","category":"page"},{"location":"tutorials/lmm_overlap/","page":"lmmERP (overlap correction)","title":"lmmERP (overlap correction)","text":"Specify a temporal basisfunction\nSpecify a formula\nFit a linear model for each channel (one model for all timepoints!)\nVisualize the results.","category":"page"},{"location":"tutorials/lmm_overlap/#1.-Specify-a-temporal-basisfunction","page":"lmmERP (overlap correction)","title":"1. Specify a temporal basisfunction","text":"","category":"section"},{"location":"tutorials/lmm_overlap/","page":"lmmERP (overlap correction)","title":"lmmERP (overlap correction)","text":"By default, we would want to use a FIR basis function.","category":"page"},{"location":"tutorials/lmm_overlap/","page":"lmmERP (overlap correction)","title":"lmmERP (overlap correction)","text":"basisfunction = firbasis(τ=(-0.4, .8), sfreq=20, name=\"stimulus\")\nnothing #hide","category":"page"},{"location":"tutorials/lmm_overlap/#2.-Specify-the-formula","page":"lmmERP (overlap correction)","title":"2. Specify the formula","text":"","category":"section"},{"location":"tutorials/lmm_overlap/","page":"lmmERP (overlap correction)","title":"lmmERP (overlap correction)","text":"Define the formula and specify a random effect.","category":"page"},{"location":"tutorials/lmm_overlap/","page":"lmmERP (overlap correction)","title":"lmmERP (overlap correction)","text":"note: Note\nWe use zerocorr to prevent the model from computing all correlations between all timepoints and factors.","category":"page"},{"location":"tutorials/lmm_overlap/","page":"lmmERP (overlap correction)","title":"lmmERP (overlap correction)","text":"f  = @formula 0 ~ 1 + A  *B + zerocorr(1 + A*B|subject);","category":"page"},{"location":"tutorials/lmm_overlap/#3.-Fit-the-model","page":"lmmERP (overlap correction)","title":"3. Fit the model","text":"","category":"section"},{"location":"tutorials/lmm_overlap/","page":"lmmERP (overlap correction)","title":"lmmERP (overlap correction)","text":"bfDict = [Any=>(f, basisfunction)]\n# Skipping this tutorial for now due to a significant error.\nm = fit(UnfoldModel, bfDict, evts, dat)\n\nresults = coeftable(m)\nfirst(results, 6)","category":"page"},{"location":"tutorials/lmm_overlap/#4.-Visualize-results","page":"lmmERP (overlap correction)","title":"4. Visualize results","text":"","category":"section"},{"location":"tutorials/lmm_overlap/","page":"lmmERP (overlap correction)","title":"lmmERP (overlap correction)","text":"plot_erp(results; mapping=(; col = :group))","category":"page"},{"location":"","page":"UnfoldMixedModels.jl Documentation","title":"UnfoldMixedModels.jl Documentation","text":"CurrentModule = UnfoldMixedModels","category":"page"},{"location":"#UnfoldMixedModels.jl-Documentation","page":"UnfoldMixedModels.jl Documentation","title":"UnfoldMixedModels.jl Documentation","text":"","category":"section"},{"location":"","page":"UnfoldMixedModels.jl Documentation","title":"UnfoldMixedModels.jl Documentation","text":"Welcome to UnfoldMixedModels.jl: a Julia package to analyse timeseries with Linear Mixed Models. This is an standalone-addon to Unfold.jl with similar syntax, but optimized for LMMs / Hierarchical Models.","category":"page"},{"location":"","page":"UnfoldMixedModels.jl Documentation","title":"UnfoldMixedModels.jl Documentation","text":"<div style=\"width:60%; margin: auto;\">\n\n<img src=\"assets/UnfoldSim_features_animation.gif\"/>\n</div>","category":"page"},{"location":"#Key-features","page":"UnfoldMixedModels.jl Documentation","title":"Key features","text":"","category":"section"},{"location":"","page":"UnfoldMixedModels.jl Documentation","title":"UnfoldMixedModels.jl Documentation","text":"MixedModels.jl: Powered by the up-to-100x faster implementation of MixedModels.jl\nSubject and Item effects: Fit the full spectrum of LMMs, random slopes and all!\nBeta: clusterpermutation: Combine with UnfoldStats.jl and fit LMM clusterpermutation tests\nAlpha: Overlap: Model overlap and LMMs (experimental!)","category":"page"},{"location":"#Installation","page":"UnfoldMixedModels.jl Documentation","title":"Installation","text":"","category":"section"},{"location":"","page":"UnfoldMixedModels.jl Documentation","title":"UnfoldMixedModels.jl Documentation","text":"julia> using Pkg; Pkg.add(\"UnfoldMixedModels\")","category":"page"},{"location":"#Usage-example","page":"UnfoldMixedModels.jl Documentation","title":"Usage example","text":"","category":"section"},{"location":"","page":"UnfoldMixedModels.jl Documentation","title":"UnfoldMixedModels.jl Documentation","text":"using UnfoldMixedModels\nusing UnfoldSim\ndata, evts = UnfoldSim.predef_eeg(10;return_epoched=true) # 10 subjects\ndata = reshape(data,size(data,1),:) # concatenate subjects\n\ntimes = range(-0.1,0.5,size(data,1)) # arbitrary time-vector\n\nfLMM = @formula 0 ~ 1 + condition + (1 + condition|subject) + (1|item)\nfit(UnfoldModel, [Any=>(f, times)], evts, data)\nnothing #hide","category":"page"},{"location":"#Where-to-start:-Learning-roadmap","page":"UnfoldMixedModels.jl Documentation","title":"Where to start: Learning roadmap","text":"","category":"section"},{"location":"#0.-First-first-steps","page":"UnfoldMixedModels.jl Documentation","title":"0. First first steps","text":"","category":"section"},{"location":"","page":"UnfoldMixedModels.jl Documentation","title":"UnfoldMixedModels.jl Documentation","text":"📌 Goal: Familiarize yourself with rERP Unfold.jl fitting & MixedModels.jl 🔗 Unfold.jl Quickstart | MixedModels.jl","category":"page"},{"location":"#1.-First-steps","page":"UnfoldMixedModels.jl Documentation","title":"1. First steps","text":"","category":"section"},{"location":"","page":"UnfoldMixedModels.jl Documentation","title":"UnfoldMixedModels.jl Documentation","text":"📌 Goal: Learn to fit a mass univariate Linear Mixed Model 🔗 Mass Univariate Linear Mixed Models","category":"page"},{"location":"#2.-Intermediate-topics","page":"UnfoldMixedModels.jl Documentation","title":"2. Intermediate topics","text":"","category":"section"},{"location":"","page":"UnfoldMixedModels.jl Documentation","title":"UnfoldMixedModels.jl Documentation","text":"📌 Goal: Extract statistics and cluster permutation tests 🔗 How To get P-Values for Mass-Univariate LMM | LMM Cluster Permutation tests","category":"page"},{"location":"","page":"UnfoldMixedModels.jl Documentation","title":"UnfoldMixedModels.jl Documentation","text":"We further recommend to skim this online book (from the authors of MixedModels.jl): embraceuncertaintybook.com/","category":"page"},{"location":"","page":"UnfoldMixedModels.jl Documentation","title":"UnfoldMixedModels.jl Documentation","text":"You should also learn about model simplification (keep it maximal ..?) and contrast codings. Enjoy!","category":"page"},{"location":"#Statement-of-need","page":"UnfoldMixedModels.jl Documentation","title":"Statement of need","text":"","category":"section"},{"location":"","page":"UnfoldMixedModels.jl Documentation","title":"UnfoldMixedModels.jl Documentation","text":"EEG researchers often analyse complex experimental procedures and want to generalize them to populations of subjects, items, schools etc. (Generalization Crisis - Yarkoni 2020). In case of hierarchical structures (e.g. repeated trials in subjects, different stimuli used) the Linear Mixed Model has become very popular. Unfortunately, fitting such models can be quite involved, especially for EEG data which require massive-modelfitting for each sensor and channel. MixedModels.jl provides a fast way for fitting, and UnfoldMixedModels.jl provides the bookkeeping to do so in a massive way. Users can easily extract fixed and random effects over time and sensors, do statistical testing and even correction for multiple comparisons (via UnfoldStats.jl / MixedModelsPermutations.jl).","category":"page"},{"location":"","page":"UnfoldMixedModels.jl Documentation","title":"UnfoldMixedModels.jl Documentation","text":"<!---\nNote: The statement of need is also used in the `README.md`. Make sure that they are synchronized.\n-->","category":"page"},{"location":"tutorials/lmm_mu/#lmm_massunivariate","page":"lmmERP (mass univariate)","title":"Mass Univariate Linear Mixed Models","text":"","category":"section"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"\nusing UnfoldMixedModels\nusing UnfoldSim\n\nusing UnfoldMakie, CairoMakie # plotting\nusing DataFrames\nusing CategoricalArrays\nnothing;#hide","category":"page"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"This notebook is similar to the Unfold.jl Mass Univariate Linear Models (no overlap correction) tutorial, but fits mass-univariate mixed models - that is, one model over all subjects, instead of one model per subject. This allows to include item effects, for example.","category":"page"},{"location":"tutorials/lmm_mu/#Mass-Univariate-**Mixed**-Models","page":"lmmERP (mass univariate)","title":"Mass Univariate Mixed Models","text":"","category":"section"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"Again we have 4 steps:","category":"page"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"Split data into epochs\nSpecify a formula\nFit a linear model to each time point & channel\nVisualize the results.","category":"page"},{"location":"tutorials/lmm_mu/#1.-Epoching","page":"lmmERP (mass univariate)","title":"1. Epoching","text":"","category":"section"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"data, evts = UnfoldSim.predef_eeg(10; return_epoched = true) # simulate 10 subjects\ndata = reshape(data, 1, size(data, 1), :) # concatenate the data into a long EEG dataset\ntimes = range(0, length = size(data, 2), step = 1 / 100)\ntransform!(evts, :subject => categorical => :subject); # :subject must be categorical, otherwise MixedModels.jl complains\nnothing #hide","category":"page"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"The events dataFrame has an additional column (besides being much taller): subject","category":"page"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"first(evts, 6)","category":"page"},{"location":"tutorials/lmm_mu/#2.-Formula-specification","page":"lmmERP (mass univariate)","title":"2. Formula specification","text":"","category":"section"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"We define the formula. Importantly, we need to specify a random effect. We use zerocorr to speed up the calculation.","category":"page"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"f = @formula 0 ~ 1 + condition * continuous + zerocorr(1 + condition * continuous | subject);\nnothing #hide","category":"page"},{"location":"tutorials/lmm_mu/#3.-Model-fitting","page":"lmmERP (mass univariate)","title":"3. Model fitting","text":"","category":"section"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"We can now run the LinearMixedModel at each time point.","category":"page"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"m = fit(UnfoldModel, f, evts, data, times)\nnothing #hide","category":"page"},{"location":"tutorials/lmm_mu/#4.-Visualization-of-results","page":"lmmERP (mass univariate)","title":"4. Visualization of results","text":"","category":"section"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"Let's start with the fixed effects. We see the condition effects and some residual overlap activity in the fixed effects.","category":"page"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"results = coeftable(m)\n\nres_fixef = results[isnothing.(results.group), :]\nplot_erp(res_fixef)","category":"page"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"And now comes the random effect:","category":"page"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"res_ranef = results[results.group .== :subject, :]\nplot_erp(res_ranef)","category":"page"},{"location":"tutorials/lmm_mu/#Statistics","page":"lmmERP (mass univariate)","title":"Statistics","text":"","category":"section"},{"location":"tutorials/lmm_mu/","page":"lmmERP (mass univariate)","title":"lmmERP (mass univariate)","text":"Check out the LMM p-value tutorial","category":"page"}]
}
