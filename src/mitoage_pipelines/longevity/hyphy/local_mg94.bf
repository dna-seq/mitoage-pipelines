/* MG94xREV with a synonymous rate (alpha) and a non-synonymous rate (beta) of its own on every branch, fitted by
   maximum likelihood; writes alpha and beta of each branch, which are proportional to the branch's dS and dN
   (both are expected substitutions per codon times one constant of the fitted model), so beta/alpha is its dN/dS.
   Written for `mitoage longevity branch-rates`, from HyPhy 2.5's FitModel.bf with the model made local.
   hyphy local_mg94.bf --alignment X.fna --tree X.nwk --code Vertebrate-mtDNA --output X.json */
RequireVersion ("2.5.34");

LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("SelectionAnalyses/modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");
LoadFunctionLibrary("libv3/convenience/regexp.bf");

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

KeywordArgument ("code",      "Which genetic code should be used", "Universal");
KeywordArgument ("alignment", "An in-frame codon alignment");
KeywordArgument ("tree",      "A phylogenetic tree", null, "Please select a tree file for the data:");

lmg.json = {terms.json.input: {}, terms.json.fits: {}, terms.json.timers: {}};

namespace lmg {
    LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
    load_file ({utility.getGlobalValue("terms.prefix"): "lmg",
                utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "selection.io.SelectAllBranches"}});
}

KeywordArgument ("output", "Write the resulting JSON to this file", lmg.codon_data_info [terms.json.json]);
lmg.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

namespace lmg {
    doGTR ("lmg");
}

lmg.results = estimators.FitCodonModel (lmg.filter_names, lmg.trees, "models.codon.MG_REV.ModelDescription",
    lmg.codon_data_info [utility.getGlobalValue("terms.code")],
    {
        terms.run_options.model_type: terms.local,
        terms.run_options.retain_lf_object: TRUE,
        terms.run_options.retain_model_object : TRUE
    },
    lmg.gtr_results);

io.ReportProgressMessageMD ("lmg", "fit", "* " + selection.io.report_fit (lmg.results, 0, lmg.codon_data_info[terms.data.sample_size]));

/* the local parameters of the fitted likelihood function, by name: <tree>.<branch>.alpha and .beta */
GetString (lmg.lf_info, ^(lmg.results[terms.likelihood_function]), -1);
lmg.locals = lmg.lf_info["Local Independent"];
lmg.branches = {};
for (lmg.i = 0; lmg.i < Columns (lmg.locals); lmg.i += 1) {
    lmg.name = lmg.locals[lmg.i];
    lmg.parts = regexp.FindSubexpressions (lmg.name, "^(.+)\\.([^.]+)\\.([^.]+)$");
    lmg.branch = lmg.parts[2];
    lmg.par = lmg.parts[3];
    if (lmg.branches / lmg.branch == FALSE) {
        lmg.branches[lmg.branch] = {};
    }
    (lmg.branches[lmg.branch])[lmg.par] = Eval (lmg.name);
}

lmg.json ["log likelihood"] = lmg.results[terms.fit.log_likelihood];
lmg.json ["parameters"] = lmg.results[terms.parameters];
lmg.json ["globals"] = utility.Map (lmg.results[terms.global], "_value_", '_value_ [terms.fit.MLE]');
lmg.json ["branches"] = lmg.branches;
lmg.json ["branch lengths"] = (lmg.results[terms.branch_length])[0];
io.SpoolJSON (lmg.json, lmg.codon_data_info [terms.json.json]);
