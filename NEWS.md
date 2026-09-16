# crossRanger 0.2.0

## New features

* `rf_clf.by_datasets()` and `rf_reg.by_datasets()`
  * accept a named list of feature tables (`df`) and metadata, harmonized to the shared features;
  * match samples by IDs (rownames) when both tables have them;
  * support group k-fold (`cv_type = "group"`) and leave-one-group-out (`cv_type = "logo"`) CV with `g_category`;
  * return a within-study performance summary `perf_summ`, and take `n_cores` and `...` (passed to `ranger`).
* `rf_clf.lodo()` and `rf_reg.lodo()`: leave-one-dataset-out validation.
* `rf_clf.rfe()` and `rf_reg.rfe()`: recursive feature elimination with a tolerance rule
  (the smallest feature set within 1% of the best performance). `plot_clf_feature_selection()` and
  `plot_reg_feature_selection()` accept their results, and `plot_topN_imp_scores()` plots any RFE step.
* `group.folds()` and `rf.cross.validation(groups = )` for group k-fold CV.
* `harmonize_features()` for multi-dataset feature harmonization.
* `plot_logfc_heatmap()`: log2 fold changes of features across datasets.
* `id_robust_markers()`: intersection of RF-selected and statistically significant features.
* `corr_datasets_by_imps()` accepts result objects and classifies features into universal and dataset-specific markers (`top_n`).
* `plot_clf_ROC()` / `plot_clf_PRC()` draw the curves of all datasets of an `rf_clf.by_datasets` object.
* `plot_perf_VS_rand()` and `plot_test_perf_VS_rand()` support classification metrics.
* `crossRanger_report()` renders an R Markdown report (HTML, PDF or Word); `rf_clf.by_datasets.summ()`
  and `rf_clf.comps.summ()` gain `report_format`.
* `rf_clf.cross_appl()` returns `prediction_list`; models from cross-validation are applied as an
  ensemble of all fold models, and test features are matched to the model features by name.
* Example dataset `crossRanger_sim` and the vignette `crossRanger--intro`.
* The package can now be installed on Windows: the parallel backend is doParallel instead of doMC,
  which is a unix-only package and therefore blocked the installation of crossRanger on Windows.
  doParallel forks the R session on Unix (the same mechanism, performance and results as doMC)
  and starts a PSOCK cluster on Windows.
* A parallel backend registered by the user (e.g., doFuture, or a cluster of their own) is now kept
  instead of being silently replaced, so the way of parallelization can be chosen by the user.
* `rf.out.of.bag()`, `plot_perf_VS_rand()` and `plot_test_perf_VS_rand()` gain a `seed` argument
  (default 123, which reproduces the results of earlier versions; NULL uses the current random number
  stream). These functions no longer leave the random number generator of the user's R session reset:
  its state is restored when they return.

## Bug fixes

* Models trained by `rf.out.of.bag()` on features with non-syntactic names (e.g., `k__Bacteria;g__Prevotella`)
  failed to predict other datasets.
* `rf_reg.cross_appl()` swapped the `R_squared` and `Spearman_rho` columns.
* `rf_reg.by_datasets()` failed with `nfolds = 1`, returned per-fold importance matrices with `nfolds = 3`,
  mislabeled results with `rf_imp_pvalues = TRUE`, and returned per-fold metrics that did not match the datasets.
* `rf_clf.by_datasets.summ()` passed arguments to `rf_clf.by_datasets()` by position into the wrong parameters and crashed;
  it and `rf_clf.comps.summ()` ignored `nfolds`, `p.adj.method` and `q_cutoff`.
* `plot_cross_appl()` failed for classification results with the default metric, and reordered datasets alphabetically.
* `plot_topN_imp_scores()` failed for `rf.out.of.bag` models and did not select the top-ranked features.
* `get.reg.performance()` computed the adjusted R squared with a squared R squared.
* `rf.cross.validation()` failed for leave-one-out CV, and computed importance p values of all folds from the last fold's data.
* `plot_perf_VS_rand()` ignored `permutation`, swapped observed and predicted values, and computed p values in the
  wrong direction for metrics where higher is better; `plot_test_perf_VS_rand()` re-predicted the data at each permutation.
* `plot_reg_feature_selection()` ignored `nfolds`.
* `rf_clf.cross_appl()` / `rf_reg.cross_appl()` did not check the lengths of the input lists.
* `rf_clf.comps()` failed when AUROC could not be computed.
* `BetweenGroup.test()` failed for non-factor `y`, for features constant within both groups, and for negative values;
  paired tests now use the two-sample interface.
* Plotting functions no longer write files when `outdir = NULL`; `boxplot_rel_predicted_train_vs_test()` saved the wrong plot
  and `plot_rel_predicted()` failed.
* The functions that write result tables used `sink()` without protecting it with `on.exit()`, so that an error
  while writing left the output of the user's R session redirected. The tables are now written through a file
  connection, with identical content.
* `R/data_trimming_util.R` no longer installs and attaches packages when the package is loaded, and no longer needs dplyr;
  `normalize_NA_in_metadata()`, `discard_uninfo_columns_in_metadata()`, `filter_features_allzero()` and `check_metadata(more_missing_values=)` are fixed.
* The number of cores is at least one on machines with four cores or fewer.
* `plot_clf_pROC()` failed whenever `outdir` was set, drew on the active device even with `outdir = NULL`,
  and let pROC choose the curve direction automatically (so the AUC of a poor model was flipped above 0.5).
  The positive class is now the case level with `direction = "<"`.
* `log.mat()` was registered only as an S3 method of `log()` and not exported, and failed for data.frames with zeros.
* `rf.out.of.bag(imp_pvalues = TRUE)` works with non-syntactic feature names.
* Broken examples of `plot_clf_ROC()`, `plot_clf_pROC()`, `rf_clf.pairwise()`, `mttest()`, `get.mislabel.scores()`,
  `plot_train_vs_test()` and `plot_rel_predicted()` are fixed.

## Other changes

* Vote counting for OOB and test probabilities is vectorized.
* `get.auroc()` uses the Mann-Whitney formula (identical values, about 16 times faster than ROCR), so ROCR is no longer a dependency;
  `get.auprc()` no longer computes the unused maximum/minimum curves.
* Models of cross-validation folds no longer keep in-bag counts (about 45% less memory; identical results).
* Within parallel workers, ranger threads are shared among the workers instead of each worker using all cores (identical results).
* `rf.cross.validation()` prints the per-fold progress only when `verbose = TRUE`.
* `mttest()` runs its features sequentially when no parallel backend is registered (as inside the workers
  of another parallel loop), which avoids nested parallelism and the warning of foreach.
* The implicit cluster that doParallel starts on Windows is stopped when the package is unloaded.
* The Title and Description fields follow the CRAN requirements (title case, and not starting with
  the package name), and the slowest examples of `rf.cross.validation()` are wrapped in `\donttest{}`.
