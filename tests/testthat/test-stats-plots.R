test_that("BetweenGroup.test is robust to degenerate and pre-transformed features", {
  set.seed(3)
  x <- data.frame(a=c(rep(0, 10), rep(5, 10)), b=rpois(20, 3), c=rpois(20, 5))
  y <- factor(rep(c("A", "B"), each=10))
  res <- BetweenGroup.test(x, y)
  expect_equal(nrow(res), 3)
  expect_equal(levels(res$Enr), c("Neutral", "A_depleted", "A_enriched"))
  expect_equal(nrow(BetweenGroup.test(x, as.character(y), paired=TRUE)), 3)
  # log fold changes of negative values are NaN (with warnings), but the tests still run
  expect_equal(nrow(suppressWarnings(BetweenGroup.test(data.frame(a=rnorm(20), b=rnorm(20)), y))), 2)
})

test_that("feature-level summaries across datasets", {
  d <- sim_clf_data()
  res <- quiet(rf_clf.by_datasets(d$df, d$md, "study", "disease", "Case", ntree=50, n_cores=2))
  consistency <- corr_datasets_by_imps(res, top_n=2)
  expect_equal(dim(consistency$corr_mat), c(3, 3))
  expect_true(all(consistency$feature_consistency$marker_class %in% c("universal", "dataset_specific", "non_marker")))
  heatmap <- plot_logfc_heatmap(res, min_sig_datasets=0)
  expect_equal(dim(heatmap$logfc_mat), c(5, 3))
  expect_s3_class(plot_clf_ROC(res), "ggplot")
  x <- d$df[1:60, ]
  y <- d$md$disease[1:60]
  model <- quiet(rf.cross.validation(x, y, nfolds=3, ntree=50))
  markers <- id_robust_markers(model, BetweenGroup.test(x, y), top_n=2)
  expect_equal(sum(markers$rf_selected), 2)
  oob <- quiet(rf.out.of.bag(x, y, ntree=50))
  top <- plot_topN_imp_scores(oob, topN=2)
  expect_setequal(as.character(top$plot$data$Feature_ID), names(sort(oob$importances, decreasing=TRUE))[1:2])
})

test_that("permutation tests handle classification and regression", {
  d <- sim_clf_data()
  x <- d$df[1:60, ]
  y <- d$md$disease[1:60]
  model <- quiet(rf.cross.validation(x, y, nfolds=3, ntree=50))
  expect_lt(plot_test_perf_VS_rand(model, x, y, metric="AUROC", permutation=50)$emp_p_value, 0.05)
  reg <- quiet(rf.cross.validation(x, d$md$age[1:60], nfolds=3, ntree=50))
  p_value <- plot_test_perf_VS_rand(reg, x, d$md$age[1:60], metric="MAE", permutation=50)$emp_p_value
  expect_true(p_value > 0 && p_value <= 1)
})

test_that("the R Markdown report renders", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available())
  d <- sim_clf_data()
  res <- quiet(rf_clf.by_datasets(d$df, d$md, "study", "disease", "Case", ntree=50, n_cores=2))
  out <- quiet(crossRanger_report(Within_study=res, outdir=tempfile("report")))
  expect_true(file.exists(out))
})
