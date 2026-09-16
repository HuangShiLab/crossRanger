test_that("rf_clf.by_datasets returns models and a performance summary", {
  d <- sim_clf_data()
  res <- quiet(rf_clf.by_datasets(d$df, d$md, s_category="study", c_category="disease",
                                  positive_class="Case", ntree=50, n_cores=2))
  expect_s3_class(res, "rf_clf.by_datasets")
  expect_equal(res$perf_summ$Dataset, c("S1", "S2", "S3"))
  expect_named(res$rf_model_list, c("S1", "S2", "S3"))
  expect_true(all(res$perf_summ$AUROC > 0.5))
})

test_that("rf_clf.by_datasets accepts a list of datasets and group k-fold CV", {
  d <- sim_clf_data()
  res <- quiet(rf_clf.by_datasets(split(d$df, d$md$study), split(d$md, d$md$study), c_category="disease",
                                  positive_class="Case", nfolds=3, cv_type="group", g_category="subject",
                                  ntree=50, n_cores=2))
  expect_equal(res$datasets, c("S1", "S2", "S3"))
  expect_true(all(res$perf_summ$Validation_type == "group_kfold_CV"))
  expect_error(rf_clf.by_datasets(d$df, d$md, "study", "disease", cv_type="group"), "g_category")
})

test_that("cross-application and LODO work for classification", {
  d <- sim_clf_data()
  res <- quiet(rf_clf.by_datasets(d$df, d$md, "study", "disease", "Case", nfolds=3, ntree=50, n_cores=2))
  cross_rf <- quiet(rf_clf.cross_appl(res$rf_model_list, res$x_list, res$y_list, positive_class="Case"))
  expect_equal(nrow(cross_rf$perf_summ), 9)
  expect_length(cross_rf$prediction_list, 9)
  expect_error(rf_clf.cross_appl(res$rf_model_list, res$x_list[1:2], res$y_list), "identical")
  expect_s3_class(plot_cross_appl(cross_rf), "ggplot")
  lodo <- quiet(rf_clf.lodo(d$df, d$md, "study", "disease", "Case", ntree=50, n_cores=2))
  expect_equal(lodo$perf_summ$Test_data, c("S1", "S2", "S3"))
  expect_s3_class(plot_cross_appl(lodo), "ggplot")
})

test_that("regression across datasets reports correct metrics", {
  d <- sim_clf_data()
  for(k in c(1, 3)){
    reg <- quiet(rf_reg.by_datasets(d$df, d$md, "study", "age", nfolds=k, ntree=50, n_cores=2))
    expect_null(dim(reg$feature_imps_list[[1]]))
    expect_length(reg$rf_MAE, 3)
  }
  cross_reg <- quiet(rf_reg.cross_appl(reg, reg$x_list, reg$y_list))
  pred <- cross_reg$predicted[[2]]
  perf <- get.reg.performance(pred$pred_y, pred$test_y)
  expect_equal(cross_reg$perf_summ$R_squared[2], perf$R_squared)
  expect_equal(cross_reg$perf_summ$Spearman_rho[2], perf$Spearman_rho)
  expect_equal(nrow(quiet(rf_reg.lodo(d$df, d$md, "study", "age", ntree=50, n_cores=2))$perf_summ), 3)
})

test_that("the summary pipeline runs", {
  d <- sim_clf_data()
  summ <- quiet(rf_clf.by_datasets.summ(d$df, d$md, "study", "disease", "Case", nfolds=3, ntree=50, n_cores=2))
  expect_equal(nrow(summ$perf_summ), 3)
})
