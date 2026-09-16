test_that("OOB votes and probabilities are consistent with a per-sample count", {
  d <- sim_clf_data()
  x <- d$df[1:60, ]
  y <- d$md$disease[1:60]
  m <- quiet(rf.out.of.bag(x, y, ntree=50))
  votes <- get.oob.votes.from.forest(m$rf.model, x)
  inbag <- do.call(cbind, m$rf.model$inbag.counts)
  preds <- predict(m$rf.model, x, predict.all=TRUE)$predictions
  naive <- t(sapply(seq_len(nrow(x)), function(i) table(factor(preds[i, inbag[i, ]==0], levels=1:2))))
  expect_equal(unname(votes), unname(naive))
  expect_equal(unname(rowSums(m$probabilities)), rep(1, 60))
})

test_that("feature names are preserved for cross-dataset prediction", {
  d <- sim_clf_data()
  x <- d$df[1:60, ]
  colnames(x) <- paste0("k__Bacteria;g__", LETTERS[1:5])
  y <- d$md$disease[1:60]
  m <- quiet(rf.out.of.bag(x, y, ntree=50))
  expect_identical(names(m$importances), colnames(x))
  expect_equal(nrow(get.predict.probability.from.forest(m$rf.model, x)), 60)
  m_p <- quiet(rf.out.of.bag(x, y, ntree=20, imp_pvalues=TRUE))
  expect_identical(rownames(m_p$importances), colnames(x))
})

test_that("cross-validation supports leave-one-out and group k-fold", {
  d <- sim_clf_data()
  x <- d$df[1:60, ]
  y <- d$md$disease[1:60]
  loo <- quiet(rf.cross.validation(x[1:20, ], factor(rep(c("a", "b"), 10)), nfolds=-1, ntree=20))
  expect_length(loo$errs, 20)
  subject <- factor(rep(1:20, each=3))
  cv <- quiet(rf.cross.validation(x, y, nfolds=5, groups=subject, ntree=20))
  expect_true(all(tapply(cv$folds, subject, function(f) length(unique(f))) == 1))
  expect_setequal(unique(group.folds(subject, nfolds=4)), 1:4)
})

test_that("AUROC and log.mat handle ties, class order and zeros", {
  y <- factor(c("A", "A", "B", "B"))
  expect_equal(get.auroc(c(0.9, 0.8, 0.2, 0.1), y, "A"), 1)
  expect_equal(get.auroc(c(0.9, 0.8, 0.2, 0.1), y, "B"), 0)
  expect_equal(get.auroc(c(0.5, 0.5, 0.5, 0.5), y, "A"), 0.5)
  expect_equal(get.auroc(c(0.9, 0.2, 0.9, 0.1), y, "A"), 0.625)
  expect_equal(log.mat(data.frame(a=c(0, 1), b=c(2, 4)))$a, log2(c(0.5, 1.5)))
})

test_that("regression performance metrics are correct", {
  y <- 1:60
  set.seed(1)
  pred <- y + rnorm(60, 0, 5)
  perf <- get.reg.performance(pred, y, n_features=5)
  expect_equal(perf$MAE, mean(abs(y - pred)))
  expect_equal(perf$Adj_R_squared, 1 - (1 - perf$R_squared) * 59 / 54)
  expect_true(is.na(get.reg.performance(pred, y, n_features=100)$Adj_R_squared))
})

test_that("recursive feature elimination selects a parsimonious feature set", {
  set.seed(7)
  pv <- gtools::rdirichlet(2, sample(100))
  x <- data.frame(rbind(t(rmultinom(30, 3000, pv[1, ])), t(rmultinom(30, 3000, pv[2, ]))))
  y <- factor(rep(c("A", "C"), each=30))
  rfe <- rf_clf.rfe(x, y, nfolds=5, ntree=50)
  expect_s3_class(rfe, "rf_clf.rfe")
  expect_equal(rfe$top_n_perf$n_features, c(2, 4, 8, 16, 32, 64, 100))
  expect_lte(rfe$optimal_n_features, rfe$best_n_features)
  expect_length(rfe$selected_features, rfe$optimal_n_features)
  expect_s3_class(plot_clf_feature_selection(rfe)$plot, "ggplot")
  expect_equal(nrow(plot_topN_imp_scores(rfe, topN=3, n_features=16)$plot$data), 3)
  reg <- rf_reg.rfe(x, 1:60, nfolds=5, ntree=50)
  expect_s3_class(plot_reg_feature_selection(reg)$plot, "ggplot")
})
