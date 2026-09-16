#' @importFrom stats cor wilcox.test residuals smooth.spline fitted median
#' @importFrom plyr ldply
#' @importFrom rlang .data
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob
#' @importFrom reshape2 melt
#' @importFrom utils write.table
#' @importFrom gtools rdirichlet

.higher_is_better <- function(metric) !metric %in% c("MSE", "RMSE", "nRMSE", "MAE", "MAPE", "MASE")

# Performance of a prediction: `prediction` is a numeric vector for regression, or a list/model
# with `predicted` classes and class `probabilities` for classification (y is a factor).
.perf_metric <- function(y, prediction, metric, positive_class=NA, n_features=NA){
  if(is.factor(y)){
    positive_class <- ifelse(is.na(positive_class), levels(y)[1], positive_class)
    .clf_perf_row(prediction$predicted, y, prediction$probabilities, positive_class)[[metric]]
  }else{
    get.reg.performance(prediction, y, n_features)[[metric]]
  }
}

.emp_p_value <- function(perf_value, rand_perf_values, higher_better){
  rand_perf_values <- rand_perf_values[!is.na(rand_perf_values)]
  n_extreme <- if(higher_better) sum(rand_perf_values >= perf_value) else sum(rand_perf_values <= perf_value)
  (n_extreme + 1)/(length(rand_perf_values) + 1)
}

.plot_perf_vs_rand <- function(perf_value, rand_perf_values, emp_p_value, metric){
  perf_values<-data.frame(perf_value, rand_perf_values)
  label = paste(metric, ": ", as.character(round(perf_value, 2)), "\np-value = ", signif(emp_p_value, 2), sep="")
  p<-ggplot(perf_values, aes(x=.data$rand_perf_values)) + geom_histogram(alpha=0.5, bins=30) +
    xlab(metric)+
    ylab("count")+
    geom_vline(xintercept = perf_value, color="red") +
    annotate(geom="text", x=perf_value, y=Inf, label=label, color="red", vjust=2, hjust=0)+theme_bw()
  list(perf_values=perf_values, plot=p)
}

#' @title plot_obs_VS_pred
#' @description Plot a scatterplot of observed and predicted values from a ranger model.
#' @param y The numeric values for labeling data.
#' @param predicted_y The predicted values for y.
#' @param metric A regression performance metric (i.e., MAE, RMSE, MSE, R_squared, Adj_R_squared, or Separman_rho) showed in the scatter plot.
#' @param SampleIDs The sample ids in the data table.
#' @param prefix The prefix for the dataset in the training or testing.
#' @param target_field A string indicating the target field of metadata for regression.
#' @param span Controls the amount of smoothing for the default loess smoother in the scatterplot.
#' Smaller numbers produce wigglier lines, larger numbers produce smoother lines.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<- 1:60
#' rf_model<-rf.out.of.bag(x, y)
#' plot_obs_VS_pred(y, predicted_y=rf_model$predicted, prefix="train", target_field="age", outdir=NULL)
#' @author Shi Huang
#' @export
plot_obs_VS_pred <- function(y, predicted_y, SampleIDs=NULL, prefix="train", target_field="value", metric="MAE", span=1, outdir=NULL){
  df<-data.frame(y, predicted_y)
  if(nrow(df)==length(SampleIDs)){
    df<-data.frame(SampleIDs, df)
  }else if(length(SampleIDs)>0 & nrow(df)!=length(SampleIDs)){
    stop("Please make sure that sample IDs match with y or predicted y.")
  }
  perf_value<-get.reg.performance(predicted_y, y)[[metric]]
  label = paste(metric, ": ", as.character(round(perf_value, 2)), sep="")
  p<-ggplot(df, aes(x=.data$y, y=.data$predicted_y))+
    ylab(paste("Predicted ",target_field,sep=""))+
    xlab(paste("Observed ",target_field,sep=""))+
    geom_point(alpha=0.1)+
    geom_smooth(method="loess", formula=y ~ x, span=span)+
    annotate(geom="text", x=Inf, y=Inf, label=label, color="grey60", vjust=2, hjust=2)+
    theme_bw()
  if(!is.null(outdir)){
    ggsave(filename=paste(outdir, prefix, ".", target_field, ".obs_vs_pred.scatterplot.pdf",sep=""), plot=p, height=4, width=4)
    .write_tsv(df, paste(outdir, prefix, ".", target_field, ".obs_vs_pred.results.xls",sep=""))
  }
  p
}
#' @title plot_residuals
#' @description Plot the residuals of observed and predicted values from a ranger model.
#' @param y The numeric values for labeling data.
#' @param predicted_y The predicted values for y.
#' @param SampleIDs The sample ids in the data table.
#' @param prefix The prefix for the dataset in the training or testing.
#' @param target_field A string indicating the target field of the metadata for regression.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<- 1:60
#' rf_model<-rf.out.of.bag(x, y)
#' plot_residuals(y, predicted_y=rf_model$predicted, prefix="train", target_field="age")
#' @author Shi Huang
#' @export
plot_residuals <- function(y, predicted_y, SampleIDs=NULL, prefix="train", target_field="value", outdir=NULL){
  df<-data.frame(y, predicted_y, rsdl=y-predicted_y)
  if(nrow(df)==length(SampleIDs)){
    df<-data.frame(SampleIDs, df)
  }else if(length(SampleIDs)>0 & nrow(df)!=length(SampleIDs)){
    stop("Please make sure that sample IDs match with y or predicted y.")
  }
  p<-ggplot(df, aes(x=.data$y, y=.data$rsdl))+
    ylab(paste("Residuals of prediceted ",target_field,sep=""))+
    xlab(paste("Observed ",target_field,sep=""))+
    geom_point(alpha=0.1)+
    geom_hline(yintercept=0)+
    theme_bw()
  if(!is.null(outdir)){
  ggsave(filename=paste(outdir, prefix, ".", target_field, ".obs_vs_residuals_of_pred.scatterplot.pdf",sep=""), plot=p, height=4, width=4)
  .write_tsv(df, paste(outdir, prefix, ".", target_field, ".obs_vs_pred.results.xls",sep=""))
  }
  invisible(p)
}

#' @title plot_perf_VS_rand
#' @description This outputs a histogram and an empirical p-value showing if the performance of a real
#' regression or classification model is significantly better than null models built on permuted labels.
#' @param x The train data.
#' @param y The numeric labeling data (regression) or a factor (classification).
#' @param nfolds The number of folds in the cross validation. If nfolds > length(y)
#' or nfolds==-1, uses leave-one-out cross-validation. If nfolds was a factor, it means customized folds
#' (e.g., leave-one-group-out cv) were set for CV.
#' @param predicted_y The predicted values for y (regression), or the rf classification model
#' (\code{rf.out.of.bag} or \code{rf.cross.validation}) whose predictions are evaluated (classification).
#' @param n_features The number of features in the training data.
#' @param prefix The prefix for the dataset in the training or testing.
#' @param target_field A string indicating the target field of the metadata.
#' @param metric The performance metric applied: MAE, RMSE, MSE, MAPE, R_squared, Adj_R_squared or Spearman_rho for regression;
#' AUROC, AUPRC, Accuracy, Kappa, F1 or Balanced_Accuracy for classification.
#' @param permutation The permutation times for a random guess of performance.
#' @param outdir The output directory.
#' @param positive_class A class of y for classification.
#' @param n_cores The number of cores for running the permutations in parallel.
#' @param seed The random seed used for reproducible models. The default (123) reproduces the
#' results of earlier versions; NULL uses the current state of the random number generator.
#' The state of the user's R session is restored when the function returns.
#' @return A list including the empirical p value, the observed and permuted performance, and the histogram.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<- 1:60
#' rf_model<-rf.out.of.bag(x, y)
#' p<-plot_perf_VS_rand(x=x, y=y, predicted_y=rf_model$predicted, prefix="train", nfolds=5,
#' permutation=20, metric="MAE", target_field="age", n_features=5)
#' p$emp_p_value
#' y_clf <- factor(rep(c("A", "B"), each=30))
#' clf_model <- rf.cross.validation(x, y_clf, nfolds=5)
#' plot_perf_VS_rand(x=x, y=y_clf, predicted_y=clf_model, nfolds=5,
#'                   permutation=20, metric="AUROC")$plot
#' @author Shi Huang
#' @export
plot_perf_VS_rand<-function(x, y, predicted_y, prefix="train", target_field="value", nfolds=5,
                            metric="MAE", permutation=100, n_features=NA, outdir=NULL, positive_class=NA,
                            n_cores=1, seed=123){
  # save the state at entry so that it can be restored, but set the seed only just before the
  # permutations, which keeps the permuted labels identical to those of earlier versions
  seed_state <- if(is.null(seed)) NULL else .save_seed_state()
  on.exit(.restore_seed(seed_state), add=TRUE)
  if(is.factor(y)){
    if(!metric %in% .clf_metric_names) stop("metric should be one of: ", paste(.clf_metric_names, collapse=", "))
    if(!is.list(predicted_y) || is.null(predicted_y$probabilities))
      stop("For classification, predicted_y should be an rf.out.of.bag/rf.cross.validation object (or a list with 'predicted' and 'probabilities').")
    observed <- list(predicted=predicted_y$predicted, probabilities=predicted_y$probabilities)
  }else{
    observed <- predicted_y
  }
  perf_value<-.perf_metric(y, observed, metric, positive_class, n_features)
  if(!is.null(seed)) set.seed(seed)
  rand_y_list <- lapply(seq_len(permutation), function(k) sample(y, replace = FALSE))
  shuffle_y_perf <- function(rand_y){
    rand_rf <- .quietly(rf.cross.validation(x, rand_y, nfolds=nfolds))
    .perf_metric(rand_y, if(is.factor(y)) rand_rf else rand_rf$predicted, metric, positive_class, n_features)
  }
  if(n_cores > 1){
    n_workers <- .register_cores(n_cores)
    if(isTRUE(attr(n_workers, "own"))) on.exit(.stop_implicit_cluster(), add=TRUE)
    rand_perf_values <- unlist(foreach(i=seq_along(rand_y_list)) %dopar% shuffle_y_perf(rand_y_list[[i]]))
  }else{
    rand_perf_values <- vapply(rand_y_list, shuffle_y_perf, numeric(1))
  }
  emp_p_value<-.emp_p_value(perf_value, rand_perf_values, .higher_is_better(metric))
  out <- .plot_perf_vs_rand(perf_value, rand_perf_values, emp_p_value, metric)
  if(!is.null(outdir)){
    ggsave(filename=paste(outdir, prefix, ".", target_field, ".", metric,
                          "_vs_rand.histogram.pdf",sep=""), plot=out$plot, height=4, width=4)
  }
  res <- list()
  res$emp_p_value <- emp_p_value
  res$perf_values <- out$perf_values
  res$plot <- out$plot
  res
}

#' @title plot_train_vs_test
#' @description Plot the observed and predicted values of both training and testing data from a ranger model.
#' @param train_y The numeric labels for training data.
#' @param predicted_train_y The predicted values for training data.
#' @param train_SampleIDs The sample ids in the train data.
#' @param test_y The numeric labels for testing data.
#' @param predicted_test_y The predicted values for test data.
#' @param test_SampleIDs The sample ids in the test data.
#' @param train_prefix The prefix for the dataset in the training data.
#' @param test_prefix The prefix for the dataset in the testing data.
#' @param train_target_field A string indicating the target field of the training metadata for regression.
#' @param test_target_field A string indicating the target field of the testing metadata for regression.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<- 1:60
#' newx <- data.frame(rbind(t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' newy <- 31:60
#' rf_model<-rf.out.of.bag(x, y)
#' predicted_newy <- predict(rf_model$rf.model, newx)$predictions
#' plot_train_vs_test(train_y=y, predicted_train_y=rf_model$predicted,
#'                    test_y=newy, predicted_test_y=predicted_newy,
#'                    train_target_field="age", test_target_field="age")
#' @author Shi Huang
#' @export
plot_train_vs_test<-function(train_y, predicted_train_y, test_y, predicted_test_y, train_SampleIDs=NULL, test_SampleIDs=NULL,
                             train_prefix="train", test_prefix="test", train_target_field="value", test_target_field="value", outdir=NULL){
  train.pred<-data.frame(value=train_y,predicted_value=predicted_train_y)
  if(nrow(train.pred)==length(train_SampleIDs)){
    train.pred<-data.frame(SampleIDs=train_SampleIDs, train.pred)
  }else if(length(train_SampleIDs)>0 & nrow(train.pred)!=length(train_SampleIDs)){
    stop("Please make sure that sample IDs match with y or predicted y in the train data.")
  }
  test.pred<-data.frame(value=test_y,predicted_value=predicted_test_y)
  if(nrow(test.pred)==length(test_SampleIDs)){
    test.pred<-data.frame(SampleIDs=test_SampleIDs, test.pred)
  }else if(length(test_SampleIDs)>0 & nrow(test.pred)!=length(test_SampleIDs)){
    stop("Please make sure that sample IDs match with y or predicted y in the test data.")
  }
  data_name<-c(rep(train_prefix,nrow(train.pred)),rep(test_prefix,nrow(test.pred)))
  pred<-data.frame(data=data_name,rbind(train.pred,test.pred))
  pred$data<-factor(pred$data,levels=c(train_prefix,test_prefix),ordered=TRUE)
  p<-ggplot(pred,aes(x=.data$value,y=.data$predicted_value))+
    ylab(paste("Predicted ",train_target_field,sep=""))+
    xlab(paste("Observed ",train_target_field,sep=""))+
    geom_point(aes(color=.data$data), alpha=0.1)+
    geom_smooth(aes(color=.data$data), method="loess", formula=y ~ x, span=1)+
    theme_bw() +
    facet_wrap(~data)+
    theme(legend.position="none")
  if(!is.null(outdir)){
  ggsave(filename=paste(outdir, train_prefix,"-",test_prefix,".",test_target_field, ".train_test_ggplot.pdf",sep=""),plot=p, height=3, width=6)
  .write_tsv(pred, paste(outdir, train_prefix,"-",test_prefix,".",test_target_field, ".train_test_results.xls",sep=""))
  }
  invisible(p)
}


#' @title plot_test_perf_VS_rand
#' @description This outputs a histogram and an empirical p-value showing if the performance of a trained
#' regression or classification model on new (e.g., independent) data is significantly better than random guesses,
#' i.e., the performance of the same predictions against permuted labels.
#' @param rf_model A trained rf model object, which should be generated from \code{rf.cross.validation} or \code{rf.out.of.bag}.
#' @param newy The data label of new data: numeric (regression) or a factor (classification).
#' @param newx A data.matrix or data.frame with the new data for rf model testing.
#' @param n_features The number of features in the training data.
#' @param prefix The prefix for the dataset in the training or testing.
#' @param target_field A string indicating the target field of the metadata for machine-learning analysis.
#' @param metric The performance metric applied: MAE, RMSE, MSE, MAPE, R_squared, Adj_R_squared or Spearman_rho for regression;
#' AUROC, AUPRC, Accuracy, Kappa, F1 or Balanced_Accuracy for classification.
#' @param permutation The permutation times for a random guess of performance.
#' @param outdir The output directory.
#' @param positive_class A class of newy for classification.
#' @param seed The random seed used for reproducible models. The default (123) reproduces the
#' results of earlier versions; NULL uses the current state of the random number generator.
#' The state of the user's R session is restored when the function returns.
#' @return A list including the empirical p value, the observed and permuted performance, and the histogram.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<- 1:60
#'
#' newx <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' newy<- 4:33
#' rf_model<-rf.out.of.bag(x, y)
#' p_test<-plot_test_perf_VS_rand(rf_model, newx, newy, permutation=100,
#'                                metric="MAE", target_field="age", n_features=5)
#' p_test$emp_p_value
#' @author Shi Huang
#' @export
plot_test_perf_VS_rand<-function(rf_model, newx, newy, prefix="test", target_field="",
                            metric="MAE", permutation=1000, n_features=NA, outdir=NULL, positive_class=NA,
                            seed=123){
  seed_state <- if(is.null(seed)) NULL else .save_seed_state()
  on.exit(.restore_seed(seed_state), add=TRUE)
  if(is.factor(newy) || is.character(newy)){
    if(!metric %in% .clf_metric_names) stop("metric should be one of: ", paste(.clf_metric_names, collapse=", "))
    prediction <- .predict_clf(rf_model, newx)
    newy <- factor(as.character(newy), levels=levels(prediction$predicted))
    if(any(is.na(newy))) stop("newy includes classes that are absent from the training data.")
  }else{
    prediction <- .predict_reg(rf_model, newx)
  }
  perf_value<-.perf_metric(newy, prediction, metric, positive_class, n_features)
  # the predictions are fixed, and only the labels are permuted
  if(!is.null(seed)) set.seed(seed)
  rand_perf_values <- vapply(seq_len(permutation), function(k)
    .perf_metric(sample(newy, replace = FALSE), prediction, metric, positive_class, n_features), numeric(1))
  emp_p_value<-.emp_p_value(perf_value, rand_perf_values, .higher_is_better(metric))
  out <- .plot_perf_vs_rand(perf_value, rand_perf_values, emp_p_value, metric)
  if(!is.null(outdir)){
    ggsave(filename=paste(outdir, prefix, ".", target_field, ".", metric,
                          "_vs_rand.histogram.pdf",sep=""), plot=out$plot, height=4, width=4)
  }
  res <- list()
  res$emp_p_value <- emp_p_value
  res$perf_values <- out$perf_values
  res$plot <- out$plot
  res
}


#' @title plot_reg_feature_selection
#' @description Plot the regression performance against the reduced number of features used in the modeling.
#' The dashed line marks the optimal parsimonious feature set: the smallest set whose performance is within
#' \code{tolerance} of the best performance.
#' @param x The data frame or data matrix for model training, or the result of \code{rf_reg.rfe}.
#' @param y The numeric values for labeling data.
#' @param nfolds The number of folds in the cross-validation for each feature set.
#' @param unit The unit of numeric metadata variables that can be printed in the output figure.
#' @param rf_reg_model An optional rf regression model from \code{rf.out.of.bag} or \code{rf.cross.validation} used to rank the features.
#' @param metric The regression performance metric applied.
#' This must be one of "MAE", "RMSE", "MSE", "MAPE", "Spearman_rho", "R_squared".
#' @param outdir The output directory.
#' @param tolerance The relative tolerance for choosing the optimal parsimonious feature set.
#' @param recursive A boolean value indicating if the feature ranking is re-computed at each step (see \code{rf_reg.rfe}).
#' @param ntree The number of trees.
#' @return A list including the performance table (\code{top_n_perf}), \code{best_n_features},
#' \code{optimal_n_features}, \code{selected_features}, the models (\code{top_n_rf}), the plot and the \code{rf_reg.rfe} object.
#' @examples
#' set.seed(123)
#' n_features <- 100
#' prob_vec <- gtools::rdirichlet(5, sample(n_features))
#' x <- data.frame(rbind(t(rmultinom(7, 7*n_features, prob_vec[1, ])),
#'             t(rmultinom(8, 8*n_features, prob_vec[2, ])),
#'             t(rmultinom(15, 15*n_features, prob_vec[3, ])),
#'             t(rmultinom(15, 15*n_features, prob_vec[4, ])),
#'             t(rmultinom(15, 15*n_features, prob_vec[5, ]))))
#' y<- 1:60
#' rf_reg_model<-rf.cross.validation(x, y, nfolds=5)
#' fs_summ <- plot_reg_feature_selection(x, y, rf_reg_model, metric="MAE", outdir=NULL)
#' fs_summ$plot
#' @author Shi Huang
#' @export
plot_reg_feature_selection <- function(x, y, rf_reg_model=NULL, nfolds=5, metric="MAE",
                                       unit=NA, outdir=NULL, tolerance=0.01, recursive=FALSE, ntree=500){
  if(inherits(x, "rf_reg.rfe")){
    rfe <- x
    if(missing(metric)) metric <- rfe$metric
  }else{
    rfe <- rf_reg.rfe(x, y, nfolds=nfolds, metric=metric, tolerance=tolerance, recursive=recursive,
                      ntree=ntree, rf_model=rf_reg_model)
  }
  top_n_perf <- rfe$top_n_perf
  if(!metric %in% colnames(top_n_perf)) stop("metric should be one of: ", paste(colnames(top_n_perf)[-1], collapse=", "))
  selection <- .rfe_select(top_n_perf, metric, higher_better=.higher_is_better(metric), tolerance=rfe$tolerance)
  y_label <- ifelse(is.na(unit), metric, paste(metric, " (", unit ,")", sep=""))
  p <- .plot_feature_selection(top_n_perf, metric, selection, y_label=y_label)
  if(!is.null(outdir)){
  ggsave(filename=paste(outdir,"rf__",metric,"__top_rankings.scatterplot.pdf",sep=""), plot=p, width=5, height=4)
  }
  res <- list()
  res$top_n_perf <- top_n_perf
  res$best_n_features <- selection$best_n_features
  res$optimal_n_features <- selection$optimal_n_features
  res$selected_features <- rfe$feature_sets[[as.character(selection$optimal_n_features)]]
  res$top_n_rf <- rfe$top_n_rf
  res$plot <- p
  res$rfe <- rfe
  res
}

#' @title calc_rel_predicted
#' @description Calculate the relative predicted values to the spline fit.
#' @param train_y The numeric labels for training data.
#' @param predicted_train_y The predicted values for training data.
#' @param train_SampleIDs The sample ids in the train data.
#' @param test_y The numeric labels for testing data.
#' @param predicted_test_y The predicted values for test data.
#' @param test_SampleIDs The sample ids in the test data.
#' @param train_prefix The prefix for the dataset in the training data.
#' @param test_prefix The prefix for the dataset in the testing data.
#' @param train_target_field A string indicating the target field of the training metadata for regression.
#' @param test_target_field A string indicating the target field of the testing metadata for regression.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' train_x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' train_y<- 1:60
#' test_x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209)))))
#' test_y<- 1:45
#' train_rf_model<-rf.out.of.bag(train_x, train_y)
#' predicted_test_y<-predict(train_rf_model$rf.model, test_x)$predictions
#' calc_rel_predicted(train_y, predicted_train_y=train_rf_model$predicted)
#' calc_rel_predicted(train_y=train_y, predicted_train_y=train_rf_model$predicted,
#'                    test_y=test_y, predicted_test_y=predicted_test_y,
#'                    train_target_field="y",  test_target_field="test_y", outdir=NULL)
#' calc_rel_predicted(train_y=train_y, predicted_train_y=train_rf_model$predicted,
#'                    train_SampleIDs=as.character(1:60),
#'                    test_y=test_y, predicted_test_y=predicted_test_y,
#'                    test_SampleIDs=as.character(1:45),
#'                    train_target_field="y",  test_target_field="test_y", outdir=NULL)
#' @author Shi Huang
#' @export
calc_rel_predicted<-function(train_y, predicted_train_y, train_SampleIDs=NULL,
                             test_y=NULL, predicted_test_y=NULL, test_SampleIDs=NULL,
                             train_prefix="train", test_prefix="test",
                             train_target_field="y", test_target_field="y", outdir=NULL){
  spl_train <- smooth.spline(train_y, predicted_train_y)
  train_relTrain <- residuals(spl_train)
  train_fittedTrain <- fitted(spl_train)
  relTrain_data<-train_relTrain_data <- data.frame(y=train_y, predicted_y=predicted_train_y,
                                                   rel_predicted_y=train_relTrain,
                                                   fitted_predicted_y=train_fittedTrain)
  if(!is.null(train_SampleIDs)){
    if(nrow(relTrain_data)!=length(train_SampleIDs)){
      stop("Please make sure that sample IDs match with y or predicted y in the train data.")
    }else{
      relTrain_data<-train_relTrain_data<-data.frame(SampleIDs=train_SampleIDs, relTrain_data)
    }
  }
  if(!is.null(outdir)){
    .write_tsv(relTrain_data, paste(outdir, train_prefix,".Relative_",train_target_field,".results.xls",sep=""))
  }

  if(!is.null(test_y) & !is.null(predicted_test_y)){
    test_fitted<-predict(spl_train, test_y)$y
    test_relTrain <- predicted_test_y - test_fitted
    test_relTrain_data <- data.frame(y=test_y, predicted_y=predicted_test_y,
                                     rel_predicted_y=test_relTrain,
                                     fitted_predicted_y=test_fitted)
    if(!is.null(test_SampleIDs)){
      if(nrow(test_relTrain_data)!=length(test_SampleIDs)){
        stop("Please make sure that sample IDs match with y or predicted y in the test data.")
        }else{
          test_relTrain_data<-data.frame(SampleIDs=test_SampleIDs, test_relTrain_data)
        }
    }
    relTrain_data<-rbind(train_relTrain_data, test_relTrain_data)
    DataSet <- factor(c(rep(train_prefix,length(train_relTrain)), rep(test_prefix,length(test_relTrain)) ))
    relTrain_data<-data.frame(relTrain_data, DataSet)

    if(!is.null(outdir)){
      .write_tsv(relTrain_data, paste(outdir, train_prefix,"-",test_prefix,".Relative_",train_target_field,".results.xls",sep=""))
    }
  }
    return(relTrain_data)
}

#' @title plot_rel_predicted
#' @description Calculate the relative predicted values to the spline fit in the training data.
#' @param relTrain_data The output dataframe of \code{calc_rel_predicted}.
#' @param prefix The prefix of a train dataset, or the prefixes of both train and test datasets.
#' @param target_field A string indicating the target field in the metadata for regression.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' train_x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' train_y<- 1:60
#' test_x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209)))))
#' test_y<- 1:45
#' train_rf_model<-rf.out.of.bag(train_x, train_y)
#' predicted_test_y<-predict(train_rf_model$rf.model, test_x)$predictions
#' relTrain_data<-calc_rel_predicted(train_y, train_rf_model$predicted,
#'                                   test_y=test_y, predicted_test_y=predicted_test_y)
#' plot_rel_predicted(relTrain_data)
#' plot_rel_predicted(relTrain_data, prefix=c("train", "test"))
#' @author Shi Huang
#' @export
plot_rel_predicted <- function(relTrain_data, prefix="train", target_field="value", outdir=NULL){
  if(length(prefix)==1){
    if("DataSet" %in% colnames(relTrain_data)) relTrain_data<-relTrain_data[relTrain_data$DataSet==prefix, , drop=FALSE]
    p<-ggplot(relTrain_data, aes(x=.data$y, y=.data$rel_predicted_y))+
      ylab(paste("Relative prediceted ",target_field,sep=""))+
      xlab(paste("Observed ",target_field,sep=""))+
      geom_point(alpha=0.1)+
      geom_hline(yintercept=0)+
      theme_bw()
  }else{
    p<-ggplot(relTrain_data, aes(x=.data$y, y=.data$rel_predicted_y, color=.data$DataSet))+
      ylab(paste("Relative prediceted ",target_field,sep=""))+
      xlab(paste("Observed ",target_field,sep=""))+
      geom_point(alpha=0.1)+
      geom_hline(yintercept=0)+
      theme_bw()
    prefix=paste(prefix, collapse = "-")
  }
  if(!is.null(outdir)){
    ggsave(filename=paste(outdir, prefix, ".", target_field, ".obs_vs_relative_pred.scatterplot.pdf",sep=""), plot=p, height=4, width=4)
  }
  invisible(p)
}

#' @title boxplot_rel_predicted_train_vs_test
#' @description Make the boxplot the relative predicted values in both train and test datasets.
#' @param relTrain_data The output dataframe of \code{calc_rel_predicted}.
#' @param train_prefix The prefix for the dataset in the training data.
#' @param test_prefix The prefix for the dataset in the testing data.
#' @param train_target_field A string indicating the target field of the training metadata for regression.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' train_x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' train_y<- 1:60
#' test_x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209)))))
#' test_y<- 1:45
#' train_rf_model<-rf.out.of.bag(train_x, train_y)
#' predicted_test_y<-predict(train_rf_model$rf.model, test_x)$predictions
#' relTrain_data<-calc_rel_predicted(train_y, predicted_train_y=train_rf_model$predicted,
#'                                   train_SampleIDs=NULL,
#'                                   test_y, predicted_test_y,
#'                                   test_SampleIDs=NULL)
#' boxplot_rel_predicted_train_vs_test(relTrain_data)
#' @author Shi Huang
#' @export
boxplot_rel_predicted_train_vs_test<-function(relTrain_data, train_target_field="value",
                                              train_prefix="train", test_prefix="test", outdir=NULL){
    NoTestDataset<-!"DataSet" %in% colnames(relTrain_data)
    if(NoTestDataset) stop("Test dataset should be included for residuals comparison between train and test datasets!")
    p_w<-formatC(stats::wilcox.test(rel_predicted_y~DataSet, data = relTrain_data)$p.value,digits=4,format="g")
    p<-ggplot(relTrain_data, aes(x=.data$DataSet, y=.data$rel_predicted_y)) +
      geom_violin(aes(color=.data$DataSet))+
      geom_boxplot(outlier.shape = NA, width=0.4)+
      geom_jitter(position=position_jitter(width=0.2),alpha=0.1) +
      geom_hline(yintercept=0)+
      theme_bw()+
      ggtitle(paste("Wilcoxon Rank Sum Test:\n P=", p_w, sep=""))+
      xlab("Data sets") +
      ylab(paste("Relative predicted", train_target_field))+
      theme(legend.position="none")
    if(!is.null(outdir)){
      ggsave(filename= paste(outdir, train_prefix,"-",test_prefix, ".Relative_",train_target_field,".boxplot.pdf",sep=""),
             plot=p, width=3, height=4)
    }
    invisible(p)
}
