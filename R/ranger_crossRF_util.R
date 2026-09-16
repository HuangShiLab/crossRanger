#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach getDoParRegistered getDoParWorkers
#' @importFrom foreach %dopar%
#' @importFrom parallel detectCores
#' @importFrom caret confusionMatrix
#' @importFrom stats complete.cases

utils::globalVariables(c("i"))

# Register a parallel backend for foreach. A backend already registered by the user
# (e.g., doFuture, or a cluster of their own) is kept, so that users can choose how to parallelize.
# Otherwise doParallel is registered: it forks the R session on Unix (as doMC did) and starts a
# PSOCK cluster on Windows. By default all but four cores are used (at least one), and at most
# two cores when the package is checked (_R_CHECK_LIMIT_CORES_).
.register_cores <- function(n_cores=NULL){
  # a backend registered by the user is used as it is and must not be stopped by us
  if(is.null(n_cores) && foreach::getDoParRegistered() && foreach::getDoParWorkers() > 1)
    return(invisible(structure(foreach::getDoParWorkers(), own=FALSE)))
  if(is.null(n_cores)){
    n_detected <- parallel::detectCores()
    n_cores <- if(is.na(n_detected)) 1L else n_detected - 4L
    # while the package is checked, run sequentially: starting a cluster would dominate the run
    # time of the examples on Windows, and CRAN limits checks to two cores anyway
    if(nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) n_cores <- 1L
  }
  n_cores <- max(1L, as.integer(n_cores))
  if(n_cores == 1L){
    # a single core needs no cluster at all: on Windows even registerDoParallel(cores=1) would
    # start a PSOCK worker, whose start-up dominates the run time of small analyses
    foreach::registerDoSEQ()
    return(invisible(structure(1L, own=FALSE)))
  }
  doParallel::registerDoParallel(cores=n_cores)
  invisible(structure(n_cores, own=TRUE))
}

# On Windows, registerDoParallel(cores=) starts an implicit PSOCK cluster whose socket connections
# stay open; they have to be closed again, otherwise they leak out of the function (R CMD check
# reports "connections left open"). On Unix the backend forks and opens no connection, so nothing
# is done there and the behaviour is unchanged. registerDoSEQ() prevents the stopped cluster from
# being picked up as an already registered backend by the next call.
.stop_implicit_cluster <- function(){
  if(.Platform$OS.type == "windows"){
    try(doParallel::stopImplicitCluster(), silent=TRUE)
    foreach::registerDoSEQ()
  }
  invisible(NULL)
}

# ranger uses all cores by default. Within parallel workers, the cores are shared among the workers
# to avoid oversubscription; the results of ranger do not depend on the number of threads.
.worker_dots <- function(n_workers, ...){
  dots <- list(...)
  if(n_workers > 1 && is.null(dots[["num.threads"]])){
    n_detected <- parallel::detectCores()
    dots[["num.threads"]] <- if(is.na(n_detected)) 1L else max(1L, n_detected %/% n_workers)
  }
  dots
}

# Match the rows of a feature table and a metadata table by sample IDs (rownames) if both have them,
# otherwise assume that they are in the same order.
.align_samples <- function(df, metadata){
  has_ids <- function(obj) !is.null(rownames(obj)) && !(is.data.frame(obj) && .row_names_info(obj) < 0)
  if(has_ids(df) && has_ids(metadata)){
    shared <- intersect(rownames(df), rownames(metadata))
    if(length(shared)==0) stop("No sample IDs (rownames) are shared by the feature table and the metadata.")
    if(length(shared) < max(nrow(df), nrow(metadata))){
      message("Keeping ", length(shared), " samples shared by the feature table and the metadata.")
    }else if(!identical(rownames(df), rownames(metadata))){
      message("The metadata were reordered to match the sample IDs of the feature table.")
    }
    return(list(df_idx=match(shared, rownames(df)), md_idx=match(shared, rownames(metadata))))
  }
  if(nrow(df)!=nrow(metadata))
    stop("The feature table and the metadata have different numbers of samples and no sample IDs (rownames) to match them.")
  list(df_idx=seq_len(nrow(df)), md_idx=seq_len(nrow(metadata)))
}

# Combine the inputs into one feature table and one metadata table:
# `df` can be a feature table or a named list of feature tables (one per dataset),
# and `metadata` a data.frame or a list of data.frames in the same order as the list of feature tables.
.prepare_datasets <- function(df, metadata, s_category=NULL, required_cols=NULL, verbose=FALSE){
  dataset_label <- NULL
  if(is.list(df) && !is.data.frame(df)){
    if(is.null(names(df))) names(df) <- paste0("dataset", seq_along(df))
    df_list <- harmonize_features(df, verbose=verbose)
    dataset_label <- rep(names(df_list), vapply(df_list, nrow, integer(1)))
    if(is.list(metadata) && !is.data.frame(metadata)){
      if(length(metadata)!=length(df_list))
        stop("The list of metadata should have the same length as the list of feature tables.")
      shared_cols <- Reduce(intersect, lapply(metadata, colnames))
      metadata <- do.call(rbind, lapply(unname(metadata), function(m) data.frame(m, check.names=FALSE)[, shared_cols, drop=FALSE]))
    }
    df <- do.call(rbind, unname(df_list))
  }
  df <- .as_feature_df(df)
  metadata <- data.frame(metadata, check.names=FALSE)
  idx <- .align_samples(df, metadata)
  df <- df[idx$df_idx, , drop=FALSE]
  metadata <- metadata[idx$md_idx, , drop=FALSE]
  if(is.null(s_category)){
    if(is.null(dataset_label)) stop("s_category is required when df is a single feature table.")
    dataset_label <- dataset_label[idx$df_idx]
    s_category <- "dataset"
    metadata[, s_category] <- factor(dataset_label, levels=unique(dataset_label))
  }
  used_cols <- c(s_category, required_cols)
  missing_cols <- setdiff(used_cols, colnames(metadata))
  if(length(missing_cols) > 0) stop("Column(s) not found in the metadata: ", paste(missing_cols, collapse=", "))
  complete <- complete.cases(metadata[, used_cols, drop=FALSE])
  if(any(!complete)){
    message("Removing ", sum(!complete), " samples with missing values in: ", paste(used_cols, collapse=", "))
    df <- df[complete, , drop=FALSE]
    metadata <- metadata[complete, , drop=FALSE]
  }
  list(df=df, metadata=metadata, s_category=s_category)
}

# Fit a RF model with the requested validation strategy.
.fit_rf <- function(x, y, nfolds=1, cv_type="stratified", groups=NULL, ntree=500, verbose=FALSE, imp_pvalues=FALSE, ...){
  if(cv_type=="logo"){
    return(rf.cross.validation(x, y, nfolds=factor(groups), ntree=ntree, verbose=verbose, imp_pvalues=imp_pvalues, ...))
  }
  if(cv_type=="group"){
    if(nfolds < 2) nfolds <- 5
    return(rf.cross.validation(x, y, nfolds=nfolds, groups=groups, ntree=ntree, verbose=verbose, imp_pvalues=imp_pvalues, ...))
  }
  if(nfolds==1){
    rf.out.of.bag(x, y, ntree=ntree, verbose=verbose, imp_pvalues=imp_pvalues, ...)
  }else{
    rf.cross.validation(x, y, nfolds=nfolds, ntree=ntree, verbose=verbose, imp_pvalues=imp_pvalues, ...)
  }
}

.validation_type <- function(nfolds, cv_type){
  if(cv_type=="logo") return("LOGO_CV")
  if(cv_type=="group") return("group_kfold_CV")
  if(nfolds==1) "OOB" else "stratified_kfold_CV"
}

# Order the columns of new data as the features of a model; features absent from the new data are set to 0.
.align_features <- function(newx, feature_names){
  newx <- .as_feature_df(newx)
  missing_features <- setdiff(feature_names, colnames(newx))
  if(length(missing_features) > 0){
    warning(length(missing_features), " of ", length(feature_names),
            " features used by the model are absent from the new data and were set to 0.", call.=FALSE)
    newx[, missing_features] <- 0
  }
  newx[, feature_names, drop=FALSE]
}

# Class probabilities and predicted classes of new data. For rf.cross.validation,
# the probabilities are averaged over the models of all folds.
.predict_clf <- function(rf_obj, newx){
  if(inherits(rf_obj, "ranger")){
    forests <- list(rf_obj)
    class_levels <- rf_obj$forest$levels
  }else{
    forests <- if(inherits(rf_obj, "rf.cross.validation")) rf_obj$rf.model else list(rf_obj$rf.model)
    class_levels <- levels(factor(rf_obj$y))
  }
  newx <- .align_features(newx, forests[[1]]$forest$independent.variable.names)
  prob_list <- lapply(forests, function(model){
    probs <- get.predict.probability.from.forest(model, newx)
    if(is.null(model$forest$levels)) colnames(probs) <- class_levels[seq_len(ncol(probs))]
    out <- matrix(0, nrow=nrow(newx), ncol=length(class_levels), dimnames=list(rownames(newx), class_levels))
    shared <- intersect(colnames(probs), class_levels)
    out[, shared] <- probs[, shared]
    out
  })
  probs <- Reduce(`+`, prob_list)/length(prob_list)
  predicted <- factor(class_levels[max.col(probs, ties.method="first")], levels=class_levels)
  list(probabilities=probs, predicted=predicted)
}

# Predicted values of new data. For rf.cross.validation (or a list of ranger models),
# the predictions are averaged over the models of all folds.
.predict_reg <- function(model, newx){
  forests <- if(inherits(model, c("rf.out.of.bag", "rf.cross.validation"))) model$rf.model else model
  if(inherits(forests, "ranger")) forests <- list(forests)
  newx <- .align_features(newx, forests[[1]]$forest$independent.variable.names)
  preds <- vapply(forests, function(m) as.numeric(predict(m, newx)$predictions), numeric(nrow(newx)))
  rowMeans(matrix(preds, nrow=nrow(newx)))
}

.dataset_names <- function(x_list, model_list){
  datasets <- names(x_list)
  if(is.null(datasets)) datasets <- names(model_list)
  if(is.null(datasets)) datasets <- paste0("dataset", seq_along(x_list))
  datasets
}

#' @title rf_clf.pairwise
#' @description Perform pairwise rf classfication for a data matrix between all pairs of group levels
#' @param df A data matrix or data.frame
#' @param f A factor with more than two levels.
#' @param nfolds The number of folds in the cross validation.
#' @param ntree The number of trees.
#' @param verbose The boolean value indicating if the computation status and estimated runtime shown.
#' @return A summary table containing the performance of Random Forest classification models
#' @seealso ranger
#' @examples
#' df <- data.frame(t(rmultinom(48, 160, c(.001,.6,.2,.3,.299))) + 0.65)
#' f<-factor(rep(c("A", "B", "C", "D"), each=12))
#' rf_clf.pairwise(df, f, ntree=500)
#' @author Shi Huang
#' @export
rf_clf.pairwise <- function (df, f, nfolds=3, ntree=5000, verbose=FALSE) {
  rf_compare_levels <- function(df, f, i=1, j=2, nfolds=1, ntree=500, verbose=FALSE) {
    df_ij <- df[which(as.integer(f) == i | as.integer(f) == j), ]
    f_ij <- factor(f[which(as.integer(f) == i | as.integer(f) == j)])
    if(nfolds==1){
      oob<-rf.out.of.bag(df_ij, f_ij, verbose=verbose, ntree=ntree)
    }else{
      oob<-rf.cross.validation(df_ij, f_ij, nfolds=nfolds, verbose=verbose, ntree=ntree)
    }
    positive_class=levels(f)[i]
    cat("\nTraining dataset: ", levels(f)[i], "-", levels(f)[j] ,"\n\n")
    conf<-caret::confusionMatrix(data=oob$predicted, f_ij, positive=positive_class)
    acc<-conf$overall[1]
    kappa_oob<-conf$overall[2]
    cat("Accuracy in the cross-validation: ", acc ,"\n")
    if(nlevels(f_ij)==2){rf_AUROC<-get.auroc(oob$probabilities[, positive_class], f_ij, positive_class)}else{rf_AUROC<-NA}
    if(nlevels(f_ij)==2){rf_AUPRC<-get.auprc(oob$probabilities[, positive_class], f_ij, positive_class)}else{rf_AUPRC<-NA}
    cat("AUROC in the cross-validation: ", rf_AUROC ,"\n")
    cat("AUPRC in the cross-validation: ", rf_AUPRC ,"\n")
    c("AUROC"=rf_AUROC, "AUPRC"=rf_AUPRC, acc, kappa_oob, conf$byClass)
  }
  if(nlevels(f)==2){
    out_summ<-rf_compare_levels(df, f, i=1, j=2, nfolds, ntree, verbose)
  }else{
    level_names<-levels(f)
    ix <- stats::setNames(seq_along(level_names), level_names)
    out_list<-outer(ix[-1L], ix[-length(ix)],
                    function(ivec, jvec) sapply(seq_along(ivec),
                                                function(k) {
                                                  i <- ivec[k]
                                                  j <- jvec[k]
                                                  if (i > j) {
                                                    rf_compare_levels(df, f, i, j, nfolds, ntree, verbose)
                                                  }else{NA}
                                                }))
    dataset_name<-outer(dimnames(out_list)[[1]],dimnames(out_list)[[2]], paste, sep="__")
    out_rownames<-dataset_name[lower.tri(dataset_name, diag = T)]
    out_summ<-do.call(rbind, out_list[lower.tri(out_list, diag = T)]); rownames(out_summ)<-out_rownames
  }
  class(out_summ)<-"rf_clf.pairwise"
  out_summ
}

#' @title rf_clf.by_datasets
#' @description It runs standard random forests with oob estimation or cross-validation for classification of
#' c_category in each the sub-datasets splited by the s_category. The output includes the rf models,
#' a within-study performance summary (accuracy, AUROC, AUPRC, Kappa statistics, sensitivity, specificity, F1 etc.)
#' and the statistics of all features.
#' @param df Training data: a data.frame, or a named list of data.frames (e.g., one feature table per study).
#' A list of feature tables is harmonized to the features shared by all datasets (see \code{harmonize_features}).
#' @param metadata Sample metadata with at least two columns: a data.frame, or a list of data.frames in the same order as the list of feature tables.
#' If both the feature table and the metadata have sample IDs as rownames, samples are matched by IDs.
#' @param s_category A string indicates the category in the sample metadata: a ‘factor’ defines the sample grouping for data spliting
#' (e.g., study, cohort, body site or timepoint). If NULL and \code{df} is a list, the datasets are defined by the list names.
#' @param c_category A indicates the category in the sample metadata as a responsive vector: if a 'factor', rf classification is performed in each of splited datasets.
#' @param positive_class A string indicates one class in the 'c_category' column of metadata.
#' @param nfolds The number of folds in the cross validation. If 1, out-of-bag estimation is used.
#' @param clr_transform A boolean value indicating if the clr-transformation applied.
#' @param rf_imp_pvalues A boolean value indicating if compute both importance score and pvalue for each feature.
#' @param verbose A boolean value indicating if show computation status and estimated runtime.
#' @param ntree The number of trees.
#' @param p.adj.method The p-value correction method, default is "BH".
#' @param q_cutoff The cutoff of q values for features, the default value is 0.05.
#' @param cv_type The cross-validation strategy within each dataset: "stratified" (stratified k-fold CV, or OOB if nfolds=1),
#' "group" (group k-fold CV) or "logo" (leave-one-group-out CV). "group" and "logo" require \code{g_category}.
#' @param g_category A string indicating the grouping column in metadata for group k-fold or leave-one-group-out CV
#' (e.g., subject ID in a longitudinal study, or family ID).
#' @param n_cores The number of cores used for parallel computation. By default, all but four cores.
#' @param ... Other parameters applicable to \code{ranger}.
#' @return An object of class \code{rf_clf.by_datasets}, a list including \code{x_list}, \code{y_list},
#' \code{datasets}, \code{sample_size}, \code{rf_model_list}, \code{rf_AUROC}, \code{rf_AUPRC},
#' \code{feature_imps_list} and the performance summary \code{perf_summ}.
#' @seealso ranger rf_clf.cross_appl rf_clf.lodo
#' @examples
#' df <- data.frame(rbind(t(rmultinom(14, 14*5, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(16, 16*5, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(30, 30*5, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(30, 30*5, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(30, 30*5, c(.001,.6,.42,.58,.299)))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 30), rep("B", 30), rep("C", 30), rep("D", 30))),
#'                      f_c=factor(c(rep("C", 14), rep("H", 16), rep("C", 14), rep("H", 16),
#'                                   rep("C", 14), rep("H", 16), rep("C", 14), rep("H", 16))),
#'                      f_d=factor(rep(c(rep("a", 10), rep("b", 10), rep("c", 10)), 4)),
#'                      subject=factor(rep(1:60, each=2)))
#' res <- rf_clf.by_datasets(df, metadata, s_category='f_s', c_category='f_c',
#'                           positive_class="C", ntree=100)
#' res$perf_summ
#' \donttest{
#' rf_clf.by_datasets(df, metadata, s_category='f_s', c_category='f_c',
#'                    positive_class="C", rf_imp_pvalues=TRUE)
#' rf_clf.by_datasets(df, metadata, s_category='f_s', c_category='f_d')
#' # group k-fold CV: samples from the same subject are kept in the same fold
#' rf_clf.by_datasets(df, metadata, s_category='f_s', c_category='f_c', positive_class="C",
#'                    nfolds=3, cv_type="group", g_category="subject")$perf_summ
#' # a list of feature tables (one per study) and a list of metadata
#' df_list <- split(df, metadata$f_s)
#' md_list <- split(metadata, metadata$f_s)
#' rf_clf.by_datasets(df_list, md_list, c_category='f_c', positive_class="C")$perf_summ
#' }
#' @author Shi Huang
#' @export
rf_clf.by_datasets<-function(df, metadata, s_category=NULL, c_category, positive_class=NA,
                             rf_imp_pvalues=FALSE, clr_transform=TRUE, nfolds=1, verbose=FALSE, ntree=500,
                             p.adj.method = "BH", q_cutoff=0.05,
                             cv_type=c("stratified", "group", "logo"), g_category=NULL, n_cores=NULL, ...){
  cv_type <- match.arg(cv_type)
  if(cv_type!="stratified" && is.null(g_category))
    stop("g_category (e.g., a subject or family ID column in metadata) is required for group k-fold or leave-one-group-out CV.")
  prep <- .prepare_datasets(df, metadata, s_category, c(c_category, g_category), verbose)
  df <- prep$df; metadata <- prep$metadata; s_category <- prep$s_category
  s <- factor(metadata[, s_category])
  y_list<-split(factor(metadata[, c_category]), s)
  x_list<-split(df, s)
  g_list<-if(is.null(g_category)) NULL else split(metadata[, g_category], s)
  datasets<-levels(s)
  L<-length(y_list)
  positive_class<-ifelse(is.na(positive_class), levels(factor(metadata[, c_category]))[1], positive_class)
  n_classes <- vapply(y_list, function(y) nlevels(droplevels(y)), integer(1))
  if(any(n_classes < 2))
    stop("Less than two classes of '", c_category, "' in dataset(s): ", paste(datasets[n_classes < 2], collapse=", "))
  has_positive <- vapply(y_list, function(y) positive_class %in% levels(droplevels(y)), logical(1))
  if(!all(has_positive))
    stop("The positive class '", positive_class, "' is absent from dataset(s): ", paste(datasets[!has_positive], collapse=", "))
  # 1. sample size of all datasets
  sample_size<-as.numeric(table(s))
  n_workers<-.register_cores(n_cores)
  if(isTRUE(attr(n_workers, "own"))) on.exit(.stop_implicit_cluster(), add=TRUE)
  rf_dots<-.worker_dots(n_workers, ...)
  oper<-foreach(i=1:L) %dopar% {
    x<-x_list[[i]]
    y<-droplevels(y_list[[i]])
    # 2. RF model and its performance
    oob<-do.call(.fit_rf, c(list(x, y, nfolds=nfolds, cv_type=cv_type, groups=g_list[[i]], ntree=ntree,
                                 verbose=verbose, imp_pvalues=rf_imp_pvalues), rf_dots))
    perf<-.clf_perf_row(oob$predicted, oob$y, oob$probabilities, positive_class)
    # 3. # of significantly differential abundant features between health and disease
    out<-BetweenGroup.test(x, y, clr_transform=clr_transform, positive_class=positive_class, p.adj.method = p.adj.method, q_cutoff=q_cutoff)
    feature_imps<-data.frame(feature=rownames(out), dataset=rep(datasets[i], ncol(x)),
                             rf_imps=.feature_importances(oob), out)
    list(oob=oob, perf=perf, feature_imps=feature_imps)
  }
  perf_summ<-data.frame(Dataset=datasets, Sample_size=sample_size, Validation_type=.validation_type(nfolds, cv_type),
                        do.call(rbind, lapply(oper, `[[`, "perf")), row.names=NULL, check.names=FALSE)
  result<-list()
  result$x_list<-x_list
  result$y_list<-y_list
  result$datasets<-datasets
  result$sample_size<-sample_size
  result$rf_model_list<-setNames(lapply(oper, `[[`, "oob"), datasets)
  result$rf_AUROC<-setNames(perf_summ$AUROC, datasets)
  result$rf_AUPRC<-setNames(perf_summ$AUPRC, datasets)
  result$feature_imps_list<-setNames(lapply(oper, `[[`, "feature_imps"), datasets)
  result$perf_summ<-perf_summ
  result$positive_class<-positive_class
  class(result)<-"rf_clf.by_datasets"
  return(result)
}

#' @title rf_reg.by_datasets
#' @description It runs standard random forests with oob estimation or cross-validation for regression of
#' \code{c_category} in each the sub-datasets splited by the \code{s_category}.
#' The output includes the rf models, predicted values and a performance summary (MSE, RMSE, MAE, MAPE, R squared etc.).
#' @param df Training data: a data.frame, or a named list of data.frames (e.g., one feature table per study).
#' @param metadata Sample metadata with at least two columns: a data.frame, or a list of data.frames in the same order as the list of feature tables.
#' @param s_category A string indicates the category in the sample metadata: a ‘factor’ defines the sample grouping for data spliting.
#' If NULL and \code{df} is a list, the datasets are defined by the list names.
#' @param c_category A string indicates the numeric variable in the sample metadata for rf regression in each of splited datasets.
#' @param rf_imp_pvalues A boolean value indicate if compute both importance score and pvalue for each feature.
#' @param ntree The number of trees.
#' @param nfolds The number of folds in the cross validation. If 1, out-of-bag estimation is used.
#' @param verbose Show computation status and estimated runtime.
#' @param cv_type The cross-validation strategy within each dataset: "stratified", "group" or "logo".
#' @param g_category A string indicating the grouping column in metadata for group k-fold or leave-one-group-out CV.
#' @param n_cores The number of cores used for parallel computation. By default, all but four cores.
#' @param ... Other parameters applicable to \code{ranger}.
#' @return An object of class \code{rf_reg.by_datasets}.
#' The performance metrics of cross-validation are computed from the pooled predictions of all folds.
#' @seealso ranger rf_reg.cross_appl rf_reg.lodo
#' @examples
#' df <- data.frame(rbind(t(rmultinom(14, 14*5, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(16, 16*5, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(30, 30*5, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(30, 30*5, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(30, 30*5, c(.001,.6,.42,.58,.299)))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 60), rep("B", 60))),
#'                      f_s1=factor(c(rep(TRUE, 60), rep(FALSE, 60))),
#'                      f_c=factor(c(rep("C", 30), rep("H", 30), rep("D", 30), rep("P", 30))),
#'                      age=c(1:60, 2:61)
#'                      )
#'
#' reg_res<-rf_reg.by_datasets(df, metadata, nfolds=5, s_category='f_s', c_category='age')
#' reg_res$perf_summ
#' @author Shi Huang
#' @export
rf_reg.by_datasets<-function(df, metadata, s_category=NULL, c_category, nfolds=5,
                             rf_imp_pvalues=FALSE, verbose=FALSE, ntree=500,
                             cv_type=c("stratified", "group", "logo"), g_category=NULL, n_cores=NULL, ...){
  cv_type <- match.arg(cv_type)
  if(cv_type!="stratified" && is.null(g_category))
    stop("g_category (e.g., a subject or family ID column in metadata) is required for group k-fold or leave-one-group-out CV.")
  prep <- .prepare_datasets(df, metadata, s_category, c(c_category, g_category), verbose)
  df <- prep$df; metadata <- prep$metadata; s_category <- prep$s_category
  y<-metadata[, c_category]
  if(!is.numeric(y)) y<-suppressWarnings(as.numeric(as.character(y)))
  if(any(is.na(y))) stop("The target variable '", c_category, "' should be numeric for regression.")
  s <- factor(metadata[, s_category])
  y_list<-split(y, s)
  x_list<-split(df, s)
  g_list<-if(is.null(g_category)) NULL else split(metadata[, g_category], s)
  datasets<-levels(s)
  L<-length(y_list)
  # sample size of all datasets
  sample_size<-as.numeric(table(s))
  n_workers<-.register_cores(n_cores)
  if(isTRUE(attr(n_workers, "own"))) on.exit(.stop_implicit_cluster(), add=TRUE)
  rf_dots<-.worker_dots(n_workers, ...)
  oper<-foreach(i=1:L) %dopar% {
    do.call(.fit_rf, c(list(x_list[[i]], y_list[[i]], nfolds=nfolds, cv_type=cv_type, groups=g_list[[i]],
                            ntree=ntree, verbose=verbose, imp_pvalues=rf_imp_pvalues), rf_dots))
  }
  reg_metrics <- c("MSE", "RMSE", "nRMSE", "MAE", "MAPE", "MASE", "Spearman_rho", "R_squared", "Adj_R_squared")
  perf_mat <- do.call(rbind, lapply(seq_len(L), function(i)
    unlist(get.reg.performance(oper[[i]]$predicted, y_list[[i]], n_features=ncol(df))[reg_metrics])))
  result<-list()
  result$x_list<-x_list
  result$y_list<-y_list
  result$datasets<-datasets
  result$sample_size<-sample_size
  result$rf_model_list<-setNames(lapply(oper, `[[`, "rf.model"), datasets)
  result$rf_predicted<-setNames(lapply(oper, `[[`, "predicted"), datasets)
  result$feature_imps_list<-setNames(lapply(oper, .feature_importances), datasets)
  for(metric in reg_metrics) result[[paste0("rf_", metric)]] <- setNames(perf_mat[, metric], datasets)
  result$perf_summ<-data.frame(Dataset=datasets, Sample_size=sample_size, Validation_type=.validation_type(nfolds, cv_type),
                               perf_mat, row.names=NULL, check.names=FALSE)
  class(result)<-"rf_reg.by_datasets"
  return(result)
}

#' @title rf_clf.cross_appl
#' @description Based on pre-computed rf models classifying 'c_category' in each the sub-datasets splited by the 's_category',
#' perform cross-datasets application of the rf models. The inputs are precalculated
#' rf models, and the outputs include accuracy, auc and Kappa statistics.
#' @param rf_model_list A list of rf.model objects from \code{rf.out.of.bag} or \code{rf.cross.validation}.
#' For \code{rf.cross.validation}, the predicted probabilities are averaged over the models of all folds.
#' @param x_list A list of training datasets usually in the format of data.frame.
#' Features are matched to the model by feature IDs; features absent from a test dataset are set to 0 with a warning.
#' @param y_list A list of responsive vector for regression in the training datasets.
#' @param positive_class A string indicates one common class in each of elements in the y_list.
#' @return A object of class rf_clf.cross_appl including a performance summary (\code{perf_summ}),
#' the predicted values of all predictions (\code{predicted}), and a named list of
#' data.frames with the observed classes, predicted classes and probabilities (\code{prediction_list}).
#' @seealso ranger
#' @examples
#' df <- data.frame(rbind(t(rmultinom(14, 14*5, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(16, 16*5, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(30, 30*5, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(30, 30*5, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(30, 30*5, c(.001,.6,.42,.58,.299)))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 30), rep("B", 30), rep("C", 30), rep("D", 30))),
#'                      f_c=factor(c(rep("C", 14), rep("H", 16), rep("C", 14), rep("H", 16),
#'                                   rep("C", 14), rep("H", 16), rep("C", 14), rep("H", 16))),
#'                      f_d=factor(rep(c(rep("a", 10), rep("b", 10), rep("c", 10)), 4)))
#' res_list<-rf_clf.by_datasets(df, metadata, s_category='f_s', nfolds=5,
#'                              c_category='f_c', positive_class="C")
#' rf_model_list<-res_list$rf_model_list
#'
#' cross_rf<-rf_clf.cross_appl(rf_model_list, res_list$x_list, res_list$y_list, positive_class="C")
#' cross_rf$perf_summ
#' #--------------------
#' comp_group="A"
#' comps_res<-rf_clf.comps(df, f=metadata[, 'f_s'], comp_group, verbose=FALSE,
#'                         ntree=500, p.adj.method = "bonferroni", q_cutoff=0.05)
#' rf_clf.cross_appl(comps_res$rf_model_list,
#'                   x_list=comps_res$x_list,
#'                   y_list=comps_res$y_list,
#'                   positive_class=comp_group)
#' @author Shi Huang
#' @export
rf_clf.cross_appl<-function(rf_model_list, x_list, y_list, positive_class=NA){
  L<-length(rf_model_list)
  if(length(x_list)!=L || length(y_list)!=L) stop("The length of x list, y list and rf model list should be identical.")
  datasets <- .dataset_names(x_list, rf_model_list)
  positive_class<-ifelse(is.na(positive_class), levels(factor(y_list[[1]]))[1], positive_class)
  perf_summ<-data.frame(matrix(NA, ncol=18, nrow=L*L))
  colnames(perf_summ)<-c("Train_data", "Test_data", "Validation_type", .clf_metric_names)
  predicted<-matrix(list(), ncol=2, nrow=L*L)
  colnames(predicted)<-c("predicted", "probabilities")
  prediction_list<-list()
  for(i in 1:L){
    oob<-rf_model_list[[i]]
    if(!inherits(oob, c("rf.out.of.bag", "rf.cross.validation")))
      stop("Each element of rf_model_list should be an object of class rf.out.of.bag or rf.cross.validation.")
    train_y<-factor(oob$y)
    #---
    #  RF Training performance
    #---
    cat("\nTraining dataset: ", datasets[i] ,"\n\n")
    self_perf<-.clf_perf_row(oob$predicted, train_y, oob$probabilities, positive_class)
    cat("Accuracy in the self-validation: ", self_perf[["Accuracy"]] ,"\n")
    cat("AUROC in the self-validation: ", self_perf[["AUROC"]] ,"\n")
    cat("AUPRC in the self-validation: ", self_perf[["AUPRC"]] ,"\n")
    a=1+(i-1)*L
    perf_summ[a, 1:3]<-c(datasets[i], datasets[i], "self_validation")
    perf_summ[a, .clf_metric_names]<-self_perf
    predicted[a, 1][[1]]<-data.frame(test_y=train_y, pred_y=oob$predicted)
    predicted[a, 2][[1]]<-oob$probabilities
    prediction_list[[paste(datasets[i], datasets[i], sep="__VS__")]]<-
      data.frame(test_y=train_y, pred_y=oob$predicted, oob$probabilities, check.names=FALSE)
    loop_num<-1
    for(j in setdiff(1:L, i)){
      if(nrow(x_list[[j]])>0){
        newy<-droplevels(factor(y_list[[j]]))
        pred<-.predict_clf(oob, x_list[[j]])
        cat("Test dataset: ", datasets[j] ,"\n")
        if(all(levels(newy) %in% levels(train_y))){
          newy_aligned<-factor(as.character(newy), levels=levels(train_y))
          test_perf<-.clf_perf_row(pred$predicted, newy_aligned, pred$probabilities, positive_class)
          cat("Accuracy in the cross-applications: ", test_perf[["Accuracy"]] ,"\n")
          cat("AUROC in the cross-applications: ", test_perf[["AUROC"]] ,"\n")
          cat("AUPRC in the cross-applications: ", test_perf[["AUPRC"]] ,"\n")
        }else{
          test_perf<-rep(NA, length(.clf_metric_names))
        }
        perf_summ[a+loop_num, 1:3]<-c(datasets[i], datasets[j], "cross_application")
        perf_summ[a+loop_num, .clf_metric_names]<-test_perf
        predicted[a+loop_num, 1][[1]]<-data.frame(test_y=newy, pred_y=pred$predicted)
        predicted[a+loop_num, 2][[1]]<-pred$probabilities
        prediction_list[[paste(datasets[i], datasets[j], sep="__VS__")]]<-
          data.frame(test_y=newy, pred_y=pred$predicted, pred$probabilities, check.names=FALSE)
        loop_num<-loop_num+1
      }
    }
  }
  res<-list()
  res$perf_summ<-perf_summ
  res$predicted<-predicted
  res$prediction_list<-prediction_list
  class(res)<-"rf_clf.cross_appl"
  res
}

#' @title rf_reg.cross_appl
#' @description Based on pre-computed rf models regressing \code{c_category} in each the sub-datasets splited by the \code{s_category},
#' perform cross-datasets application of the rf models. The inputs are precalculated
#' rf regression models, x_list and y_list.
#' @param rf_list The output of \code{rf_reg.by_datasets}. For models from cross-validation,
#' the predictions are averaged over the models of all folds.
#' @param x_list A list of training datasets usually in the format of data.frame.
#' Features are matched to the model by feature IDs; features absent from a test dataset are set to 0 with a warning.
#' @param y_list A list of responsive vector for regression in the training datasets.
#' @return An object of class \code{rf_reg.cross_appl} including a performance summary (\code{perf_summ})
#' and the predicted values of all predictions (\code{predicted}).
#'
#' @seealso ranger
#' @examples
#' df <- data.frame(rbind(t(rmultinom(14, 14*5, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(16, 16*5, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(30, 30*5, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(30, 30*5, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(30, 30*5, c(.001,.6,.42,.58,.299)))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 60), rep("B", 60))),
#'                      f_s1=factor(c(rep(TRUE, 60), rep(FALSE, 60))),
#'                      f_c=factor(c(rep("C", 30), rep("H", 30), rep("D", 30), rep("P", 30))),
#'                      age=c(1:60, 2:61)
#'                      )
#'
#' table(metadata[, 'f_c'])
#' reg_res<-rf_reg.by_datasets(df, metadata, s_category='f_c', c_category='age', nfolds=1)
#' rf_reg.cross_appl(reg_res, x_list=reg_res$x_list, y_list=reg_res$y_list)
#' reg_res<-rf_reg.by_datasets(df, metadata, nfolds=5, s_category='f_c', c_category='age')
#' rf_reg.cross_appl(reg_res, x_list=reg_res$x_list, y_list=reg_res$y_list)
#' @author Shi Huang
#' @export
rf_reg.cross_appl<-function(rf_list, x_list, y_list){
  L<-length(rf_list$rf_model_list)
  if(length(x_list)!=L || length(y_list)!=L) stop("The length of x list, y list and rf model list should be identical.")
  if(!all(vapply(y_list, is.numeric, logical(1)))) stop("All elements in the y list should be numeric for regression.")
  datasets <- .dataset_names(x_list, rf_list$rf_model_list)
  perf_summ<-data.frame(matrix(NA, ncol=14, nrow=L*L))
  colnames(perf_summ)<-c("Train_data", "Test_data", "Validation_type", "Sample_size", "Min_acutal_value", "Max_acutal_value", "Min_predicted_value", "Max_predicted_value",
                         "MSE", "RMSE", "MAE", "MAPE", "Spearman_rho", "R_squared")
  perf_values <- function(y, pred_y){
    perf <- get.reg.performance(pred_y, y)
    c(length(y), range(y), range(pred_y), perf$MSE, perf$RMSE, perf$MAE, perf$MAPE, perf$Spearman_rho, perf$R_squared)
  }
  report <- function(values, type){
    cat("MSE in the ", type, ": ", values[5], "\n", sep="")
    cat("RMSE in the ", type, ": ", values[6], "\n", sep="")
    cat("MAE in the ", type, ": ", values[7], "\n", sep="")
    cat("MAE percentage in the ", type, ": ", values[8], "\n", sep="")
    cat("Spearman_rho in the ", type, ": ", values[9], "\n", sep="")
    cat("R squared in the ", type, ": ", values[10], "\n", sep="")
  }
  predicted<-list()
  for(i in 1:L){
    y<-y_list[[i]]
    rf_model<-rf_list$rf_model_list[[i]]
    train_pred<-rf_list$rf_predicted[[i]]
    #---  RF Training performance: MSE, MAE and R_squared
    cat("\nTraining dataset: ", datasets[i] ,"\n\n")
    a=1+(i-1)*L
    train_values<-perf_values(y, train_pred)
    report(train_values[-1] , "self-validation")
    perf_summ[a, 1:3]<-c(datasets[i], datasets[i], "self_validation")
    perf_summ[a, 4:14]<-train_values
    predicted[[paste(datasets[i], datasets[i], sep="__VS__")]]<-data.frame(test_y=y, pred_y=train_pred)
    loop_num<-1
    for(j in setdiff(1:L, i)){
      if(nrow(x_list[[j]])>0){
        newy<-y_list[[j]]
        pred_newy<-.predict_reg(rf_model, x_list[[j]])
        #---  RF test performance: MSE, MAE and R_squared
        cat("Test dataset: ", datasets[j] ,"\n")
        test_values<-perf_values(newy, pred_newy)
        report(test_values[-1], "cross-applications")
        perf_summ[a+loop_num, 1:3]<-c(datasets[i], datasets[j], "cross_application")
        perf_summ[a+loop_num, 4:14]<-test_values
        predicted[[paste(datasets[i], datasets[j], sep="__VS__")]]<-data.frame(test_y=newy, pred_y=pred_newy)
        loop_num<-loop_num+1
      }
    }
  }
  res<-list()
  res$perf_summ<-perf_summ
  res$predicted<-predicted
  class(res)<-"rf_reg.cross_appl"
  res
}

#' @title rf_clf.lodo
#' @description Leave-one-dataset-out (LODO) validation for random forest classification.
#' For each dataset defined by \code{s_category}, a model is trained on the samples pooled from all other datasets
#' and evaluated on the held-out dataset. This is the gold standard for meta-analytic performance estimation,
#' and, compared with \code{rf_clf.cross_appl} (single-dataset training), shows whether pooling multiple datasets
#' during training improves the generalizability.
#' @param df A data.frame, or a named list of data.frames (one feature table per dataset).
#' @param metadata A data.frame, or a list of data.frames in the same order as the list of feature tables.
#' @param s_category A string indicating the dataset (e.g., study or cohort) column in the metadata.
#' If NULL and \code{df} is a list, the datasets are defined by the list names.
#' @param c_category A string indicating the class column in the metadata.
#' @param positive_class A string indicates one class in the 'c_category' column of metadata.
#' @param ntree The number of trees.
#' @param verbose A boolean value indicating if show computation status.
#' @param n_cores The number of cores used for parallel computation. By default, all but four cores.
#' @param ... Other parameters applicable to \code{ranger}.
#' @return An object of class \code{rf_clf.lodo} including the performance on each held-out dataset (\code{perf_summ}),
#' the mean performance over held-out datasets (\code{perf_mean}), the predictions and the models.
#' @seealso rf_clf.by_datasets rf_clf.cross_appl plot_cross_appl
#' @examples
#' df <- data.frame(rbind(t(rmultinom(40, 200, c(.21,.6,.12,.38,.099))),
#'                        t(rmultinom(40, 200, c(.001,.6,.42,.58,.299))),
#'                        t(rmultinom(40, 200, c(.011,.6,.22,.28,.289)))))
#' metadata <- data.frame(study=factor(rep(c("S1", "S2", "S3"), each=40)),
#'                        disease=factor(rep(rep(c("Case", "Control"), each=20), 3)))
#' lodo <- rf_clf.lodo(df, metadata, s_category="study", c_category="disease", positive_class="Case")
#' lodo$perf_summ
#' @export
rf_clf.lodo <- function(df, metadata, s_category=NULL, c_category, positive_class=NA, ntree=500,
                        verbose=FALSE, n_cores=NULL, ...){
  prep <- .prepare_datasets(df, metadata, s_category, c_category, verbose)
  df <- prep$df; metadata <- prep$metadata; s_category <- prep$s_category
  s <- factor(metadata[, s_category])
  y <- factor(metadata[, c_category])
  datasets <- levels(s)
  if(length(datasets) < 2) stop("At least two datasets are required for leave-one-dataset-out validation.")
  positive_class <- ifelse(is.na(positive_class), levels(y)[1], positive_class)
  n_workers<-.register_cores(n_cores)
  if(isTRUE(attr(n_workers, "own"))) on.exit(.stop_implicit_cluster(), add=TRUE)
  rf_dots<-.worker_dots(n_workers, ...)
  oper <- foreach(i=seq_along(datasets)) %dopar% {
    test_idx <- which(s==datasets[i])
    train_y <- droplevels(y[-test_idx])
    model <- do.call(rf.out.of.bag, c(list(df[-test_idx, , drop=FALSE], train_y, ntree=ntree, verbose=verbose), rf_dots))
    pred <- .predict_clf(model, df[test_idx, , drop=FALSE])
    test_y <- droplevels(y[test_idx])
    if(all(levels(test_y) %in% levels(train_y))){
      perf <- .clf_perf_row(pred$predicted, factor(as.character(test_y), levels=levels(train_y)), pred$probabilities, positive_class)
    }else{
      perf <- setNames(rep(NA_real_, length(.clf_metric_names)), .clf_metric_names)
    }
    list(model=model, perf=perf, n_train=length(train_y), n_test=length(test_idx),
         prediction=data.frame(SampleID=rownames(df)[test_idx], test_y=test_y, pred_y=pred$predicted,
                               pred$probabilities, check.names=FALSE, row.names=NULL))
  }
  perf_mat <- do.call(rbind, lapply(oper, `[[`, "perf"))
  res <- list()
  res$perf_summ <- data.frame(Train_data=paste0("All_except_", datasets), Test_data=datasets, Validation_type="LODO",
                              Train_size=vapply(oper, `[[`, numeric(1), "n_train"),
                              Test_size=vapply(oper, `[[`, numeric(1), "n_test"),
                              perf_mat, row.names=NULL, check.names=FALSE)
  res$perf_mean <- colMeans(perf_mat, na.rm=TRUE)
  res$predicted <- setNames(lapply(oper, `[[`, "prediction"), datasets)
  res$rf_model_list <- setNames(lapply(oper, `[[`, "model"), paste0("All_except_", datasets))
  res$datasets <- datasets
  res$positive_class <- positive_class
  class(res) <- "rf_clf.lodo"
  res
}

#' @title rf_reg.lodo
#' @description Leave-one-dataset-out (LODO) validation for random forest regression.
#' For each dataset defined by \code{s_category}, a model is trained on the samples pooled from all other datasets
#' and evaluated on the held-out dataset.
#' @param df A data.frame, or a named list of data.frames (one feature table per dataset).
#' @param metadata A data.frame, or a list of data.frames in the same order as the list of feature tables.
#' @param s_category A string indicating the dataset (e.g., study or cohort) column in the metadata.
#' If NULL and \code{df} is a list, the datasets are defined by the list names.
#' @param c_category A string indicating the numeric target column in the metadata.
#' @param ntree The number of trees.
#' @param verbose A boolean value indicating if show computation status.
#' @param n_cores The number of cores used for parallel computation. By default, all but four cores.
#' @param ... Other parameters applicable to \code{ranger}.
#' @return An object of class \code{rf_reg.lodo} including the performance on each held-out dataset (\code{perf_summ}),
#' the mean performance over held-out datasets (\code{perf_mean}), the predictions and the models.
#' @seealso rf_reg.by_datasets rf_reg.cross_appl plot_cross_appl
#' @examples
#' df <- data.frame(rbind(t(rmultinom(40, 200, c(.21,.6,.12,.38,.099))),
#'                        t(rmultinom(40, 200, c(.001,.6,.42,.58,.299))),
#'                        t(rmultinom(40, 200, c(.011,.6,.22,.28,.289)))))
#' metadata <- data.frame(study=factor(rep(c("S1", "S2", "S3"), each=40)), age=rep(1:40, 3))
#' rf_reg.lodo(df, metadata, s_category="study", c_category="age")$perf_summ
#' @export
rf_reg.lodo <- function(df, metadata, s_category=NULL, c_category, ntree=500, verbose=FALSE, n_cores=NULL, ...){
  prep <- .prepare_datasets(df, metadata, s_category, c_category, verbose)
  df <- prep$df; metadata <- prep$metadata; s_category <- prep$s_category
  s <- factor(metadata[, s_category])
  y <- metadata[, c_category]
  if(!is.numeric(y)) y <- suppressWarnings(as.numeric(as.character(y)))
  if(any(is.na(y))) stop("The target variable '", c_category, "' should be numeric for regression.")
  datasets <- levels(s)
  if(length(datasets) < 2) stop("At least two datasets are required for leave-one-dataset-out validation.")
  reg_metrics <- c("MSE", "RMSE", "nRMSE", "MAE", "MAPE", "MASE", "Spearman_rho", "R_squared")
  n_workers<-.register_cores(n_cores)
  if(isTRUE(attr(n_workers, "own"))) on.exit(.stop_implicit_cluster(), add=TRUE)
  rf_dots<-.worker_dots(n_workers, ...)
  oper <- foreach(i=seq_along(datasets)) %dopar% {
    test_idx <- which(s==datasets[i])
    model <- .quietly(do.call(rf.out.of.bag, c(list(df[-test_idx, , drop=FALSE], y[-test_idx], ntree=ntree, verbose=verbose), rf_dots)))
    pred_y <- .predict_reg(model, df[test_idx, , drop=FALSE])
    list(model=model, n_train=length(y)-length(test_idx), n_test=length(test_idx),
         perf=unlist(get.reg.performance(pred_y, y[test_idx])[reg_metrics]),
         prediction=data.frame(SampleID=rownames(df)[test_idx], test_y=y[test_idx], pred_y=pred_y, row.names=NULL))
  }
  perf_mat <- do.call(rbind, lapply(oper, `[[`, "perf"))
  res <- list()
  res$perf_summ <- data.frame(Train_data=paste0("All_except_", datasets), Test_data=datasets, Validation_type="LODO",
                              Train_size=vapply(oper, `[[`, numeric(1), "n_train"),
                              Test_size=vapply(oper, `[[`, numeric(1), "n_test"),
                              perf_mat, row.names=NULL, check.names=FALSE)
  res$perf_mean <- colMeans(perf_mat, na.rm=TRUE)
  res$predicted <- setNames(lapply(oper, `[[`, "prediction"), datasets)
  res$rf_model_list <- setNames(lapply(oper, `[[`, "model"), paste0("All_except_", datasets))
  res$datasets <- datasets
  class(res) <- "rf_reg.lodo"
  res
}

#' @title generate.comps_datalist
#' @description To generate sub-datasets and sub-metadata by pairing
#' one level and each of other levels of the specified category in the metadata.
#' @param df Training data: a data.frame.
#' @param f A factor in the metadata with at least two levels (groups).
#' @param comp_group A string indicates the group in the f
#' @return A object of list with elements: a list of data.frame and a list of factor
#' @author Shi Huang
#'
generate.comps_datalist<-function(df, f, comp_group){
  all_other_groups<-levels(f)[which(levels(f)!=comp_group)]
  L<-length(all_other_groups)
  f_list<-list()
  df_list<-list()
  for(i in 1:L){
    f_list[[i]]<-factor(f[which(f==comp_group | f==all_other_groups[i])])
    df_list[[i]]<-df[which(f==comp_group | f==all_other_groups[i]), ]
  }
  names(f_list)<-names(df_list)<-paste(comp_group, all_other_groups, sep="_VS_")
  res<-list()
  res$f_list<-f_list
  res$df_list<-df_list
  res
}

#' @title rf_clf.comps
#' @description It runs standard random forests with oob estimation for classification of
#' one level VS all other levels of one category in the datasets.
#' The output includes a list of rf models for the sub datasets
#' and all important statistics for each of features.
#'
#' @param df Training data: a data.frame.
#' @param f A factor in the metadata with at least two levels (groups).
#' @param comp_group A string indicates the group in the f
#' @param clr_transform A boolean value indicating if the clr-transformation applied.
#' @param rf_imp_values A boolean value indicating if compute both importance score and pvalue for each feature.
#' @param verbose A boolean value indicating if show computation status and estimated runtime.
#' @param ntree The number of trees.
#' @param p.adj.method The p-value correction method, default is "bonferroni".
#' @param q_cutoff The cutoff of q values for features, the default value is 0.05.
#' @param nfolds The number of folds in the cross validation. If 1, out-of-bag estimation is used.
#' @param n_cores The number of cores used for parallel computation. By default, all but four cores.
#' @param ... Other parameters applicable to \code{ranger}.
#' @return An object of class \code{rf_clf.comps}.
#' @seealso ranger
#' @examples
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' f=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' comp_group="A"
#' comps_res<-rf_clf.comps(df, f, comp_group, verbose=FALSE, ntree=500,
#'                         p.adj.method = "bonferroni", q_cutoff=0.05)
#' comps_res$perf_summ
#' @author Shi Huang
#' @export
rf_clf.comps<-function(df, f, comp_group, verbose=FALSE, clr_transform=TRUE,
                       rf_imp_values=FALSE,ntree=500, p.adj.method = "bonferroni",
                       q_cutoff=0.05, nfolds=1, n_cores=NULL, ...){
  f<-factor(f)
  df<-.as_feature_df(df)
  if(!comp_group %in% levels(f)) stop("comp_group '", comp_group, "' is not a level of f.")
  all_other_groups<-levels(f)[which(levels(f)!=comp_group)]
  L<-length(all_other_groups)
  n_workers<-.register_cores(n_cores)
  if(isTRUE(attr(n_workers, "own"))) on.exit(.stop_implicit_cluster(), add=TRUE)
  rf_dots<-.worker_dots(n_workers, ...)
  oper<-foreach(i=1:L) %dopar% {
    idx<-which(f==comp_group | f==all_other_groups[i])
    sub_f<-factor(f[idx])
    sub_df<-df[idx, , drop=FALSE]
    dataset<-paste(comp_group, all_other_groups[i], sep="_VS_")
    oob<-do.call(.fit_rf, c(list(sub_df, sub_f, nfolds=nfolds, ntree=ntree, verbose=verbose, imp_pvalues=rf_imp_values), rf_dots))
    perf<-.clf_perf_row(oob$predicted, oob$y, oob$probabilities, comp_group)
    # # of significantly differential abundant features between two groups
    out<-BetweenGroup.test(sub_df, sub_f, clr_transform=clr_transform, q_cutoff=q_cutoff,
                           positive_class=comp_group, p.adj.method = p.adj.method)
    stats_out<-data.frame(feature=rownames(out),
                          dataset=rep(dataset, ncol(sub_df)), rf_imps=.feature_importances(oob),
                          out)
    list(x=sub_df, y=sub_f, sample_size=length(sub_f), dataset=dataset,
         oob=oob, perf=perf, stats_out=stats_out)
  }
  datasets<-vapply(oper, `[[`, character(1), "dataset")
  sample_size<-vapply(oper, `[[`, numeric(1), "sample_size")
  perf_mat<-do.call(rbind, lapply(oper, `[[`, "perf"))
  result<-list()
  result$x_list<-setNames(lapply(oper, `[[`, "x"), datasets)
  result$y_list<-setNames(lapply(oper, `[[`, "y"), datasets)
  result$sample_size<-sample_size
  result$datasets<-datasets
  result$rf_model_list<-setNames(lapply(oper, `[[`, "oob"), datasets)
  result$rf_AUROC<-setNames(perf_mat[, "AUROC"], datasets)
  result$rf_AUPRC<-setNames(perf_mat[, "AUPRC"], datasets)
  result$feature_imps_list<-setNames(lapply(oper, `[[`, "stats_out"), datasets)
  result$perf_summ<-data.frame(Dataset=datasets, Sample_size=sample_size, Validation_type=.validation_type(nfolds, "stratified"),
                               perf_mat, row.names=NULL, check.names=FALSE)
  result$positive_class<-comp_group
  class(result)<-"rf_clf.comps"
  return(result)
}
