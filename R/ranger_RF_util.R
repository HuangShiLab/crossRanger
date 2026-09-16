#' @importFrom stats model.matrix predict quantile sd cor cor.test median setNames
#' @importFrom Matrix Matrix
#' @importFrom ranger ranger importance_pvalues
#' @importFrom PRROC pr.curve roc.curve

.clf_metric_names <- c("Accuracy", "AUROC", "AUPRC", "Kappa", "Sensitivity", "Specificity",
                       "Pos_Pred_Value", "Neg_Pred_Value", "Precision", "Recall", "F1", "Prevalence",
                       "Detection_Rate", "Detection_Prevalence", "Balanced_Accuracy")

# Keep feature names as they are (e.g., "k__Bacteria;g__Prevotella"), so that models trained on
# one dataset can be applied to another dataset with the same feature IDs.
.as_feature_df <- function(x){
  if(is.null(colnames(x))) colnames(x) <- paste0("X", seq_len(ncol(x)))
  data.frame(x, check.names = FALSE)
}

.quietly <- function(expr){
  utils::capture.output(value <- expr)
  value
}

# Altmann's permutation p values of importance scores. ranger's formula interface requires syntactic
# column names, so the features are renamed for the permuted models; the scores are matched by position.
.altmann_pvalues <- function(model, data){
  names(data) <- make.names(names(data), unique=TRUE)
  importance_pvalues(model, method = "altmann", formula = y ~ ., data = data)
}

#' @title rf.out.of.bag
#' @description It runs standard random forests with out-of-bag error estimation for both classification and regression using \code{ranger}.
#' This is merely a wrapper that extracts relevant info from \code{ranger} output.
#' @param x Training data: data.matrix or data.frame.
#' @param y A response vector. If a factor, classification is assumed, otherwise regression is assumed.
#' @param ntree The number of trees.
#' @param sparse A boolean value indicates if the input matrix transformed into sparse matrix for rf modeling.
#' @param verbose A boolean value indicates if showing computation status and estimated runtime.
#' @param imp_pvalues If compute both importance score and pvalue for each feature.
#' @param ... Other parameters applicable to \code{ranger} (e.g., \code{mtry}, \code{num.threads}).
#' @return Object of class \code{rf.out.of.bag} with elements including the ranger object and critical metrics for model evaluation.
#'
#' @details
#' Ranger is a fast implementation of random forests (Breiman 2001) or recursive partitioning,
#' particularly suited for high dimensional data.
#' Classification, regression, and survival forests are supported.
#' @references
#' Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for
#' high dimensional data in C++ and R. Journal of Statistical Software 77:1-17.
#' @seealso ranger
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("A", 20), rep("B", 20), rep("C", 20)))
#' rf.out.of.bag(x, y, imp_pvalues=FALSE)
#' rf.out.of.bag(x, y, imp_pvalues=TRUE)
#' y0<-factor(c(rep("old", 30), rep("young", 30)))
#' rf.out.of.bag(x, y0, imp_pvalues=FALSE)
#' y<- 1:60
#' rf.out.of.bag(x, y, imp_pvalues=FALSE)
#' \donttest{
#' x_ <- data.frame(rbind(t(rmultinom(7, 7500, rep(c(.201,.5,.02,.18,.099), 1000))),
#'             t(rmultinom(8, 750, rep(c(.201,.4,.12,.18,.099), 1000))),
#'             t(rmultinom(15, 750, rep(c(.011,.3,.22,.18,.289), 1000))),
#'             t(rmultinom(15, 750, rep(c(.091,.2,.32,.18,.209), 1000))),
#'             t(rmultinom(15, 750, rep(c(.001,.1,.42,.18,.299), 1000)))))
#' y_<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' rf.out.of.bag(x_, y_, imp_pvalues=FALSE)
#' }
#' @author Shi Huang
#' @export
"rf.out.of.bag" <-function(x, y, ntree=500, verbose=FALSE, sparse = FALSE, imp_pvalues=FALSE, ...){
  set.seed(123)
  x <- .as_feature_df(x)
  if(length(y)!=nrow(x)) stop("The target varible doesn't match the shape of train data matrix, please check!")
  if("y" %in% colnames(x)) stop("The feature table should not contain a column named 'y'.")
  data<-data.frame(y=y, x, check.names = FALSE)
  if(sparse){
    sparse_data <- Matrix::Matrix(data.matrix(data), sparse = sparse)
    rf.model <- ranger(dependent.variable.name="y", data=sparse_data, keep.inbag=TRUE, importance='permutation',
                               classification=is.factor(y), num.trees=ntree, verbose=verbose, probability = FALSE, ...)
  }else{
    rf.model<-ranger(dependent.variable.name="y", data=data, keep.inbag=TRUE, importance='permutation',
                             classification=is.factor(y), num.trees=ntree, verbose=verbose, probability = FALSE, ...)
  }
  result <- list()
  result$rf.model <- rf.model
  result$y <- y
  if(is.factor(y)){
    if(sparse){
      y_numeric<-as.factor(sparse_data[,'y'])
      result$predicted <- factor(rf.model$predictions,levels=levels(y_numeric)); levels(result$predicted)<-levels(y)
    }else{
      result$predicted <- factor(rf.model$predictions,levels=levels(y))
    }
    result$probabilities <- get.oob.probability.from.forest(rf.model, x); colnames(result$probabilities)<-levels(result$y)
    result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
    result$errs <- mean(result$predicted != result$y)
  }else{
    result$predicted <- rf.model$predictions
    reg_perf<-get.reg.oob.performance(rf.model, y)
    result<-append(result, reg_perf)
    cat("Mean squared residuals: ", reg_perf$MSE, "\n")
    cat("Mean absolute error: ", reg_perf$MAE, "\n")
    cat("pseudo R-squared (%explained variance): ", reg_perf$R_squared, "\n")
  }
  if(imp_pvalues==FALSE){
    result$importances <- rf.model$variable.importance
    result$importances[colSums(x)==0]<-NA
  }else{
      result$importances <- .altmann_pvalues(rf.model, data)
  }
  result$params <- list(ntree=ntree)
  result$error.type <- "oob"
  class(result) <- "rf.out.of.bag"
  return(result)
}

#' @title balanced.folds
#' @description Get balanced folds where each fold has close to overall class ratio
#' @param y A response vector. If a factor, classification is assumed, otherwise regression is assumed.
#' @param nfolds The number of folds in the cross validation.
#' @return folds
#' @examples
#' y<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' print(balanced.folds(y, nfolds=3))
#' y<- 1:60
#' print(balanced.folds(y, nfolds=3))
#' @author Shi Huang
#' @export
balanced.folds <- function(y, nfolds=3){
    folds = rep(0, length(y))
    if(is.factor(y)){
      classes<-levels(y)
    }else{
      y<-factor(findInterval(y, quantile(y, seq(0, 1, by=1/nfolds), type=5), rightmost.closed=TRUE))
      classes<-levels(y)
    }
    # size of each class
    Nk = table(y)
    # -1 or nfolds = len(y) means leave-one-out
    if (nfolds == -1 || nfolds == length(y)){
        invisible(1:length(y))
    }else{
    # Can't have more folds than there are items per class
    nfolds = min(nfolds, max(Nk))
    # Assign folds evenly within each class, then shuffle within each class
        for (k in 1:length(classes)){
            ixs <- which(y==classes[k])
            folds_k <- rep(1:nfolds, ceiling(length(ixs) / nfolds))
            folds_k <- folds_k[1:length(ixs)]
            folds_k <- sample(folds_k)
            folds[ixs] = folds_k
        }
        invisible(folds)
    }
}

#' @title group.folds
#' @description Get folds for group k-fold cross-validation: all samples from the same group
#' (e.g., the same individual in a longitudinal study, the same family, or the same study)
#' are assigned to the same fold, preventing information leakage from correlated samples.
#' Groups are assigned greedily to the currently smallest fold so that fold sizes are balanced.
#' @param groups A vector or factor indicating the group of each sample.
#' @param nfolds The number of folds. It cannot exceed the number of groups.
#' @return An integer vector of fold assignments.
#' @examples
#' g <- factor(rep(paste0("subject", 1:12), each=5))
#' folds <- group.folds(g, nfolds=3)
#' table(folds, g)
#' @export
group.folds <- function(groups, nfolds=5){
  groups <- factor(groups)
  if(nlevels(groups) < 2) stop("At least two groups are required for group k-fold cross-validation.")
  nfolds <- min(nfolds, nlevels(groups))
  group_sizes <- table(groups)
  group_sizes <- group_sizes[sample.int(length(group_sizes))] # random tie-breaking
  group_sizes <- group_sizes[order(group_sizes, decreasing=TRUE)]
  fold_sizes <- rep(0, nfolds)
  group_fold <- setNames(integer(length(group_sizes)), names(group_sizes))
  for(g in names(group_sizes)){
    k <- which.min(fold_sizes)
    group_fold[g] <- k
    fold_sizes[k] <- fold_sizes[k] + group_sizes[[g]]
  }
  as.integer(group_fold[as.character(groups)])
}

#' @title rf.cross.validation
#' @description It runs standard random forests with n-folds cross-validation error estimation for both classification and regression using rf.out.of.bag.
#' @param x Training data: data.matrix or data.frame.
#' @param y A response vector. If a factor, classification is assumed, otherwise regression is assumed.
#' @param ntree The number of trees.
#' @param nfolds The number of folds in the cross validation. If nfolds > length(y)
#' or nfolds==-1, uses leave-one-out cross-validation. If nfolds was a factor, it means customized folds (e.g., leave-one-group-out cv) were set for CV.
#' @param sparse A boolean value indicates if the input matrix transformed into sparse matrix for rf modeling.
#' @param verbose A boolean value indicates if showing computation status and estimated runtime.
#' @param imp_pvalues If compute both importance score and pvalue for each feature.
#' @param groups An optional vector indicating the group (e.g., subject ID) of each sample.
#' If provided (and nfolds is a number), group k-fold cross-validation is performed using \code{group.folds}.
#' @param ... Other parameters applicable to `ranger`.
#' @return Object of class \code{rf.cross.validation} with elements including a \code{ranger} object and mutiple metrics for model evaluations.
#' @seealso ranger
#' @examples
#'
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' s<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' y<-factor(rep(c("Y", "N"), 30))
#' y0<-factor(c(rep("A", 10), rep("B", 30), rep("C", 5), rep("D", 15)))
#' rf.cross.validation(x, y, imp_pvalues=FALSE)
#' rf.cross.validation(x, y0, imp_pvalues=FALSE)
#' y_n<- 1:60
#' rf.cross.validation(x, y_n, nfolds=5, imp_pvalues=FALSE)
#' # group k-fold cv: samples from the same subject are never split across folds
#' subject <- factor(rep(1:20, each=3))
#' rf.cross.validation(x, y, nfolds=5, groups=subject)
#' \donttest{
#' # the permutation-based importance p values (Altmann's method) take much longer
#' rf.cross.validation(x, y, imp_pvalues=TRUE)
#' rf.cross.validation(x, y_n, nfolds=5, imp_pvalues=TRUE)
#' # when nfolds is a factor, it actually run a leave-one-group-out cv
#' rf.cross.validation(x, y, nfolds=s, imp_pvalues=TRUE)
#' rf.cross.validation(x, y_n, nfolds=s, imp_pvalues=FALSE)
#' }
#' @author Shi Huang
#' @export
"rf.cross.validation" <- function(x, y, nfolds=3, ntree=500, verbose=FALSE, sparse = FALSE, imp_pvalues=FALSE, groups=NULL, ...){
  x <- .as_feature_df(x)
  if(length(y)!=nrow(x)) stop("The target varible doesn't match the shape of train data matrix, please check!")
  if("y" %in% colnames(x)) stop("The feature table should not contain a column named 'y'.")
  if(is.factor(nfolds)){
    if(length(nfolds)!=length(y)) stop("The customized folds (a factor) should have the same length as y.")
    folds = factor(nfolds) # it means the customized folds were set
  }else if(!is.null(groups)){
    if(length(groups)!=length(y)) stop("The length of groups should match the length of y.")
    folds <- group.folds(groups, nfolds=nfolds)
  }else{
    if(nfolds==-1 || nfolds > length(y)) nfolds <- length(y)
    folds <- balanced.folds(y, nfolds=nfolds)
  }
  fold_names <- as.character(sort(unique(folds)))
  nfolds <- length(fold_names)
  if(nfolds < length(y) && any(table(folds)<3)) stop("Less than 3 samples in at least one fold!\n")
  result <- list()
  result$rf.model<-list()
  reg_metrics <- c("MSE", "RMSE", "nRMSE", "MAE", "MAPE", "MASE", "Spearman_rho", "R_squared", "Adj_R_squared")
  if(is.factor(y)){
    result$y <- y <- factor(y)
    result$predicted <- y
    result$probabilities <- matrix(0, nrow=length(y), ncol=length(levels(y)))
    rownames(result$probabilities) <- rownames(x)
    colnames(result$probabilities) <- levels(y)
    result$errs <- setNames(numeric(nfolds), fold_names)
  }else{
    result$y <- y
    result$predicted <- y
    for(metric in reg_metrics) result[[metric]] <- setNames(rep(NA_real_, nfolds), fold_names)
  }
  result$importances <- matrix(0, nrow=ncol(x), ncol=nfolds, dimnames=list(colnames(x), fold_names))
  if(imp_pvalues) result$importance_pvalues <- result$importances
  # K-fold cross-validation
  for(fold in fold_names){
    if(verbose) cat("The # of folds: ", fold, "\n")
    foldix <- which(as.character(folds)==fold)
    if(is.factor(y)) y_tr<-factor(result$y[-foldix]) else y_tr<-result$y[-foldix]
    data<-data.frame(y=y_tr, x[-foldix, , drop=FALSE], check.names = FALSE)
    if(sparse){
      sparse_data <- Matrix::Matrix(data.matrix(data), sparse = TRUE)
      model <- ranger(dependent.variable.name="y", data=sparse_data, classification=is.factor(y_tr),
                      keep.inbag=FALSE, importance='permutation', verbose=verbose, num.trees=ntree, ...)
    }else{
      model <- ranger(dependent.variable.name="y", data=data, keep.inbag=FALSE, importance='permutation',
                      classification=is.factor(y_tr), num.trees=ntree, verbose=verbose, ...)
    }
    result$rf.model[[fold]] <- model
    newx <- x[foldix, , drop=FALSE]
    predicted_foldix<-predict(model, newx)$predictions
    if(is.factor(y)){
      if(sparse){
        y_numeric<-as.factor(sparse_data[,'y'])
        predicted_foldix <- factor(predicted_foldix, levels=levels(y_numeric))
        levels(predicted_foldix) <- levels(y_tr)
      }
      result$predicted[foldix] <- factor(as.character(predicted_foldix), levels=levels(y))
      probs <- get.predict.probability.from.forest(model, newx)
      if(sparse) colnames(probs) <- levels(y_tr)
      result$probabilities[foldix, colnames(probs)] <- probs
      result$errs[fold] <- mean(result$predicted[foldix] != result$y[foldix])
      if(verbose) cat("Error rate: ", result$errs[fold], "\n")
    }else{
      result$predicted[foldix] <- predicted_foldix
      reg_perf<-get.reg.performance(predicted_foldix, result$y[foldix], n_features=ncol(x))
      for(metric in reg_metrics) result[[metric]][fold] <- reg_perf[[metric]]
      if(verbose){
        cat("Mean squared residuals: ", reg_perf$MSE, "\n")
        cat("Mean absolute error: ", reg_perf$MAE, "\n")
        cat("pseudo R-squared (%explained variance): ", reg_perf$R_squared, "\n")
      }
    }
    # importance scores (and p values) are computed from the model trained in this fold
    if(imp_pvalues){
      imp<-.altmann_pvalues(model, data)
      result$importances[, fold] <- imp[, 1]
      result$importance_pvalues[, fold] <- imp[, 2]
    }else{
      result$importances[, fold] <- model$variable.importance
    }
  }
  if(is.factor(y)) result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
  result$params <- list(ntree=ntree, nfolds=nfolds)
  result$error.type <- "cv"
  result$nfolds <- nfolds
  result$folds <- folds
  class(result) <- "rf.cross.validation"
  return(result)
}

.class_codes <- function(model){
  if(!is.null(model$forest$levels)) seq_along(model$forest$levels) else seq_along(model$forest$class.values)
}

.class_names <- function(model){
  if(!is.null(model$forest$levels)) model$forest$levels else as.character(.class_codes(model))
}

#' @title get.oob.probability.from.forest
#' @description Get probability of each class using only out-of-bag predictions from RF
#' @param model A \code{ranger} classification model trained with \code{keep.inbag=TRUE}.
#' @param x The training data.
#' @return A table of probabilities \code{probs}
#' @author Shi Huang
#' @export
"get.oob.probability.from.forest" <- function(model, x){
    # get aggregated class votes for each sample using only OOB trees
    votes <- get.oob.votes.from.forest(model,x)
    # convert to probs
    probs <- votes/rowSums(votes)
    rownames(probs) <- rownames(x)
    return(invisible(probs))
}

#' @title get.oob.votes.from.forest
#' @description get votes for each class using only out-of-bag predictions from RF
#' @param model A \code{ranger} classification model trained with \code{keep.inbag=TRUE}.
#' @param x The training data.
#' @return votes
#' @author Shi Huang
#' @export
"get.oob.votes.from.forest" <- function(model, x){
    if(is.null(model$inbag.counts)) stop("OOB votes require a ranger model trained with keep.inbag=TRUE.")
    rf.pred <- predict(model, data.frame(x, check.names = FALSE), type="response", predict.all=TRUE)
    preds <- matrix(rf.pred$predictions, nrow=nrow(x))
    # a samples-by-trees matrix indicating which trees did not use each sample
    outofbag <- do.call(cbind, model$inbag.counts) == 0
    votes <- vapply(.class_codes(model), function(k) rowSums(preds == k & outofbag), numeric(nrow(x)))
    votes <- matrix(votes, nrow=nrow(x), dimnames=list(rownames(x), .class_names(model)))
    return(invisible(votes))
}

#' @title get.predict.probability.from.forest
#' @description to get votes from the forest in the external test. Last update: Apr. 7, 2019
#' @param model A \code{ranger} classification model.
#' @param newx The external test data.
#' @return \code{probs}
#' @author Shi Huang
#' @export
"get.predict.probability.from.forest" <- function(model, newx){
  rf.pred <- predict(model, data.frame(newx, check.names = FALSE), type="response", predict.all=TRUE)
  preds <- matrix(rf.pred$predictions, nrow=nrow(newx))
  votes <- vapply(.class_codes(model), function(k) rowSums(preds == k), numeric(nrow(newx)))
  votes <- matrix(votes, nrow=nrow(newx), dimnames=list(rownames(newx), .class_names(model)))
  probs <- votes/rowSums(votes)
  return(invisible(probs))
}

#' @title get.reg.performance
#' @description Get regression performance
#' @param pred_y The predicted y in the numeric format.
#' @param y A response vector for regression.
#' @param n_features The number of features.
#' @return A list of regression performance metrics.
#' @details The adjusted R squared is NA when the number of features is not smaller than the number of samples minus one.
#' @author Shi Huang
#' @export
"get.reg.performance" <- function(pred_y, y, n_features=NA){
  pred_y<-as.numeric(pred_y)
  y<-as.numeric(y)
  mean_y<-mean(y)
  MSE <- mean((y-pred_y)^2)
  RMSE <- sqrt(MSE)
  nRMSE <- RMSE/mean_y # <https://en.wikipedia.org/wiki/Root-mean-square_deviation#Normalized_root-mean-square_deviation>
  MAE <- mean(abs(y-pred_y))
  MAPE <- mean(abs(y-pred_y)/abs(y)) # <https://en.wikipedia.org/wiki/Mean_absolute_percentage_error>
  MASE <- MAE/mean(abs(diff(y))) # <https://en.wikipedia.org/wiki/Mean_absolute_scaled_error>
  Spearman_rho <- if(length(y) > 1) suppressWarnings(stats::cor(y, pred_y, method = "spearman")) else NA
  R2 <- function(y, pred_y){1 - (sum((y-pred_y)^2) / sum((y-mean(y))^2))}
  R_squared <- R2(y, pred_y)
  adj.R2<-function(y, pred_y, k){
    n <- length(y)
    if(n-k-1 <= 0) return(NA)
    1-(1-R2(y, pred_y))*(n-1)/(n-k-1) # k is # of predictors
  }
  if(!is.na(n_features)){
    Adj_R_squared <- adj.R2(y, pred_y, k = n_features)
  }else{
    Adj_R_squared <- NA
  }
  perf<-list()
  perf$MSE<-MSE
  perf$RMSE<-RMSE
  perf$nRMSE<-nRMSE
  perf$MAE<-MAE
  perf$MAPE<-MAPE
  perf$MASE<-MASE
  perf$Spearman_rho <- Spearman_rho
  perf$R_squared<-R_squared
  perf$Adj_R_squared<-Adj_R_squared
  return(invisible(perf))
}

#' @title get.reg.oob.performance
#' @description Get oob performance from forests for training data. Last update: Apr. 7, 2019
#' @param model A object of class \code{rf.out.of.bag} with elements.
#' @param y The responsive variable for regression in training data.
#' @return \code{perf}
#' @author Shi Huang
#' @export
"get.reg.oob.performance" <- function(model, y){
  pred_y<-as.numeric(model$predictions)
  y<-as.numeric(y)
  perf<-get.reg.performance(pred_y, y, n_features=model$num.independent.variables)
  return(invisible(perf))
}


#' @title get.reg.predict.performance
#' @description Get prediction performance from forests for the external test data. Last update: Apr. 7, 2019
#' @param model A object of random forest out-of-bag regression model.
#' @param newx The external test data.
#' @param newy The dependent variable in the external test data.
#' @return \code{perf}
#' @author Shi Huang
#' @export
"get.reg.predict.performance" <- function(model, newx, newy){
  rf.pred <- predict(model, data.frame(newx, check.names = FALSE))
  pred_y <- as.numeric(rf.pred$predictions)
  y <- as.numeric(newy)
  perf<-get.reg.performance(pred_y, y, n_features=ncol(newx))
  return(invisible(perf))
}

#' @title get.auroc
#' @description calculates the area under the ROC curve
#' @param y A binary factor vector indicates observed values.
#' @param predictor A predictor. For example, one column in the probability table, indicating
#'            the probability of a given observation being the positive class in the factor y.
#' @param positive_class A class of the factor y.
#' @return The auroc value.
#' @examples
#' y<-factor(c(rep("A", 31), rep("B", 29)))
#' y1<-factor(c(rep("A", 18), rep("B", 20), rep("C", 22)))
#' y0<-factor(rep("A", 60))
#' pred <- c(runif(30, 0.5, 0.9), runif(30, 0, 0.6))
#' prob <-data.frame(A=pred, B=1-pred)
#' positive_class="A"
#' get.auroc(predictor=prob[, positive_class], y, positive_class="A")
#' get.auroc(predictor=prob[, positive_class], y, positive_class="B")
#' get.auroc(predictor=prob[, positive_class], y0, positive_class="A")
#' get.auroc(predictor=prob[, positive_class], y1, positive_class="A")
#' @author Shi Huang
#' @export
get.auroc <- function(predictor, y, positive_class) {
  if(nlevels(factor(y))==1){
    cat("All y values are identical!\n")
    return(NA)
  }
  # the normalized Mann-Whitney U statistic equals the area under the ROC curve (ties count 1/2)
  keep <- !is.na(y) & !is.na(predictor)
  is_pos <- as.character(y[keep])==positive_class
  n_pos <- sum(is_pos)
  n_neg <- sum(!is_pos)
  if(n_pos==0 || n_neg==0) return(NA)
  ranks <- rank(predictor[keep])
  auroc <- (sum(ranks[is_pos]) - n_pos*(n_pos+1)/2)/(n_pos*n_neg)
  return(auroc)
}

#' @title get.auprc
#' @description calculates the area under Precision-recall curve (AUPRC).
#' @param y A binary factor vector indicates observed values.
#' @param predictor A predictor. For example, one column in the probability table, indicating
#'            the probability of a given observation being the positive class in the factor y.
#' @param positive_class A class of the factor y.
#' @details It’s a bit trickier to interpret AUPRC than it is to interpret AUROC.
#' The AUPRC of a random classifier is equal to the fraction of positives (Saito et al.),
#' where the fraction of positives is calculated as (# positive examples / total # examples).
#' That means that _different_ classes have _different_ AUPRC baselines. A class with 12% positives
#' has a baseline AUPRC of 0.12, so obtaining an AUPRC of 0.40 on this class is great. However
#' a class with 98% positives has a baseline AUPRC of 0.98, which means that obtaining an AUPRC
#' of 0.40 on this class is bad.
#' @return auprc
#' @examples
#' y<-factor(c(rep("A", 10), rep("B", 50)))
#' y1<-factor(c(rep("A", 18), rep("B", 20), rep("C", 22)))
#' y0<-factor(rep("A", 60))
#' pred <- c(runif(10, 0.4, 0.9), runif(50, 0, 0.6))
#' prob <-data.frame(A=pred, B=1-pred)
#' positive_class="A"
#' get.auprc(predictor=prob[, positive_class], y, positive_class="A")
#' get.auprc(predictor=prob[, positive_class], y, positive_class="B")
#' get.auprc(predictor=prob[, positive_class], y0, positive_class="A")
#' get.auprc(predictor=prob[, positive_class], y1, positive_class="A")
#' @author Shi Huang
#' @export
get.auprc<- function(predictor, y, positive_class){
  if(nlevels(factor(y))==1){
    cat("All y values are identical!\n")
    auprc <- NA # if y contains all identical values, no threshold can be used for auprc calculation
  }else{
    df<-data.frame(y, predictor)
    prob_pos<-df[df$y==positive_class, "predictor"]
    prob_neg<-df[df$y!=positive_class,"predictor"]
    pr<-pr.curve(prob_pos, prob_neg, curve=FALSE, rand.compute = TRUE)
    auprc<-pr$auc.integral
    if(pr$auc.integral==pr$rand$auc.integral) cat("Note: ", auprc, "is the auprc of a random classifier!\n")
  }
  auprc
}


#' @title get.mislabel.scores
#' @description Get mislabelled scores of samples
#' @param y A response vector. If a factor, classification is assumed, otherwise regression is assumed.
#' @param y.prob The probability output from the RF model.
#' @return mislabelled scores
#' @examples
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' rf_model<-rf.out.of.bag(x, y)
#' get.mislabel.scores(rf_model$y, rf_model$probabilities)
#' @author Shi Huang
#' @export
"get.mislabel.scores" <- function(y, y.prob){
  result <- matrix(0,nrow=length(y),ncol=3)
  # get matrices containing only p(other classes), and containing only p(class)
  mm <- model.matrix(~0 + y)
  y.prob.other.max <- apply(y.prob * (1-mm),1,max)
  y.prob.alleged <- apply(y.prob * mm, 1, max)
  result <- cbind(y.prob.alleged, y.prob.other.max, y.prob.alleged - y.prob.other.max)
  rownames(result) <- rownames(y.prob)
  colnames(result) <- c('P(alleged label)','P(second best)','P(alleged label)-P(second best)')
  return(result)
}

# Classification performance of predicted labels/probabilities against the observed classes.
.clf_perf_row <- function(predicted, y, probabilities, positive_class){
  out <- setNames(rep(NA_real_, length(.clf_metric_names)), .clf_metric_names)
  if(!positive_class %in% levels(y) || nlevels(y) < 2) return(out)
  predicted <- factor(as.character(predicted), levels=levels(y))
  conf <- caret::confusionMatrix(data=predicted, reference=y, positive=positive_class)
  by_class <- conf$byClass
  if(is.matrix(by_class)) by_class <- by_class[paste0("Class: ", positive_class), ]
  out[] <- c(conf$overall[["Accuracy"]],
             get.auroc(probabilities[, positive_class], y, positive_class),
             get.auprc(probabilities[, positive_class], y, positive_class),
             conf$overall[["Kappa"]], by_class)
  out
}

# Mean importance score of each feature (averaged over folds for rf.cross.validation).
.feature_importances <- function(rf_model){
  imps <- rf_model$importances
  if(inherits(rf_model, "rf.cross.validation")) return(rowMeans(imps, na.rm=TRUE))
  if(is.matrix(imps)) return(setNames(imps[, 1], rownames(imps)))
  imps
}

# Feature ranks (1 = most important); for rf.cross.validation, the median rank over folds.
.feature_rank <- function(rf_model){
  if(inherits(rf_model, "rf.cross.validation")){
    rank_mat <- apply(rf_model$importances, 2, function(x){rank(-x, na.last = "keep")})
    rank(apply(rank_mat, 1, median), na.last = "keep", ties.method = "first")
  }else if(inherits(rf_model, "rf.out.of.bag")){
    rank(-.feature_importances(rf_model), na.last = "keep", ties.method = "first")
  }else{
    stop("The class of input rf model should be rf.out.of.bag or rf.cross.validation.")
  }
}

.rfe_sizes <- function(n_total, n_features=NULL){
  if(is.null(n_features)) n_features <- 2^(1:10)
  n_features <- sort(unique(as.integer(n_features)))
  c(n_features[n_features >= 1 & n_features < n_total], n_total)
}

.rfe_select <- function(top_n_perf, metric, higher_better, tolerance){
  values <- top_n_perf[[metric]]
  if(all(is.na(values))) stop("All performance values of ", metric, " are NA.")
  best_idx <- if(higher_better) which.max(values) else which.min(values)
  best <- values[best_idx]
  within <- if(higher_better) values >= best - tolerance*abs(best) else values <= best + tolerance*abs(best)
  list(best_n_features=top_n_perf$n_features[best_idx],
       optimal_n_features=top_n_perf$n_features[min(which(within))])
}

.rfe <- function(x, y, folds, rf_model, n_features, metric, higher_better, tolerance, recursive,
                 ntree, verbose, perf_fun, class_name, extra=list(), ...){
  n_total <- ncol(x)
  sizes <- .rfe_sizes(n_total, n_features)
  full_rank <- .feature_rank(rf_model)
  key <- as.character(n_total)
  models <- setNames(list(rf_model), key)
  feature_sets <- setNames(list(colnames(x)), key)
  perf_list <- setNames(list(perf_fun(rf_model)), key)
  prev_rank <- full_rank
  for(n in rev(sizes[sizes < n_total])){
    ranks <- if(recursive) prev_rank else full_rank
    ranks <- ranks[!is.na(ranks)]
    keep <- names(ranks)[order(ranks)][seq_len(min(n, length(ranks)))]
    if(verbose) cat("Feature selection: modeling with the top", length(keep), "features\n")
    fit <- function() rf.cross.validation(x[, keep, drop=FALSE], y, nfolds=folds, ntree=ntree, verbose=verbose, ...)
    model <- if(verbose) fit() else .quietly(fit())
    key <- as.character(n)
    models[[key]] <- model
    feature_sets[[key]] <- keep
    perf_list[[key]] <- perf_fun(model)
    prev_rank <- .feature_rank(model)
  }
  top_n_perf <- data.frame(n_features=as.numeric(names(perf_list)), do.call(rbind, perf_list), check.names=FALSE)
  top_n_perf <- top_n_perf[order(top_n_perf$n_features), ]
  rownames(top_n_perf) <- top_n_perf$n_features
  selection <- .rfe_select(top_n_perf, metric, higher_better, tolerance)
  ordered_keys <- as.character(top_n_perf$n_features)
  res <- c(list(top_n_perf=top_n_perf,
                best_n_features=selection$best_n_features,
                optimal_n_features=selection$optimal_n_features,
                selected_features=feature_sets[[as.character(selection$optimal_n_features)]],
                feature_sets=feature_sets[ordered_keys],
                top_n_rf=models[ordered_keys],
                importances=.feature_importances(rf_model),
                metric=metric, tolerance=tolerance, recursive=recursive, folds=folds),
           extra)
  class(res) <- class_name
  res
}

#' @title rf_clf.rfe
#' @description Recursive feature elimination (RFE) for random forest classification.
#' RF models are built with cross-validation on successively reduced feature sets
#' (by default the top 2, 4, 8, 16, 32, 64, 128, 256, 512 and 1024 features by importance, plus all features).
#' With \code{recursive=TRUE}, the importance ranking is re-computed from the model at each step;
#' otherwise the ranking of the model with all features is reused.
#' The same cross-validation folds are used at every step, so that feature sets are compared on identical splits.
#' @param x The data frame or data matrix for model training.
#' @param y A factor related to the responsive vector for training data.
#' @param nfolds The number of folds in the cross-validation, or a factor defining customized folds (e.g., leave-one-group-out).
#' @param n_features A numeric vector of feature set sizes. By default powers of two up to 1024.
#' @param metric The classification performance metric used for selection: one of "AUROC", "AUPRC", "Accuracy", "Kappa", "F1".
#' @param positive_class A class of the y.
#' @param tolerance The relative tolerance: the optimal parsimonious model is the smallest feature set
#' whose performance is within \code{tolerance} (default 1\%) of the best performance.
#' @param recursive A boolean value indicating if the feature ranking is re-computed at each step.
#' @param ntree The number of trees.
#' @param verbose A boolean value indicating if showing computation status.
#' @param rf_model An optional pre-computed \code{rf.cross.validation} or \code{rf.out.of.bag} model on all features.
#' @param ... Other parameters applicable to \code{ranger}.
#' @return An object of class \code{rf_clf.rfe} including the performance table (\code{top_n_perf}),
#' the best and optimal number of features, the \code{selected_features}, and the models at each step (\code{top_n_rf}).
#' @details As the features are ranked using all samples, the cross-validated performance of the selected
#' feature sets is optimistic. The generalizability of the selected biomarkers should be assessed in
#' independent datasets, e.g., using \code{rf_clf.cross_appl} or \code{rf_clf.lodo}.
#' @seealso plot_clf_feature_selection plot_topN_imp_scores id_robust_markers
#' @examples
#' set.seed(123)
#' n_features <- 100
#' prob_vec <- gtools::rdirichlet(5, sample(n_features))
#' x <- data.frame(rbind(t(rmultinom(7, 7*n_features, prob_vec[1, ])),
#'             t(rmultinom(8, 8*n_features, prob_vec[2, ])),
#'             t(rmultinom(15, 15*n_features, prob_vec[3, ])),
#'             t(rmultinom(15, 15*n_features, prob_vec[4, ])),
#'             t(rmultinom(15, 15*n_features, prob_vec[5, ]))))
#' y<-factor(c(rep("A", 30), rep("C", 30)))
#' rfe <- rf_clf.rfe(x, y, nfolds=5, metric="AUROC", ntree=200)
#' rfe$top_n_perf
#' rfe$selected_features
#' @export
rf_clf.rfe <- function(x, y, nfolds=5, n_features=NULL, metric="AUROC", positive_class=NA,
                       tolerance=0.01, recursive=TRUE, ntree=500, verbose=FALSE, rf_model=NULL, ...){
  metrics <- c("AUROC", "AUPRC", "Accuracy", "Kappa", "F1")
  if(!metric %in% metrics) stop("metric should be one of: ", paste(metrics, collapse=", "))
  x <- .as_feature_df(x)
  y <- factor(y)
  positive_class <- ifelse(is.na(positive_class), levels(y)[1], positive_class)
  folds <- if(is.factor(nfolds)) nfolds else factor(balanced.folds(y, nfolds=nfolds))
  if(is.null(rf_model)){
    fit <- function() rf.cross.validation(x, y, nfolds=folds, ntree=ntree, verbose=verbose, ...)
    rf_model <- if(verbose) fit() else .quietly(fit())
  }
  perf_fun <- function(model) .clf_perf_row(model$predicted, factor(model$y), model$probabilities, positive_class)[metrics]
  .rfe(x, y, folds, rf_model, n_features, metric, higher_better=TRUE, tolerance, recursive,
       ntree, verbose, perf_fun, class_name="rf_clf.rfe", extra=list(positive_class=positive_class), ...)
}

#' @title rf_reg.rfe
#' @description Recursive feature elimination (RFE) for random forest regression.
#' See \code{rf_clf.rfe} for details of the procedure.
#' @param x The data frame or data matrix for model training.
#' @param y The numeric values for labeling data.
#' @param nfolds The number of folds in the cross-validation, or a factor defining customized folds.
#' @param n_features A numeric vector of feature set sizes. By default powers of two up to 1024.
#' @param metric The regression performance metric used for selection: one of "MAE", "RMSE", "MSE", "MAPE", "Spearman_rho", "R_squared".
#' @param tolerance The relative tolerance for choosing the smallest feature set close to the best performance.
#' @param recursive A boolean value indicating if the feature ranking is re-computed at each step.
#' @param ntree The number of trees.
#' @param verbose A boolean value indicating if showing computation status.
#' @param rf_model An optional pre-computed \code{rf.cross.validation} or \code{rf.out.of.bag} model on all features.
#' @param ... Other parameters applicable to \code{ranger}.
#' @return An object of class \code{rf_reg.rfe}.
#' @seealso plot_reg_feature_selection
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
#' rfe <- rf_reg.rfe(x, y, nfolds=5, metric="MAE", ntree=200)
#' rfe$top_n_perf
#' @export
rf_reg.rfe <- function(x, y, nfolds=5, n_features=NULL, metric="MAE", tolerance=0.01,
                       recursive=TRUE, ntree=500, verbose=FALSE, rf_model=NULL, ...){
  metrics <- c("MSE", "RMSE", "MAE", "MAPE", "Spearman_rho", "R_squared")
  if(!metric %in% metrics) stop("metric should be one of: ", paste(metrics, collapse=", "))
  if(!is.numeric(y)) stop("y should be numeric for regression.")
  x <- .as_feature_df(x)
  folds <- if(is.factor(nfolds)) nfolds else factor(balanced.folds(y, nfolds=nfolds))
  if(is.null(rf_model)){
    fit <- function() rf.cross.validation(x, y, nfolds=folds, ntree=ntree, verbose=verbose, ...)
    rf_model <- if(verbose) fit() else .quietly(fit())
  }
  perf_fun <- function(model) unlist(get.reg.performance(model$predicted, y)[metrics])
  .rfe(x, y, folds, rf_model, n_features, metric, higher_better=metric %in% c("Spearman_rho", "R_squared"),
       tolerance, recursive, ntree, verbose, perf_fun, class_name="rf_reg.rfe", ...)
}
