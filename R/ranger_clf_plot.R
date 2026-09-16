#' @importFrom stats median reorder
#' @importFrom pROC plot.roc ci.se
#' @importFrom grDevices dev.off hcl pdf
#' @importFrom graphics plot text
#' @importFrom rlang .data
#' @importFrom utils data write.table head
#' @importFrom ggplot2 ggplot aes xlab ylab theme_bw coord_flip ggsave geom_boxplot geom_histogram geom_violin geom_jitter geom_point geom_line geom_hline geom_vline xlim ylim scale_x_continuous scale_color_manual theme element_rect element_line element_text element_blank
#' @importFrom PRROC roc.curve pr.curve
#' @importFrom gridExtra arrangeGrob
#' @importFrom reshape2 melt
#' @importFrom gtools rdirichlet
#' @import viridis

.clf_res_classes <- c("rf_clf.by_datasets", "rf_clf.comps")

# ROC or precision-recall curves of all within-study models in an rf_clf.by_datasets/rf_clf.comps object.
.plot_clf_curves <- function(res_list, type=c("ROC", "PRC"), positive_class=NA, prefix="train", outdir=NULL){
  type <- match.arg(type)
  if(is.na(positive_class)){
    positive_class <- if(!is.null(res_list$positive_class)) res_list$positive_class else levels(factor(res_list$y_list[[1]]))[1]
  }
  datasets <- if(!is.null(res_list$datasets)) res_list$datasets else names(res_list$rf_model_list)
  auc_name <- if(type=="ROC") "AUROC" else "AUPRC"
  curve_list <- lapply(seq_along(res_list$rf_model_list), function(i){
    model <- res_list$rf_model_list[[i]]
    y <- factor(model$y)
    if(!positive_class %in% levels(y) || nlevels(y) < 2) return(NULL)
    score <- model$probabilities[, positive_class]
    fg <- score[y==positive_class]
    bg <- score[y!=positive_class]
    if(type=="ROC"){
      crv <- PRROC::roc.curve(fg, bg, curve=TRUE)
      auc <- crv$auc
    }else{
      crv <- PRROC::pr.curve(fg, bg, curve=TRUE)
      auc <- crv$auc.integral
    }
    data.frame(x=crv$curve[, 1], y=crv$curve[, 2],
               dataset=sprintf("%s (%s = %.3f)", datasets[i], auc_name, auc))
  })
  curves <- do.call(rbind, curve_list)
  if(is.null(curves)) stop("The positive class '", positive_class, "' is absent from all datasets.")
  curves$dataset <- factor(curves$dataset, levels=unique(curves$dataset))
  p <- ggplot(curves, aes(x=.data$x, y=.data$y, color=.data$dataset)) +
    geom_path(linewidth=0.8) +
    xlab(if(type=="ROC") "False positive rate (1 - specificity)" else "Recall") +
    ylab(if(type=="ROC") "True positive rate (sensitivity)" else "Precision") +
    labs(color=NULL, title=paste("Positive class:", positive_class)) +
    coord_equal(xlim=c(0, 1), ylim=c(0, 1)) +
    theme_bw()
  if(type=="ROC") p <- p + geom_abline(intercept=0, slope=1, linetype="dashed", color="grey60")
  if(!is.null(outdir)){
    ggsave(filename=paste(outdir, prefix, ".rf_clf_", type, ".datasets.pdf", sep=""), plot=p, width=6, height=4.5)
  }
  p
}

#' @title plot_clf_pROC
#' @param y A factor of classes to be used as the true results
#' @param rf_clf_model A list object of a random forest model.
#' @param positive_class an optional character string for the factor level that corresponds to a "positive" result (if that makes sense for your data).
#' If there are only two factor levels, the first level will be used as the "positive" result.
#' @param prefix The prefix of data set.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("A", 30), rep("C", 30)))
#' rf_clf_model<-rf.out.of.bag(x, y)
#' plot_clf_pROC(y, rf_clf_model, positive_class="A", outdir=paste0(tempdir(), "/"))
#' @author Shi Huang
#' @export
plot_clf_pROC<-function(y, rf_clf_model, positive_class=NA, prefix="train", outdir=NULL){
  if(!inherits(rf_clf_model, c("rf.out.of.bag", "rf.cross.validation"))) stop("The rf_clf_model is not rf.out.of.bag/rf.cross.validation.")
  if(nlevels(y)!=2) stop("pROC only support for ROC analysis in the binary classification!")
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  predictor<-rf_clf_model$probabilities[, positive_class]
  # the positive class is the "case" level, and higher probabilities predict it (no automatic direction flipping)
  rocobj <- pROC::roc(response=y, predictor=predictor, levels=c(setdiff(levels(y), positive_class), positive_class),
                      direction="<", percent=TRUE, ci=TRUE) # the AUC with its CI
  ciobj <- pROC::ci.se(rocobj, specificities=seq(0, 100, 5)) # over a select set of specificities
  if(!is.null(outdir)){
  pdf(paste(outdir, prefix, ".rf_clf_pROC.ci.pdf",sep=""), width=4, height=4)
  pROC::plot.roc(rocobj, main="")
  plot(ciobj, type="shape", col="#1c61b6AA") # plot as a blue shape
  text(70,25, paste0(levels(y), collapse = " VS "), pos=4)
  text(70,15, paste("AUC = ",formatC(rocobj$auc,digits=2,format="f"),sep=""),pos=4)
  ci.lower<-formatC(rocobj$ci[1],digits=2,format="f")
  ci.upper<-formatC(rocobj$ci[3],digits=2,format="f")
  text(70,5, paste("95% CI: ",ci.lower,"-",ci.upper,sep=""),pos=4)
  dev.off()
  }
  result<-list()
  result$rocobj<-rocobj
  result$ciobj<-ciobj
  invisible(result)
}

#' @title plot_clf_PRC
#' @description Precision-recall curve of a classification model, or of all within-study models
#' when an \code{rf_clf.by_datasets}/\code{rf_clf.comps} object is given.
#' @param y A factor of classes to be used as the true results, or an object of class
#' \code{rf_clf.by_datasets} or \code{rf_clf.comps} to plot the curves of all within-study models in one figure.
#' @param rf_clf_model A list object of a random forest model.
#' @param positive_class an optional character string for the factor level that corresponds to a "positive" result (if that makes sense for your data).
#' If there are only two factor levels, the first level will be used as the "positive" result.
#' @param prefix The prefix of data set.
#' @param outdir The output directory.
#' @return A \code{PRROC} object (invisibly), or a ggplot object for multiple datasets.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' x0 <- data.frame(rbind(t(rmultinom(7, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(8, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289)))))
#' y<-factor(c(rep("A", 30), rep("B", 30)))
#' y0<-factor(c(rep("A", 5), rep("B", 55)))
#' rf_clf_model<-rf.out.of.bag(x, y)
#' plot_clf_PRC(y, rf_clf_model, positive_class="A")
#' plot_clf_PRC(y, rf_clf_model, positive_class="B")
#' rf_clf_model0<-rf.out.of.bag(x0, y0)
#' plot_clf_PRC(y0, rf_clf_model0, positive_class="A")
#' @author Shi Huang
#' @export
plot_clf_PRC<-function(y, rf_clf_model, positive_class=NA, prefix="train", outdir=NULL){
  if(inherits(y, .clf_res_classes)) return(.plot_clf_curves(y, type="PRC", positive_class=positive_class, prefix=prefix, outdir=outdir))
  if(!inherits(rf_clf_model, c("rf.out.of.bag", "rf.cross.validation"))) stop("The rf_clf_model is not rf.out.of.bag/rf.cross.validation.")
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  predictor<-rf_clf_model$probabilities[, positive_class]
  df<-data.frame(y, predictor)
  prob_pos<-df[df$y==positive_class, "predictor"]
  prob_neg<-df[df$y!=positive_class,"predictor"]
  pr<-pr.curve(prob_pos, prob_neg, curve=TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
  auprc<-round(pr$auc.integral, 3)
  rel_auprc<-round((pr$auc.integral-pr$min$auc.integral)/(pr$max$auc.integral-pr$min$auc.integral), 3)
  if(!is.null(outdir)){
    pdf(paste(outdir, prefix, ".rf_clf_PRC.pdf",sep=""), width=5, height=5)
    plot(pr, max.plot = TRUE, min.plot = TRUE, rand.plot = TRUE, fill.area = TRUE,
         auc.main = FALSE, main = paste("AUPRC=",auprc,"\n", "Relative AUPRC=", rel_auprc, sep=""))
    dev.off()
  }
  invisible(pr)
}

#' @title plot_clf_ROC
#' @description ROC curve of a classification model, or of all within-study models
#' when an \code{rf_clf.by_datasets}/\code{rf_clf.comps} object is given.
#' @param y A factor of classes to be used as the true results, or an object of class
#' \code{rf_clf.by_datasets} or \code{rf_clf.comps} to plot the ROC curves of all within-study models in one figure.
#' @param rf_clf_model A list object of a random forest model.
#' @param positive_class an optional character string for the factor level that corresponds to a "positive" result (if that makes sense for your data).
#' If there are only two factor levels, the first level will be used as the "positive" result.
#' @param prefix The prefix of data set.
#' @param outdir The output directory.
#' @return A \code{PRROC} object (invisibly), or a ggplot object for multiple datasets.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("A", 20), rep("B", 20), rep("C", 20)))
#' rf_clf_model<-rf.out.of.bag(x, y)
#' plot_clf_ROC(y, rf_clf_model, positive_class="A")
#' plot_clf_ROC(y, rf_clf_model, positive_class="B")
#' pred_df<-data.frame(y=y, prediction=rf_clf_model$probabilities[, "A"])
#' ggplot2::ggplot(pred_df, ggplot2::aes(prediction, fill = y)) +
#' ggplot2::geom_histogram(alpha = 0.5, position = 'identity', bins = 30)
#' # ROC curves of all within-study models
#' metadata <- data.frame(study=factor(rep(c("S1", "S2"), each=30)),
#'                        disease=factor(rep(rep(c("Case", "Control"), each=15), 2)))
#' res <- rf_clf.by_datasets(x, metadata, s_category="study", c_category="disease",
#'                           positive_class="Case", n_cores=1)
#' plot_clf_ROC(res)
#' @author Shi Huang
#' @export
plot_clf_ROC<-function(y, rf_clf_model, positive_class=NA, prefix="train", outdir=NULL){
  if(inherits(y, .clf_res_classes)) return(.plot_clf_curves(y, type="ROC", positive_class=positive_class, prefix=prefix, outdir=outdir))
  if(!inherits(rf_clf_model, c("rf.out.of.bag", "rf.cross.validation"))) stop("The rf_clf_model is not rf.out.of.bag/rf.cross.validation.")
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  predictor<-rf_clf_model$probabilities[, positive_class]
  df<-data.frame(y, predictor)
  prob_pos<-df[df$y==positive_class, "predictor"]
  prob_neg<-df[df$y!=positive_class,"predictor"]
  roc<-roc.curve(prob_pos, prob_neg, curve=TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
  if(!is.null(outdir)){
    pdf(paste(outdir, prefix, ".rf_clf_ROC.pdf",sep=""), width=5, height=5)
    plot(roc, rand.plot = TRUE, fill.area = TRUE)
    dev.off()
  }
  invisible(roc)
}



#' @title plot_clf_probabilities
#' @param y A factor of classes to be used as the true results
#' @param rf_clf_model The rf classification model from \code{rf.out.of.bag}
#' @param positive_class an optional character string for the factor level that corresponds to a "positive" result (if that makes sense for your data).
#' If there are only two factor levels, the first level will be used as the "positive" result.
#' @param prefix The prefix of data set.
#' @param outdir The output directory.
#' @return A ggplot object.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("1", 20), rep("A", 20), rep("C", 20)))
#' rf_clf_model<-rf.out.of.bag(x, y)
#' positive_class="1"
#' plot_clf_probabilities(y, rf_clf_model, positive_class)
#' @author Shi Huang
#' @export
plot_clf_probabilities<-function(y, rf_clf_model, positive_class=NA, prefix="train", outdir=NULL){
  if(!inherits(rf_clf_model, c("rf.out.of.bag", "rf.cross.validation"))) stop("The rf_clf_model is not rf.out.of.bag/rf.cross.validation.")
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  l<-levels(y); l_sorted<-sort(levels(y))
  Mycolor <- rep(c("#D55E00", "#0072B2"), length.out=length(l))
  if(identical(order(l), order(l_sorted))){
    Mycolor=Mycolor; l_ordered=l
  }else{Mycolor=rev(Mycolor); l_ordered=l_sorted}
  y_prob<-data.frame(y, predictor=rf_clf_model$probabilities[,positive_class])
  p<-ggplot(y_prob, aes(x=.data$y, y=.data$predictor)) +
    geom_violin()+
    geom_jitter(position=position_jitter(width=0.2), alpha=0.1) +
    geom_boxplot(outlier.shape = NA, width=0.4, alpha=0.01)+
    geom_hline(yintercept=0.5, linetype="dashed")+
    ylim(0, 1)+
    theme_bw()+
    xlab("") +
    ylab(paste("Probability of ", positive_class))+
    theme(legend.position="none")+
    theme(axis.line = element_line(color="black"),
          strip.background = element_rect(colour = "white"),
          panel.border = element_blank())+
    scale_color_manual(values = Mycolor, labels=l_ordered)
  if(!is.null(outdir)){
  ggsave(filename=paste(outdir,prefix,".probability_",positive_class,".boxplot.pdf",sep=""), plot=p, width=3, height=4)
  }
  p
}

#' @title plot_topN_imp_scores
#' @description Bar plot of the top N features ranked by random forest importance scores.
#' @param rf_model The rf model from \code{rf.out.of.bag} or \code{rf.cross.validation}, or the feature selection
#' result from \code{rf_clf.rfe}/\code{rf_reg.rfe} (the model of the step with \code{n_features} features is used).
#' @param topN A number indicating how many top important feature need to visualize in the barplot.
#' @param feature_md an optional data.frame including feature IDs and feature annotations.
#' @param feature_id_col The Feature_ID column in the feature metadata.
#' @param bar_color The column name in the feature metadata for coloring the bar plot.
#' @param plot_width plotting parameter.
#' @param plot_height plotting parameter.
#' @param outdir The output directory.
#' @param n_features For a feature selection result, the size of the feature set whose model is plotted.
#' By default, the optimal number of features.
#' @return A list including the importance table of all features ordered by rank (\code{imps_df}) and the plot.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("1", 20), rep("A", 20), rep("C", 20)))
#' feature_md <- data.frame(Feature_ID=c("X1", "X2", "X3", "X4", "X5"),
#'                          Taxon=c("AA", "Bd", "CC", "CC", "AA"))
#' rf_model<-rf.out.of.bag(x, y)
#' plot_topN_imp_scores(rf_model, topN=4)
#' rf_model<-rf.cross.validation(x, y, nfolds=5)
#' plot_topN_imp_scores(rf_model, topN=4, feature_md, outdir=NULL, bar_color="Taxon")
#' @author Shi Huang
#' @export
plot_topN_imp_scores<-function(rf_model, topN=4, feature_md=NULL, feature_id_col="Feature_ID", bar_color=NA,
                               outdir=NULL, plot_height=8, plot_width=5, n_features=NULL){
  if(inherits(rf_model, c("rf_clf.rfe", "rf_reg.rfe"))){
    key <- as.character(if(is.null(n_features)) rf_model$optimal_n_features else n_features)
    if(!key %in% names(rf_model$top_n_rf))
      stop("n_features should be one of: ", paste(names(rf_model$top_n_rf), collapse=", "))
    rf_model <- rf_model$top_n_rf[[key]]
  }
  if(!inherits(rf_model, c("rf.out.of.bag", "rf.cross.validation")))
    stop("The rf_model is not rf.out.of.bag/rf.cross.validation.")
  feature_rank <- .feature_rank(rf_model)
  if(inherits(rf_model, "rf.cross.validation")){
    imps_df<-data.frame(Feature_ID=rownames(rf_model$importances),
                     rf_model$importances,
                     Imps=apply(rf_model$importances, 1, median),
                     ImpsSD=apply(rf_model$importances, 1, sd),
                     ImpsSE=apply(rf_model$importances, 1, sd)/sqrt(ncol(rf_model$importances)),
                     Rank=feature_rank)
  }else{
    imps <- .feature_importances(rf_model)
    imps_df<-data.frame(Feature_ID=names(imps),
                    Imps=imps,
                    ImpsSD=0,
                    ImpsSE=0,
                    Rank=feature_rank)
  }
  # add feature metadata
  merge_feature_md<-function(df, feature_md, df_id_col=1, fmd_id_col=1){
    matched_idx<-which(feature_md[, fmd_id_col] %in% df[, df_id_col])
    uniq_features_len<-length(unique(df[, df_id_col]))
    if(uniq_features_len  > length(matched_idx)){
      warning("# of features has no matching IDs in the feature metadata file: ", uniq_features_len-length(matched_idx), "\n")
    }
    feature_md_matched<-feature_md[matched_idx, , drop=FALSE]
    out<-merge(df, feature_md_matched, by.x=df_id_col, by.y=fmd_id_col, all.x=TRUE)
    out
  }
  if(!is.null(feature_md)){
    imps_df<- merge_feature_md(imps_df, feature_md, df_id_col = "Feature_ID", fmd_id_col = feature_id_col)
  }
  imps_df<-imps_df[order(imps_df$Rank), ]
  rownames(imps_df)<-NULL
  top_n_imps_df<-head(imps_df[!is.na(imps_df$Rank), ], topN)
  # bar plot
  p <- ggplot(top_n_imps_df, aes(x=reorder(.data$Feature_ID, .data$Imps), y=.data$Imps))
  if(!is.na(bar_color)){
    p <- p + geom_bar(aes(fill=.data[[bar_color]]), stat="identity", alpha=0.7) +
      scale_fill_viridis(discrete=TRUE) + labs(fill=bar_color)
  }else{
    p <- p + geom_bar(stat="identity", alpha=0.5)
  }
  p <- p +
    ylab("RF importance score") +
    xlab("Feature ID") +
    theme_minimal() +
    coord_flip()
  if(!is.null(outdir)){
    ggsave(filename=paste(outdir, "Top_", topN,".imps.barplot.pdf",sep=""), plot=p, width=plot_width, height=plot_height)
  }

  res <- list()
  res$imps_df <- imps_df
  res$plot <- p
  res

}

.plot_feature_selection <- function(perf, metric, selection, y_label=metric){
  optimal_value <- perf[[metric]][perf$n_features==selection$optimal_n_features]
  ggplot(perf, aes(x=.data$n_features, y=.data[[metric]])) +
    xlab("# of features used")+
    ylab(y_label)+
    scale_x_continuous(trans = "log2", breaks=perf$n_features)+
    geom_point() + geom_line()+
    geom_vline(xintercept = selection$optimal_n_features, linetype="dashed", color="blue")+
    annotate(geom="text", x=selection$optimal_n_features, y=optimal_value, label=selection$optimal_n_features,
             color="blue", vjust=-1)+
    theme_bw()+
    theme(axis.line = element_line(color="black"),
          axis.title = element_text(size=18),
          strip.background = element_rect(colour = "white"),
          panel.border = element_blank())
}

#' @title plot_clf_feature_selection
#' @description Plot the classification performance against the gradually reduced number of features used in the modeling.
#' The dashed line marks the optimal parsimonious feature set: the smallest set whose performance is within
#' \code{tolerance} of the best performance.
#' @param x The data frame or data matrix for model training, or the result of \code{rf_clf.rfe}.
#' @param y A factor related to the responsive vector for training data.
#' @param nfolds The number of folds in the cross-validation for each feature set.
#' @param rf_clf_model An optional rf classification model from \code{rf.out.of.bag} or \code{rf.cross.validation}
#' used to rank the features.
#' @param positive_class A class of the y.
#' @param metric The classification performance metric applied.
#' If binary classification, this must be one of "AUROC", "AUPRC", "Accuracy", "Kappa", "F1".
#' If multi-class classification, this must be one of "Accuracy", "Kappa".
#' @param outdir The output directory.
#' @param tolerance The relative tolerance for choosing the optimal parsimonious feature set.
#' @param recursive A boolean value indicating if the feature ranking is re-computed at each step (see \code{rf_clf.rfe}).
#' @param ntree The number of trees.
#' @return A list including the performance table (\code{top_n_perf}), \code{best_n_features},
#' \code{optimal_n_features}, \code{selected_features}, the models (\code{top_n_rf}), the plot and the \code{rf_clf.rfe} object.
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
#' s<-factor(rep(c("B1", "B2", "B3", "B4"), 15))
#' rf_model<-rf.cross.validation(x, y, nfolds=5)
#' summ <- plot_clf_feature_selection(x, y, nfolds=5, rf_model, metric="AUROC", outdir=NULL)
#' summ$plot
#' summ <- plot_clf_feature_selection(x, y, nfolds=s, rf_model, metric="AUROC", outdir=NULL)
#' # recursive feature elimination
#' rfe <- rf_clf.rfe(x, y, nfolds=5, metric="AUROC", ntree=200)
#' plot_clf_feature_selection(rfe)$plot
#' @author Shi Huang
#' @export
plot_clf_feature_selection <- function(x, y, nfolds=5, rf_clf_model=NULL,
                                       metric="AUROC", positive_class=NA, outdir=NULL,
                                       tolerance=0.01, recursive=FALSE, ntree=500){
  if(inherits(x, "rf_clf.rfe")){
    rfe <- x
    if(missing(metric)) metric <- rfe$metric
  }else{
    rfe <- rf_clf.rfe(x, y, nfolds=nfolds, metric=metric, positive_class=positive_class,
                      tolerance=tolerance, recursive=recursive, ntree=ntree, rf_model=rf_clf_model)
  }
  top_n_perf <- rfe$top_n_perf
  if(!metric %in% colnames(top_n_perf)) stop("metric should be one of: ", paste(colnames(top_n_perf)[-1], collapse=", "))
  selection <- .rfe_select(top_n_perf, metric, higher_better=TRUE, tolerance=rfe$tolerance)
  p <- .plot_feature_selection(top_n_perf, metric, selection)
  if(!is.null(outdir)){
  ggsave(filename=paste(outdir,"train.",metric,"_VS_top_ranking_features.scatterplot.pdf",sep=""), plot=p, width=5, height=4)
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
