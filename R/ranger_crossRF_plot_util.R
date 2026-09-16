#' @importFrom corrplot corrplot
#' @importFrom grDevices colorRampPalette
#' @importFrom stats cor reorder hclust dist
#' @importFrom utils write.table
#' @importFrom plyr ldply
#' @importFrom rlang .data
#' @import ggplot2
#' @import viridis
#' @importFrom gridExtra arrangeGrob
#' @importFrom reshape2 melt

#' @title corr_datasets_by_imps
#' @description The correlation of importance scores between datasets, and the consistency of important
#' features across datasets: features ranked in the top N of all (or at least \code{min_datasets}) datasets
#' are candidate universal markers, whereas features top-ranked in fewer datasets are dataset-specific.
#' @param feature_imps_list a list of feature importance scores from the output of \code{rf_reg.by_datasets} or \code{rf_clf.by_datasets},
#' or the \code{rf_clf.by_datasets}/\code{rf_reg.by_datasets}/\code{rf_clf.comps} object itself.
#' Importance scores are matched across datasets by feature names.
#' @param ranked if transform importance scores into rank for correlation analysis
#' @param plot if plot the correlation matrix
#' @param top_n An optional number: if provided, the features ranked in the top N of each dataset are
#' summarized in \code{feature_consistency}.
#' @param min_datasets The minimum number of datasets in which a feature is top-ranked to be a universal marker.
#' By default, all datasets.
#' @return A list including the correlation matrix (\code{corr_mat}), its p values (\code{p_mat}),
#' the importance (\code{imp_mat}) and rank (\code{rank_mat}) matrices, and, if \code{top_n} is given,
#' the \code{feature_consistency} table.
#' @seealso ranger rf_clf.by_datasets rf_reg.by_datasets
#' @examples
#'
#' df <- data.frame(rbind(t(rmultinom(14, 14*5, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(16, 16*5, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(30, 30*5, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(30, 30*5, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(30, 30*5, c(.001,.6,.42,.58,.299)))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15),
#'                                   rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15))),
#'                      f_c=factor(c(rep("C", 60), rep("D", 60))),
#'                      f_d=factor(c(rep("A", 30), rep("B", 30), rep("C", 30), rep("D", 30))),
#'                      age=c(1:60, 2:61)
#'                      )
#' reg_res<-rf_reg.by_datasets(df, metadata, s_category='f_d', c_category='age')
#' corr_datasets_by_imps(reg_res$feature_imps_list, plot=TRUE)
#' clf_res<-rf_clf.by_datasets(df, metadata, s_category='f_s', c_category='f_c')
#' corr_datasets_by_imps(clf_res, plot=TRUE, top_n=2)
#' @author Shi Huang
#' @export
corr_datasets_by_imps<-function(feature_imps_list, ranked=TRUE, plot=FALSE, top_n=NULL, min_datasets=NULL){
  cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], ...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  if(inherits(feature_imps_list, c("rf_clf.by_datasets", "rf_reg.by_datasets", "rf_clf.comps"))){
    datasets <- feature_imps_list$datasets
    feature_imps_list <- feature_imps_list$feature_imps_list
    if(is.null(names(feature_imps_list))) names(feature_imps_list) <- datasets
  }
  feature_imps_list <- lapply(feature_imps_list, function(x){
    if(is.data.frame(x)) return(setNames(x[, "rf_imps"], as.character(x[, "feature"])))
    x
  })
  if(length(feature_imps_list) < 2) stop("At least two datasets are required.")
  if(is.null(names(feature_imps_list))) names(feature_imps_list) <- paste0("dataset", seq_along(feature_imps_list))
  feature_names <- lapply(feature_imps_list, names)
  if(any(vapply(feature_names, is.null, logical(1)))){
    if(length(unique(lengths(feature_imps_list))) > 1) stop("Importance vectors without feature names should have the same length.")
    imp_mat <- do.call(cbind, feature_imps_list)
  }else{
    shared <- Reduce(intersect, feature_names)
    if(length(shared) < 3) stop("Less than three features are shared by all datasets.")
    imp_mat <- do.call(cbind, lapply(feature_imps_list, function(x) x[shared]))
  }
  if(is.null(rownames(imp_mat))) rownames(imp_mat) <- seq_len(nrow(imp_mat))
  rank_mat <- apply(imp_mat, 2, function(x) rank(-x, na.last="keep"))
  mat <- if(ranked) rank_mat else imp_mat
  corr_mat<-stats::cor(mat, use="pairwise.complete.obs")
  # matrix of the p-value of the correlation
  p_mat <- cor.mtest(mat)

  if(plot){
    col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
    corrplot(corr = corr_mat, method="color", col=col(200),
             type="upper", addCoef.col = "black",
             tl.col="black", tl.srt=45,
             p.mat = p_mat, sig.level = 0.05, insig = "pch",
             diag = FALSE
             )
  }
  res <-list()
  res$corr_mat <- corr_mat
  res$p_mat <- p_mat
  res$imp_mat <- imp_mat
  res$rank_mat <- rank_mat
  if(!is.null(top_n)){
    if(is.null(min_datasets)) min_datasets <- ncol(rank_mat)
    in_top <- !is.na(rank_mat) & rank_mat <= top_n
    n_top <- rowSums(in_top)
    feature_consistency <- data.frame(feature=rownames(rank_mat),
                                      n_datasets_top=n_top,
                                      mean_rank=rowMeans(rank_mat, na.rm=TRUE),
                                      marker_class=ifelse(n_top >= min_datasets, "universal",
                                                          ifelse(n_top > 0, "dataset_specific", "non_marker")),
                                      top_in_datasets=apply(in_top, 1, function(v) paste(colnames(rank_mat)[v], collapse="|")),
                                      row.names=NULL)
    res$feature_consistency <- feature_consistency[order(-feature_consistency$n_datasets_top, feature_consistency$mean_rank), ]
  }
  res
}

#' @title plot_clf_res_list
#' @description Plot a summary of the performance of a list of classification models
#' output from \code{rf_clf.by_datasets}.
#' @param clf_res_list A list, the output of the function \code{rf_clf.by_datasets}.
#' @param q_cutoff A number indicating the cutoff of q values after fdr correction.
#' @param p.adj.method A string indicating the p-value correction method.
#' @param p_cutoff A number indicating the cutoff of p values.
#' @param outdir The output directory. If NULL, no file is written.
#' @param plot_height The plot height.
#' @seealso ranger rf_clf.by_datasets
#' @examples
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 15), rep("B", 15), rep("A", 15), rep("B", 15))),
#'                      f_c=factor(c(rep("C", 30), rep("D", 30))),
#'                      f_d=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15))),
#'                      age=c(1:30, 2:31)
#'                      )
#' clf_res_list<-rf_clf.by_datasets(df, metadata, s_category="f_c", c_category="f_s")
#' plot_clf_res_list(clf_res_list)
#' @author Shi Huang
#' @export
plot_clf_res_list<-function(clf_res_list, p_cutoff=0.05, p.adj.method = "bonferroni", q_cutoff=0.05, outdir=NULL, plot_height=NULL){
  datasets<-clf_res_list$datasets
  sample_size<-clf_res_list$sample_size
  rf_AUROC<-clf_res_list$rf_AUROC
  rf_AUPRC<-clf_res_list$rf_AUPRC
  feature_imps_list<-clf_res_list$feature_imps_list
  # Statistics summary of biomarkers discovery in multiple datasets
  num_sig_p.adj<-sapply(feature_imps_list, function(x) sum(x[, "non.param.test_p.adj"] < q_cutoff, na.rm=TRUE))
  num_sig_p<-sapply(feature_imps_list, function(x) sum(x[, "non.param.test_p"]< p_cutoff, na.rm=TRUE))
  num_enriched<-sapply(feature_imps_list, function(x) length(grep("enriched", x[,"Enr"])))
  num_depleted<-sapply(feature_imps_list, function(x) length(grep("depleted", x[,"Enr"])))
  summ<-data.frame(datasets=datasets, sample_size=sample_size, AUROC=rf_AUROC, AUPRC=rf_AUPRC, num_sig_p, num_sig_p.adj, num_enriched, num_depleted, row.names=NULL)
  names(summ)[1:2]<-c("Data_sets", "Sample_size")
  p_a<- ggplot(summ, aes(x=.data$Data_sets, y=.data$Sample_size)) + xlab("Data sets") + ylab("Sample size")+
    geom_bar(stat="identity", alpha=0.5, width=0.5)+
    coord_flip()+ # if want to filp coordinate
    theme_bw()
  p_b<- ggplot(summ, aes(x=.data$Data_sets, y=.data$AUROC)) + xlab("") + ylab("AUROC")+
    geom_point(shape="diamond", size=4)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    geom_hline(yintercept=0.5, linetype="dashed")+
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    scale_y_continuous(limits = c(0,1)) +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p_c<- ggplot(summ, aes(x=.data$Data_sets, y=.data$num_sig_p.adj)) + xlab("") + ylab("# of sig. features")+
    geom_point(size=2)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  summ_Enr<-reshape2::melt(summ[, c("Data_sets","num_enriched","num_depleted")], id.vars="Data_sets")
  summ_Enr[summ_Enr$variable=="num_depleted",]$value<--summ_Enr[summ_Enr$variable=="num_depleted",]$value
  enr_lim<-max(1, abs(summ_Enr$value))
  p_d<- ggplot(summ_Enr, aes(x=.data$Data_sets, y=.data$value, colour=.data$variable)) + xlab("") + ylab("Enrichment")+
    ylim(-enr_lim, enr_lim)+
    geom_hline(yintercept=0)+
    geom_point(show.legend=F, size=2)+
    geom_bar(show.legend=F, stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p1<-gridExtra::arrangeGrob(p_a, p_b, p_c, p_d, ncol = 4, nrow = 1, widths = c(4, 2, 2, 2))
  if(!is.null(outdir)){
    ggsave(filename=paste(outdir,"Datasets_AUROC.ggplot.pdf",sep=""),p1,
           width=9, height=ifelse(is.null(plot_height), 3+nrow(summ)*0.2, plot_height))
  }
  # boxplot indicating p and p.adj values of sig. features
  feature_res<-plyr::ldply(unname(feature_imps_list))
  feature_res_m<-reshape2::melt(feature_res[, c("feature","dataset", "Enr", "AUROC", "AUPRC","rf_imps", "non.param.test_p", "non.param.test_p.adj")],
                                id.vars=c("feature", "dataset", "Enr"))
  # customised colors for Enr's 3 factors
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  my3cols<-c("grey60", rev(gg_color_hue(2)))
  p2<-ggplot(feature_res_m, aes(x=.data$dataset, y=.data$value)) +
    geom_boxplot(outlier.shape = NA)+
    facet_grid(.~variable, scales="free") + scale_color_manual(values = my3cols) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    geom_jitter(aes(color=.data$Enr), position=position_jitter() ,size=1, alpha=0.4)+ #jitter
    theme(axis.line = element_line(color="black"),
          strip.background = element_rect(colour = "white"),
          panel.border = element_blank())
  if(!is.null(outdir)){
    ggsave(filename=paste(outdir,"Datasets_feature_res.ggplot.pdf",sep=""),plot=p2, width=10, height=4)
  }

  result<-list()
  result$summ<-summ
  result$summ_plot<-p1
  result$feature_res<-feature_res
  result$feature_res_plot<-p2
  return(result)
}

#' @title plot_reg_res_list
#' @description Plot a summary of the performance of a list of regression models
#' output from \code{rf_reg.by_datasets}.
#' @param reg_res_list A list object, the output of the function \code{rf_reg.by_datasets}.
#' @param plot_height The plot height.
#' @param outdir The output directory. If NULL, no file is written.
#' @seealso ranger rf_reg.by_datasets
#' @examples
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 15), rep("B", 15), rep("A", 15), rep("B", 15))),
#'                      f_c=factor(c(rep("C", 30), rep("D", 30))),
#'                      f_d=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15))),
#'                      age=c(1:30, 2:31)
#'                      )
#' reg_res_list<-rf_reg.by_datasets(df, metadata, s_category="f_c", c_category="age")
#' plot_reg_res_list(reg_res_list)
#' @author Shi Huang
#' @export
plot_reg_res_list<-function(reg_res_list, outdir=NULL, plot_height=NULL){
  datasets<-reg_res_list$datasets
  sample_size<-reg_res_list$sample_size
  summ<-data.frame(datasets=datasets, sample_size=sample_size,
                   MSE=reg_res_list$rf_MSE, RMSE=reg_res_list$rf_RMSE, MAE=reg_res_list$rf_MAE,
                   MAPE=reg_res_list$rf_MAPE, Spearman_rho=reg_res_list$rf_Spearman_rho,
                   R_squared=reg_res_list$rf_R_squared, Adj_R_squared=reg_res_list$rf_Adj_R_squared, row.names=NULL)
  names(summ)[1:2]<-c("Data_sets", "Sample_size")
  p_a<- ggplot(summ, aes(x=.data$Data_sets, y=.data$Sample_size)) + xlab("Data sets") + ylab("Sample size")+
    geom_bar(stat="identity", alpha=0.5, width=0.5)+
    coord_flip()+ # if want to filp coordinate
    theme_bw()
  p_b<- ggplot(summ, aes(x=.data$Data_sets, y=.data$RMSE)) + xlab("") + ylab("RMSE")+
    geom_point(size=2)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p_c<- ggplot(summ, aes(x=.data$Data_sets, y=.data$MAE)) + xlab("") + ylab("MAE")+
    geom_point(size=2)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p_d<- ggplot(summ, aes(x=.data$Data_sets, y=.data$R_squared)) + xlab("") + ylab("R_squared")+
    geom_point(size=2)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  summ_plot<-arrangeGrob(p_a, p_b, p_c, p_d, ncol = 4, nrow = 1, widths = c(3, 2, 2, 2))
  if(!is.null(outdir)){
    ggsave(filename=paste(outdir,"Datasets_perfs.ggplot.pdf",sep=""), summ_plot,
           width=9, height=ifelse(is.null(plot_height), 3+nrow(summ)*0.2, plot_height))
  }
  feature_res<-do.call(cbind, reg_res_list$feature_imps_list)
  result<-list()
  result$summ<-summ
  result$summ_plot<-summ_plot
  result$feature_res<-feature_res
  return(result)
}


#' @title rf_clf.by_datasets.summ
#' @description A integrated pipeline for \code{rf_clf.by_datasets}. It runs standard random forests with oob estimation for classification of
#' c_category in each the sub-datasets splited by the s_category. The output includes a summary of rf models in the sub datasets
#' and all important statistics for each of features, and optionally an R Markdown report.
#' @param df Training data: a data.frame, or a named list of data.frames.
#' @param metadata A metadata with at least two categorical variables.
#' @param s_category A string indicates the category in the sample metadata: a ‘factor’ defines the sample grouping for data spliting.
#' @param c_category A indicates the category in the sample metadata: a 'factor' used as sample label for rf classification in each of splited datasets.
#' @param positive_class A string indicates one class in the 'c_category' column of metadata.
#' @param ntree The number of trees.
#' @param nfolds The number of folds in the cross validation.
#' @param p.adj.method The p-value correction method, default is "bonferroni".
#' @param verbose Show computation status and estimated runtime.
#' @param rf_imp_pvalues If compute both importance score and pvalue for each feature.
#' @param p_cutoff The cutoff of p values for features, the default value is 0.05.
#' @param q_cutoff The cutoff of q values for features, the default value is 0.05.
#' @param outdir The output directory. If NULL, no file is written.
#' @param clr_transform A boolean value indicating if the clr-transformation applied.
#' @param cv_type The cross-validation strategy: "stratified", "group" or "logo" (see \code{rf_clf.by_datasets}).
#' @param g_category The grouping column in metadata for group k-fold or leave-one-group-out CV.
#' @param n_cores The number of cores used for parallel computation.
#' @param report_format An optional R Markdown output format ("html_document", "pdf_document" or "word_document").
#' If provided, a report is rendered by \code{crossRanger_report} into \code{outdir}.
#' @return A list includes a summary of rf models in the sub datasets,
#' all important statistics for each of features, and plots.
#' @seealso ranger rf_clf.by_datasets crossRanger_report
#' @examples
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 15), rep("B", 15), rep("A", 15), rep("B", 15))),
#'                      f_c=factor(c(rep("C", 30), rep("D", 30))))
#' res <- rf_clf.by_datasets.summ(df, metadata, s_category="f_c", c_category="f_s", nfolds=3)
#' res$perf_summ
#' @export
rf_clf.by_datasets.summ<-function(df, metadata, s_category, c_category, positive_class=NA,
                                  rf_imp_pvalues=FALSE, nfolds=3, verbose=FALSE, ntree=500,
                                  p_cutoff=0.05, p.adj.method = "bonferroni", q_cutoff=0.05,
                                  outdir=NULL, clr_transform=TRUE, cv_type="stratified", g_category=NULL,
                                  n_cores=NULL, report_format=NULL){
  res_list<-rf_clf.by_datasets(df, metadata, s_category=s_category, c_category=c_category, positive_class=positive_class,
                               rf_imp_pvalues=rf_imp_pvalues, clr_transform=clr_transform, nfolds=nfolds, verbose=verbose,
                               ntree=ntree, p.adj.method=p.adj.method, q_cutoff=q_cutoff,
                               cv_type=cv_type, g_category=g_category, n_cores=n_cores)
  plot_res_list<-plot_clf_res_list(res_list, p_cutoff=p_cutoff,
                                   p.adj.method = p.adj.method, q_cutoff=q_cutoff,
                                   outdir=outdir)
  result<-list()
  result$rf_models<-res_list$rf_model_list
  result$perf_summ<-res_list$perf_summ
  result$summ<-plot_res_list$summ
  result$summ_plot<-plot_res_list$summ_plot
  result$feature_res<-plot_res_list$feature_res
  result$feature_res_plot<-plot_res_list$feature_res_plot
  result$res_list<-res_list
  class(result)<-"rf_clf.by_datasets.summ"
  if(!is.null(report_format)){
    result$report<-crossRanger_report(Within_study_classification=result, output_format=report_format,
                                      outdir=outdir, output_file=NULL)
  }
  return(result)
}


#' @title rf_clf.comps.summ
#' @description Runs standard random forests with oob estimation for classification of
#' one level VS all other levels of one category in the datasets.
#' The output includes a summary of rf models in the sub datasets,
#' all important statistics for each of features, plots, and optionally an R Markdown report.
#' @param df Training data: a data.frame.
#' @param f A factor in the metadata with at least two levels (groups).
#' @param comp_group A string indicates the group in the f
#' @param clr_transform A boolean value indicating if the clr-transformation applied.
#' @param ntree The number of trees.
#' @param nfolds The number of folds in the cross validation. If 1, out-of-bag estimation is used.
#' @param p.adj.method The p-value correction method, default is "bonferroni".
#' @param q_cutoff The cutoff of q values for features, the default value is 0.05.
#' @param p_cutoff The cutoff of p values for features, the default value is 0.05.
#' @param verbose A boolean value indicates if showing computation status and estimated runtime.
#' @param outdir The output directory. If NULL, no file is written.
#' @param n_cores The number of cores used for parallel computation.
#' @param report_format An optional R Markdown output format ("html_document", "pdf_document" or "word_document").
#' If provided, a report is rendered by \code{crossRanger_report} into \code{outdir}.
#' @return A list includes a summary of rf models in the sub datasets,
#' all important statistics for each of features, and plots.
#' @seealso ranger
#' @examples
#' df <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' f=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' comp_group="A"
#' rf_clf.comps.summ(df, f, comp_group, verbose=FALSE, ntree=500, p_cutoff=0.05,
#'                   p.adj.method = "bonferroni", q_cutoff=0.05, outdir=NULL)
#'@export
rf_clf.comps.summ<-function(df, f, comp_group, clr_transform=TRUE, nfolds=1, verbose=FALSE, ntree=5000,
                           p_cutoff=0.05, p.adj.method = "bonferroni", q_cutoff=0.05, outdir=NULL,
                           n_cores=NULL, report_format=NULL){
  clf_comps_list<-rf_clf.comps(df, f, comp_group, verbose=verbose, ntree=ntree, clr_transform=clr_transform,
                               p.adj.method=p.adj.method, q_cutoff=q_cutoff, nfolds=nfolds, n_cores=n_cores)
  plot_res_list<-plot_clf_res_list(clf_res_list=clf_comps_list, p_cutoff=p_cutoff,
                                   p.adj.method = p.adj.method, q_cutoff=q_cutoff, outdir=outdir)
  result<-list()
  result$rf_models<-clf_comps_list$rf_model_list
  result$perf_summ<-clf_comps_list$perf_summ
  result$summ<-plot_res_list$summ
  result$summ_plot<-plot_res_list$summ_plot
  result$feature_res<-plot_res_list$feature_res
  result$feature_res_plot<-plot_res_list$feature_res_plot
  result$res_list<-clf_comps_list
  class(result)<-"rf_clf.comps.summ"
  if(!is.null(report_format)){
    result$report<-crossRanger_report(Pairwise_group_comparisons=result, output_format=report_format,
                                      outdir=outdir, output_file=NULL)
  }
  return(result)
}


#' @title plot_cross_appl
#' @description The heatmap indicating ML performance (e.g., AUROC or MAE) in the self-validation and cross-applications.
#' The datasets are shown in their original order (e.g., the order of timepoints), which reveals temporal transition patterns.
#' @param cross_appl_res The output object from \code{rf_clf.cross_appl}, \code{rf_reg.cross_appl},
#' \code{rf_clf.lodo} or \code{rf_reg.lodo}.
#' @param metric The performance metric applied. By default, "AUROC" for classification and "MAE" for regression.
#' @param plot_height The height (inches) of heatmap.
#' @param plot_width The width (inches) of heatmap
#' @param outdir The outputh directory, default is NULL.
#' @return A heat map showing ML performance in self-validation and cross-applications.
#' @examples
#' df <- data.frame(rbind(t(rmultinom(14, 14*5, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(16, 16*5, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(30, 30*5, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(30, 30*5, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(30, 30*5, c(.001,.6,.42,.58,.299)))))
#' metadata<-data.frame(f_s=factor(rep(c("A", "B"), 60)),
#'                      f_s1=factor(c(rep(TRUE, 60), rep(FALSE, 60))),
#'                      f_c=factor(c(rep("C", 30), rep("H", 30), rep("D", 30), rep("P", 30))),
#'                      age=c(1:60, 2:61)
#'                      )
#'
#' table(metadata[, c('f_s', 'f_c')])
#' clf_res<-rf_clf.by_datasets(df, metadata, nfolds=5, s_category='f_c', c_category='f_s')
#' clf_cross_appl_res <- rf_clf.cross_appl(clf_res$rf_model_list,
#'                                         x_list=clf_res$x_list,
#'                                         y_list=clf_res$y_list)
#' plot_cross_appl(clf_cross_appl_res)
#'
#' reg_res<-rf_reg.by_datasets(df, metadata, nfolds=5, s_category='f_c', c_category='age')
#' reg_cross_appl_res <- rf_reg.cross_appl(reg_res,
#'                                         x_list=reg_res$x_list,
#'                                         y_list=reg_res$y_list)
#' plot_cross_appl(reg_cross_appl_res, metric="MAE")
#'@export
plot_cross_appl<-function(cross_appl_res, metric=NULL, outdir=NULL, plot_width=8, plot_height=7){
  clf_classes <- c("rf_clf.cross_appl", "rf_clf.lodo")
  lodo_classes <- c("rf_clf.lodo", "rf_reg.lodo")
  if(!inherits(cross_appl_res, c(clf_classes, "rf_reg.cross_appl", "rf_reg.lodo")))
    stop("The class of cross_appl_res should be 'rf_clf.cross_appl', 'rf_reg.cross_appl', 'rf_clf.lodo' or 'rf_reg.lodo'")
  is_clf <- inherits(cross_appl_res, clf_classes)
  if(is.null(metric)) metric <- if(is_clf) "AUROC" else "MAE"
  perf_summ<-cross_appl_res$perf_summ
  if(!metric %in% colnames(perf_summ))
    stop("metric should be one of: ", paste(setdiff(colnames(perf_summ), c("Train_data", "Test_data", "Validation_type")), collapse=", "))
  perf_summ<-perf_summ[!is.na(perf_summ$Train_data), ]
  if(inherits(cross_appl_res, lodo_classes)) perf_summ$Train_data <- "LODO model"
  perf_summ$Train_data<-factor(perf_summ$Train_data, levels=unique(perf_summ$Train_data))
  perf_summ$Test_data<-factor(perf_summ$Test_data, levels=unique(perf_summ$Test_data))
  perf_summ[[metric]]<-as.numeric(perf_summ[[metric]])
  p<-ggplot(perf_summ, aes(x=.data$Test_data, y=.data$Train_data)) +
    xlab("Test data")+ylab("Train data")+
    geom_tile(aes(fill = .data[[metric]], color = .data$Validation_type), width=0.9, height=0.9, linewidth=1) +
    scale_color_manual(values=c("white","grey80"))+
    geom_text(aes(label = round(.data[[metric]], 2)), color = "white") +
    scale_fill_viridis()+
    labs(fill=metric, color="Validation type") +
    theme_bw() + theme_classic() +
    theme(axis.line = element_blank(), axis.text.x = element_text(angle = 90),
          axis.ticks = element_blank())
  if(!is.null(outdir)){
    write.table(cross_appl_res$perf_summ, file=paste(outdir, if(is_clf) "crossRF_clf_perf_summ.xls" else "crossRF_reg_perf_summ.xls", sep=""),
                quote=FALSE, sep="\t", row.names = F)
    ggsave(filename=paste(outdir, metric, "_cross_appl_matrix.heatmap.pdf",sep=""),
           plot=p, width=plot_width, height=plot_height)
  }
  p
}

#' @title id_non_spcf_markers
#' @description Non-specific features across datasets (at least present in two of datasets): Last update: 20190130
#' @param feature_res the inheritant output from the function of plot.res_list, rf_clf.comps.summ or rf_clf.by_dataset.summ
#' @param positive_class A string indicates one class in the 'c_category' column of metadata.
#' @param other_class A string indicates the other class in the factor, such as 'health'.
#' @param p.adj.method The p-value correction method, default is "bonferroni".
#' @param outdir The outputh directory. If NULL, no file is written.
#' @return A list including the fraction and number of all-shared, non-specific and specific markers in each dataset,
#' the feature table annotated with the marker specificity, and plots.
#' @seealso ranger
#' @examples
#' df <- data.frame(rbind(t(rmultinom(15, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' f=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' comp_group="A"
#' res <-rf_clf.comps.summ(df, f, comp_group, verbose=FALSE, ntree=500,
#'                         p_cutoff=0.05, p.adj.method = "bonferroni",
#'                         q_cutoff=0.05, outdir=NULL)
#' feature_res<-res$feature_res
#' id_non_spcf_markers(feature_res, positive_class="disease",
#'                     other_class="health", p.adj.method= "BH", outdir=NULL)
#' @author Shi Huang
#' @export
id_non_spcf_markers <- function(feature_res, positive_class="disease",
                                other_class="health", p.adj.method = "BH", outdir=NULL){
  feature_res<-feature_res[order(feature_res$dataset, feature_res$feature), ]
  if(all(feature_res$Enr=='Neutral', na.rm=TRUE)) stop("No feature significantly changed in any sub-datasets!")
  features<-as.character(feature_res$feature)
  enr<-as.character(feature_res$Enr)
  n_depleted<-tapply(grepl("depleted$", enr), features, sum)
  n_enriched<-tapply(grepl("enriched$", enr), features, sum)
  n_total<-tapply(!is.na(enr), features, sum)
  shared_global<-ifelse(n_depleted==n_total | n_enriched==n_total, "all_shared", "unshared")
  non_spcf_global<-ifelse(n_depleted==n_total | n_enriched==n_total, "all_shared",
                          ifelse(n_depleted>=2 | n_enriched>=2, "non_spcf",
                                 ifelse(n_enriched==0 & n_depleted==0, "non_marker", "spcf")))
  # customised colors for Enr's 3 factors
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  #---- statistics summary of non_spcf markers aross all datasets
  non_spcf_disease_global<-ifelse(n_enriched>=2 & n_depleted==0, paste("non_spcf_", positive_class, sep=""),
                                  ifelse(n_enriched==0 & n_depleted>=2, paste("non_spcf_", other_class, sep=""),
                                         ifelse(n_enriched>=1 & n_depleted>=1, "non_spcf_mixed",
                                                ifelse(n_enriched==0 & n_depleted==0, "non_marker", "spcf"))))
  #---- calculate the overlap fraction of non-specific markers in each of datasets
  tmp<-data.frame(feature_res, shared_global=shared_global[features],
                  non_spcf_global=non_spcf_global[features],
                  non_spcf_disease_global=non_spcf_disease_global[features],
                  non_spcf_local=non_spcf_global[features],
                  non_spcf_disease_local=non_spcf_disease_global[features],
                  stringsAsFactors =F)
  tmp[which(tmp[, "Enr"]=="Neutral"),  c("non_spcf_local", "non_spcf_disease_local")]<-"non_marker"
  count_spcf<-do.call(rbind, tapply(tmp$non_spcf_local, tmp$dataset, function(x)
    c(sum(x=="all_shared"),
      sum(x=="non_spcf"),
      sum(x=="spcf"),
      sum(x!="all_shared" & x!="non_spcf" & x!="spcf")
    )
  )
  )
  zero_count <- tapply(tmp$mean_all, tmp$dataset, function(x) sum(x==0))
  count_spcf[, 4]<-count_spcf[, 4]-zero_count
  fac_spcf<-sweep(count_spcf, 1, rowSums(count_spcf), "/")
  colnames(fac_spcf)<-colnames(count_spcf)<-c("All_shared_markers","Non_specific_markers",
                                              "Specific_markers", "Others")
  # remove the columns with total zero values
  fac_spcf<-fac_spcf[, which(!apply(fac_spcf, 2, function(x) all(x==0))), drop=FALSE]
  count_spcf<-count_spcf[, which(!apply(count_spcf, 2, function(x) all(x==0))), drop=FALSE]
  fac_spcf<-data.frame(dataset=rownames(fac_spcf), fac_spcf)
  if(!is.null(outdir)){
    # output the fac_spcf table
    .write_tsv(fac_spcf, paste(outdir,"Markers_fraction_overlap_with_non_specific_",p.adj.method,".txt",sep=""), row.names=TRUE)
    # output the count_spcf table
    .write_tsv(count_spcf, paste(outdir,"Markers_number_overlap_with_non_specific_",p.adj.method,".txt",sep=""), row.names=TRUE)
  }
  # ggplot the fraction overlap of non-specific markers in all datasets
  fac_spcf_m<-reshape2::melt(fac_spcf, id.vars="dataset")
  fac_spcf_m$variable<-factor(fac_spcf_m$variable, levels = rev(levels(fac_spcf_m$variable)))
  mycolors<-c("grey80", gg_color_hue(nlevels(fac_spcf_m$variable)-1))
  p_summ<-ggplot(fac_spcf_m, aes(x=.data$dataset, y=.data$value, fill=.data$variable)) +
    geom_bar(stat="identity")+
    scale_fill_manual(values=mycolors)+
    ylab("Fraction overlap with non-specific markers")+
    coord_flip()+
    theme_bw()
  #---- summary of the abundance and occurence rate of
  #---- non-specific disease, non-specific health, non-specific mixed and non markers
  #---- across all patients in all datasets
  feature_res_spcf<-tmp
  p_spcf_abd<-ggplot(feature_res_spcf, aes(x=.data$non_spcf_disease_global, y=log10(.data$mean_all)))  +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position=position_jitter(width = 0.2) ,size=1, alpha=0.4) + #jitter
    xlab("")+ ylab("log10(mean abundance)")+
    coord_flip()+
    theme_bw()
  p_spcf_OccRate<-ggplot(feature_res_spcf, aes(x=.data$non_spcf_disease_global, y=.data$OccRate_all))  +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position=position_jitter(width = 0.2) ,size=1, alpha=0.4) + #jitter
    xlab("")+ ylab("Ubiquity")+
    coord_flip()+
    theme_bw()
  if(!is.null(outdir)){
    ggsave(filename=paste(outdir,"Markers_fraction_overlap_with_non_specific_",p.adj.method,".pdf",sep=""),
           plot=p_summ, width=6, height=4)
    ggsave(filename=paste(outdir,"Markers_specific_VS_mean_",p.adj.method,".pdf",sep=""),
           plot=p_spcf_abd, width=5, height=3)
    ggsave(filename=paste(outdir,"Markers_specific_VS_ubiquity_",p.adj.method, ".pdf",sep=""),
           plot=p_spcf_OccRate, width=5, height=3)
  }

  result<-list()
  result$fac_spcf<-fac_spcf
  result$count_spcf<-count_spcf
  result$feature_res_spcf<-feature_res_spcf
  result$plot_spcf_summ<-p_summ
  result$plot_spcf_abd<-p_spcf_abd
  result$plot_spcf_OccRate<-p_spcf_OccRate
  return(result)
}

#' @title plot_logfc_heatmap
#' @description Heatmap of the log2 fold changes of features across datasets (e.g., cohorts, body sites or timepoints),
#' with significantly enriched or depleted features marked by asterisks. It facilitates the identification of
#' generalizable versus study-specific microbial signatures.
#' @param feature_res A data.frame of feature statistics stacked over datasets with columns "feature", "dataset",
#' "Enr" and a log2 fold change column (e.g., the \code{feature_res} of \code{plot_clf_res_list},
#' \code{rf_clf.by_datasets.summ} or \code{rf_clf.comps.summ}), or an \code{rf_clf.by_datasets}/\code{rf_clf.comps} object.
#' @param logfc_col The column of log2 fold changes: "mean_logfc", "median_logfc" or "generalized_logfc".
#' @param features An optional vector of features to show. By default, the features significant in
#' at least \code{min_sig_datasets} datasets.
#' @param min_sig_datasets The minimum number of datasets in which a feature is significant to be shown.
#' @param feature_md An optional data.frame of feature annotations (e.g., taxonomy) used to label the features.
#' @param feature_id_col The feature ID column in \code{feature_md}.
#' @param label_col The column in \code{feature_md} used as feature labels.
#' @param outdir The output directory. If NULL, no file is written.
#' @param plot_width The width (inches) of heatmap.
#' @param plot_height The height (inches) of heatmap. By default, it scales with the number of features.
#' @return A list including the matrix of log2 fold changes (features by datasets) and the heatmap.
#' @examples
#' df <- data.frame(rbind(t(rmultinom(30, 300, c(.21,.6,.12,.38,.099))),
#'                        t(rmultinom(30, 300, c(.001,.6,.42,.58,.299))),
#'                        t(rmultinom(30, 300, c(.21,.6,.12,.38,.099))),
#'                        t(rmultinom(30, 300, c(.001,.6,.42,.58,.299)))))
#' metadata <- data.frame(study=factor(rep(c("S1", "S2"), each=60)),
#'                        disease=factor(rep(rep(c("Case", "Control"), each=30), 2)))
#' res <- rf_clf.by_datasets(df, metadata, s_category="study", c_category="disease",
#'                           positive_class="Case", n_cores=1)
#' plot_logfc_heatmap(res)$plot
#' @export
plot_logfc_heatmap <- function(feature_res, logfc_col="mean_logfc", features=NULL, min_sig_datasets=1,
                               feature_md=NULL, feature_id_col="Feature_ID", label_col=NULL,
                               outdir=NULL, plot_width=6, plot_height=NULL){
  if(inherits(feature_res, c("rf_clf.by_datasets", "rf_clf.comps"))){
    feature_res <- do.call(rbind, unname(feature_res$feature_imps_list))
  }
  required_cols <- c("feature", "dataset", "Enr", logfc_col)
  missing_cols <- setdiff(required_cols, colnames(feature_res))
  if(length(missing_cols) > 0) stop("Column(s) not found in feature_res: ", paste(missing_cols, collapse=", "))
  plot_df <- data.frame(feature=as.character(feature_res$feature),
                        dataset=as.character(feature_res$dataset),
                        logFC=feature_res[[logfc_col]],
                        Enr=as.character(feature_res$Enr),
                        stringsAsFactors=FALSE)
  plot_df$significant <- !is.na(plot_df$Enr) & plot_df$Enr!="Neutral"
  if(is.null(features)){
    n_sig <- tapply(plot_df$significant, plot_df$feature, sum)
    features <- names(n_sig)[n_sig >= min_sig_datasets]
  }
  if(length(features)==0) stop("No feature is significant in at least ", min_sig_datasets, " dataset(s).")
  plot_df <- plot_df[plot_df$feature %in% features, ]
  logfc_mat <- tapply(plot_df$logFC, list(plot_df$feature, plot_df$dataset), mean)
  logfc_mat <- logfc_mat[, unique(plot_df$dataset), drop=FALSE]
  # cluster features by their fold change profiles
  feature_order <- rownames(logfc_mat)
  if(nrow(logfc_mat) > 2){
    mat <- logfc_mat
    mat[is.na(mat)] <- 0
    feature_order <- rownames(mat)[hclust(dist(mat))$order]
  }
  plot_df$feature <- factor(plot_df$feature, levels=feature_order)
  plot_df$dataset <- factor(plot_df$dataset, levels=colnames(logfc_mat))
  plot_df$label <- ifelse(plot_df$significant, "*", "")
  feature_labels <- feature_order
  if(!is.null(feature_md) && !is.null(label_col)){
    matched <- match(feature_order, as.character(feature_md[, feature_id_col]))
    feature_labels <- ifelse(is.na(matched), feature_order, as.character(feature_md[matched, label_col]))
  }
  lim <- max(abs(plot_df$logFC), na.rm=TRUE)
  p <- ggplot(plot_df, aes(x=.data$dataset, y=.data$feature, fill=.data$logFC)) +
    geom_tile(color="white") +
    geom_text(aes(label=.data$label), vjust=0.75) +
    scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0, limits=c(-lim, lim),
                         name=expression(log[2]~fold~change)) +
    scale_y_discrete(labels=setNames(feature_labels, feature_order)) +
    xlab("Data sets") + ylab("") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank())
  if(!is.null(outdir)){
    height <- if(is.null(plot_height)) 2 + 0.2*length(feature_order) else plot_height
    ggsave(filename=paste(outdir, logfc_col, ".heatmap.pdf", sep=""), plot=p, width=plot_width, height=height)
  }
  res <- list()
  res$logfc_mat <- logfc_mat[feature_order, , drop=FALSE]
  res$plot <- p
  res
}

#' @title id_robust_markers
#' @description Intersect the machine-learning-derived feature rankings with the statistically significant features
#' from univariate tests to obtain a robust biomarker set with both predictive power and biological interpretability.
#' @param rf_features The features selected by random forests: the result of \code{rf_clf.rfe}/\code{rf_reg.rfe}
#' (its selected features), a model from \code{rf.out.of.bag}/\code{rf.cross.validation} (its \code{top_n} features),
#' or a character vector of feature IDs.
#' @param bg_test The output data.frame of \code{BetweenGroup.test} (rownames are feature IDs).
#' @param top_n The number of top-ranked features used when \code{rf_features} is a model.
#' @param q_cutoff An optional cutoff of adjusted p values of the non-parametric test. By default, the enrichment
#' (\code{Enr}) of \code{BetweenGroup.test} defines the significant features.
#' @return A data.frame of all features with their RF importance rank (if available), the adjusted p value,
#' log2 fold change, enrichment, and whether they are robust markers (selected by RF and statistically significant).
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(30, 300, c(.21,.6,.12,.38,.099))),
#'                       t(rmultinom(30, 300, c(.001,.6,.42,.58,.299)))))
#' y <- factor(rep(c("Case", "Control"), each=30))
#' rf_model <- rf.cross.validation(x, y, nfolds=5)
#' bg_test <- BetweenGroup.test(x, y, clr_transform=TRUE, p.adj.method="fdr", q_cutoff=0.05)
#' id_robust_markers(rf_model, bg_test, top_n=3)
#' @export
id_robust_markers <- function(rf_features, bg_test, top_n=20, q_cutoff=NULL){
  rf_rank <- NULL
  if(inherits(rf_features, c("rf_clf.rfe", "rf_reg.rfe"))){
    selected <- rf_features$selected_features
    rf_rank <- rank(-rf_features$importances, na.last="keep", ties.method="first")
  }else if(inherits(rf_features, c("rf.out.of.bag", "rf.cross.validation"))){
    rf_rank <- .feature_rank(rf_features)
    selected <- names(rf_rank)[!is.na(rf_rank) & rf_rank <= top_n]
  }else if(is.character(rf_features)){
    selected <- rf_features
  }else{
    stop("rf_features should be an rf_clf.rfe/rf_reg.rfe result, an rf.out.of.bag/rf.cross.validation model, or a character vector.")
  }
  if(is.null(q_cutoff)){
    significant <- !is.na(bg_test$Enr) & bg_test$Enr!="Neutral"
  }else{
    significant <- !is.na(bg_test$non.param.test_p.adj) & bg_test$non.param.test_p.adj < q_cutoff
  }
  features <- rownames(bg_test)
  out <- data.frame(feature=features,
                    rf_rank=if(is.null(rf_rank)) NA else as.numeric(rf_rank[features]),
                    rf_selected=features %in% selected,
                    non.param.test_p.adj=bg_test$non.param.test_p.adj,
                    mean_logfc=bg_test$mean_logfc,
                    Enr=bg_test$Enr,
                    significant=significant,
                    robust_marker=features %in% selected & significant,
                    row.names=NULL)
  out[order(!out$robust_marker, out$rf_rank, out$non.param.test_p.adj), ]
}
