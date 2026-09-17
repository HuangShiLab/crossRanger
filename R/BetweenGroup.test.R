#' @import foreach
#' @importFrom parallel detectCores
#' @importFrom rlang .data
#' @importFrom stats sd var t.test wilcox.test bartlett.test oneway.test kruskal.test p.adjust


#' @title BetweenGroup.test
#' @description It runs standard univarate statistical tests (such as t.test, wilcox.test,
#' oneway.test, kruskal.test)  for all variates (columns) in the data.frame/data.matrix.
#' @param x A data.martrix or data.frame including multiple numeric vectors.
#' @param y A factor with two or more levels.
#' @param clr_transform A logical value indicates if clr transformation before statistical analysis of compositional microbiome data.
#' @param p.adj.method A string indicating the p-value correction method.
#' @param positive_class A string indicating the specified class in the factor y.
#' @param q_cutoff A number indicating the cutoff of q values after fdr correction.
#' @param paired A logical indicating if paired between-group comparison is desired.
#' For paired tests, the samples of the two groups should be in the same order of pairs.
#' @param pseudocount The value added before the log transformation of the fold changes. By default
#' (\code{NULL}) it is derived per feature from the data: half of the smallest positive value of that
#' feature, and only if the feature contains a zero. Pass a fixed value (e.g. \code{1e-6} for relative
#' abundances) when the fold changes of several datasets or strata are to be compared, so that the
#' shrinkage of features containing zeros does not differ between them.
#' @return A data.frame including descriptive statistics of all samples and each group,
#' the log2 fold changes (\code{mean_logfc}, \code{median_logfc}, \code{generalized_logfc}),
#' p values and adjusted p values of parametric and non-parametric tests, and the enrichment (\code{Enr}) of each feature.
#' @seealso clr, t.test, wilcox.test, oneway.test, kruskal.test
#' @examples
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' # clr_x<-compositions::clr(x)
#' y<-factor(c(rep("A", 30), rep("B", 30)))
#' y1<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' system.time(BetweenGroup.test(x, y, clr_transform=FALSE))
#' system.time(BetweenGroup.test(x, y, clr_transform=TRUE))
#' system.time(BetweenGroup.test(x, y1, clr_transform=TRUE))
#' x_ <- data.frame(rbind(t(rmultinom(7, 7500, rep(c(.201,.5,.02,.18,.099), 100))),
#'             t(rmultinom(8, 7500, rep(c(.201,.4,.12,.18,.099), 100))),
#'             t(rmultinom(15, 7500, rep(c(.011,.3,.22,.18,.289), 100))),
#'             t(rmultinom(15, 7500, rep(c(.091,.2,.32,.18,.209), 100))),
#'             t(rmultinom(15, 7500, rep(c(.001,.1,.42,.18,.299), 100)))))
#' y_<-factor(c(rep("A", 30), rep("B", 30)))
#' y_1<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' system.time(BetweenGroup.test(x_, y_))
#' system.time(BetweenGroup.test(x_, y_1))
#' @author Shi Huang
#' @export
BetweenGroup.test <-function(x, y, clr_transform=FALSE, p.adj.method="bonferroni", positive_class=NA, q_cutoff=0.2, paired=FALSE, pseudocount=NULL){
  # p.adjust.methods
  # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  y<-factor(y)
  if(length(y)!=nrow(x)) stop("The length of y should match the number of rows of x.")
  n_group<-nlevels(y)
  if(n_group < 2)
    stop("y should be a factor with at least two levels\n")
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  if(!positive_class %in% levels(y)) stop("The positive_class '", positive_class, "' is not a level of y.")
  if(clr_transform){ xx<-compositions::clr(x) }else{xx<-x}
  test_output<-mttest(xx, y, p.adj.method=p.adj.method, paired=paired)
  #  descriptive statistics of each feature (column) for all samples
  desc_stats_all_df<-desc_stats_all(x, y, positive_class=positive_class, clr_transform=clr_transform)
  #  descriptive statistics of each feature (column) for samples (rows) grouped by y
  desc_stats_by_group_df<-desc_stats_by_group(x, y, clr_transform=clr_transform, positive_class=positive_class,
                                              pseudocount=pseudocount)
  #-------------------------------Enrichment
  Enr_all<-Enr_by_q_cutoff(test_output, desc_stats_by_group_df, positive_class, q_cutoff=q_cutoff)
  #-------------------------------
  output1<-data.frame(desc_stats_all_df, desc_stats_by_group_df, test_output,  Enr_all)
  output1
}


#' @title log.mat
#' @description Log transformation of a data matrix with pseudo count.
#' @param x A data.martrix, data.frame or numeric vector.
#' @param base The base in the log-transformation.
#' @param pseudocount The value added to \code{x} before the log transformation. By default (\code{NULL})
#' it is chosen from the data: half of the minimum positive value, and only if \code{x} contains a zero.
#' A fixed value makes the transformation independent of the data at hand, which is needed when the results
#' of several datasets are to be compared; \code{0} disables the pseudo count.
#' @details If any value is zero, half of the minimum positive value is added to all values as a pseudo count.
#' Missing values are ignored when the pseudo count is derived and are returned as \code{NA}.
#' @examples
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' log.mat(x, base=10)
#' @rdname log.mat
#' @author Shi Huang
#' @rawNamespace export(log.mat)
#' @rawNamespace S3method(log, mat)
log.mat<-function(x, base=2, pseudocount=NULL){
  if(is.null(pseudocount)){
    pseudocount <- 0
    if(any(x == 0, na.rm=TRUE)) {
      v <- unlist(x)
      v <- v[!is.na(v) & v > 0]
      pseudocount <- if(length(v) == 0) 1e-5 else min(v)/2
    }
  }
  if(pseudocount != 0) x <- x + pseudocount
  out<-log(x, base)
  return(out)
}


#' @title mttest
#' @description Perform the univariate test for all features in the compositional microbiome data.
#' Two-sample t test and Wilcoxon rank-sum (or signed-rank if paired) test for two groups,
#' and Welch's one-way ANOVA and Kruskal-Wallis test for more than two groups.
#' @param x A data.martrix or data.frame including multiple numeric vectors.
#' @param y A factor with two or more levels.
#' @param p.adj.method A string indicating the p-value correction method.
#' @param paired A logical indicating if paired between-group comparison is desired.
#' @return A matrix of p values and adjusted p values. The p value is NA if a test cannot be computed
#' (e.g., the values are constant within both groups).
#' @examples
#' y <-factor(c(rep("A", 30), rep("B", 30)))
#' y1 <-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' mttest(x, y)
#' mttest(x, y1)
#' @author Shi Huang
#' @export
mttest<-function(x, y, p.adj.method="bonferroni", paired=FALSE){
  y<-factor(y)
  test_output<-matrix(NA, ncol=4, nrow=ncol(x))
  rownames(test_output)<-colnames(x)
  colnames(test_output)<- c("param.test_p","non.param.test_p","param.test_p.adj","non.param.test_p.adj")
  # comb function for parallelization using foreach
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  safe_p <- function(test) tryCatch(test$p.value, error=function(e) NA_real_)
  # run the features in parallel only if a backend is registered; inside the workers of another
  # parallel loop there is none, and %do% avoids both nested parallelism and its warning
  `%run%` <- if(foreach::getDoParRegistered()) foreach::`%dopar%` else foreach::`%do%`
  oper<-foreach::foreach(i=1:ncol(x), .combine='comb', .multicombine=TRUE, .init=list(c(), c())) %run% {
    xi<-as.numeric(x[, i])
    if(stats::var(xi)==0){test_out1<-test_out2<-1
    }else{
      if(nlevels(y)==2){
        g1<-xi[y==levels(y)[1]]
        g2<-xi[y==levels(y)[2]]
        test_out1<-safe_p(t.test(g1, g2, paired=paired))
        test_out2<-safe_p(wilcox.test(g1, g2, paired=paired, exact=FALSE, correct=FALSE))
      }else{
        test_out1<-safe_p(oneway.test(xi~y, var.equal=FALSE))
        test_out2<-safe_p(kruskal.test(xi~y))
      }
    }
    out<-c(test_out1, test_out2)
  }
  test_output[,1]<-unlist(oper[[1]])
  test_output[,2]<-unlist(oper[[2]])
  test_output[,3]<-p.adjust(test_output[,1], method = p.adj.method, n = ncol(x))
  test_output[,4]<-p.adjust(test_output[,2], method = p.adj.method, n = ncol(x))
  test_output
}



#' @title desc_stats_all
#' @description The descriptive statistics summary of all features of compositional microbiome data.
#' @param x A data.martrix or data.frame including multiple numeric vectors.
#' @param y A factor with two or more levels.
#' @param clr_transform A logical value indicates if clr transformation before statistical analysis of compositional microbiome data.
#' @param positive_class A string indicating the specified class in the factor y. The positive class should be specified in the case-control design.
#' @examples
#' y <-factor(c(rep("A", 30), rep("B", 30)))
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' desc_stats_all(x, y, clr_transform=FALSE)
#' desc_stats_all(x, y, clr_transform=TRUE)
#' @author Shi Huang
#' @export
desc_stats_all<-function(x, y, positive_class=NA, clr_transform=FALSE){
  y<-factor(y)
  OccRate<-function(x) sum(x!=0)/length(x)
  func_all_list<-c("mean_all", "var_all", "sd_all", "OccRate_all", "AUROC", "AUPRC")
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  desc_stats_all_df<-data.frame(t(apply(x,2,function(a)
    c(mean(a), stats::var(a), stats::sd(a), OccRate(a), get.auroc(a,y, positive_class), get.auprc(a,y, positive_class)))))
  colnames(desc_stats_all_df)<-func_all_list
  if(clr_transform){
    clr_x<-compositions::clr(x)
    func_all_list<-c("clr_mean_all", "clr_var_all", "clr_sd_all", "clr_AUROC", "clr_AUPRC")
    clr_desc_stats_all_df<-data.frame(t(apply(clr_x, 2, function(a)
      c(mean(a), stats::var(a), stats::sd(a), get.auroc(a,y, positive_class), get.auprc(a,y, positive_class)))))
    colnames(clr_desc_stats_all_df)<-func_all_list
    desc_stats_all_df<-data.frame(desc_stats_all_df, clr_desc_stats_all_df)
  }
  desc_stats_all_df
}


#' @title desc_stats_by_group_using_mean_logfc
#' @description The descriptive statistics summary of all features by grouping of compositional microbiome data.
#' @param x A data.martrix or data.frame including multiple numeric vectors.
#' @param y A factor with two or more levels.
#' @param clr_transform A logical value indicates if clr transformation before statistical analysis of compositional microbiome data.
#' @param positive_class A string indicating the specified class in the factor y. The positive class should be specified in the case-control design.
#' @examples
#' y <-factor(c(rep("A", 30), rep("B", 30)))
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' desc_stats_by_group(x, y, clr_transform=FALSE)
#' desc_stats_by_group(x, y, clr_transform=TRUE)
#' @param pseudocount The value added to the group summaries before the log transformation of the fold
#' changes. By default (\code{NULL}) it is derived per feature from the data (see \code{\link{log.mat}}).
#' Pass a fixed value when the fold changes of several datasets are to be compared, so that features
#' containing zeros are shrunk equally in all of them.
#' @author Shi Huang
#' @export
desc_stats_by_group<-function(x, y, clr_transform=FALSE, positive_class=NA, pseudocount=NULL){
  y<-factor(y)
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  OccRate<-function(x) sum(x!=0)/length(x)
  log10_median<-function(x, base=10) log.mat(stats::median(x), base=10)

  # the group summaries are log-transformed one feature at a time: with a pseudo count taken from the
  # data, transforming the whole feature-by-group matrix at once would let the smallest value anywhere
  # in the table set the pseudo count of every feature, so that the fold change of one feature would
  # depend on the other features the table happens to contain.
  logfc_by <- function(x, y, FUN, base=2){
    if(nlevels(y)>2) levels(y)[levels(y)!=positive_class] <- "Others"
    abd <- t(apply(x, 2, function(v) tapply(v, y, FUN)))
    logAbd <- abd
    for(i in seq_len(nrow(abd)))
      logAbd[i, ] <- log.mat(abd[i, ], base=base, pseudocount=pseudocount)
    logAbd[, positive_class] - logAbd[, colnames(logAbd)!=positive_class]
  }

  mean_logfc <- function(x, y, base=2) logfc_by(x, y, mean, base=base)

  median_logfc <- function(x, y, base=2) logfc_by(x, y, stats::median, base=base)

  new_quantile <- function(x) {
    r <- quantile(x, probs = seq(0.05, 0.95, 0.05), na.rm = TRUE) # NaN from log of negative (e.g., pre-transformed) values
    return (r)
  }

  generalized_logfc <- function(x, y, base=2){
    if(nlevels(y)>2) levels(y)[levels(y)!=positive_class] <- "Others"
    logQuantileAbd <- apply(x,2,function(x) tapply(log.mat(x, base = base, pseudocount = pseudocount),
                                                  y, new_quantile))
    results <- rep(0, ncol(x))
    names(results) <- colnames(x)
    for(i in 1:length(logQuantileAbd)) {
      logQuantileAbd_i <- logQuantileAbd[[i]]
      out<- logQuantileAbd_i[[positive_class]]-logQuantileAbd_i[[which(names(logQuantileAbd_i)!=positive_class)]]
      out<-mean(out)
      results[i] <- out
    }
    results
  }

  func_by_group_list<-c("mean", "sd", "median", "log10_median", "OccRate")
  tmp<-apply(x,2,function(x) tapply(x, y, function(x) c(mean(x), stats::sd(x), stats::median(x), log10_median(x), OccRate(x))));
  desc_stats_by_group_df<-t(sapply(tmp, unlist));
  colnames(desc_stats_by_group_df)<-unlist(lapply(levels(y), function(x) paste(func_by_group_list, x, sep="__")))
  desc_stats_by_group_df<-data.frame(desc_stats_by_group_df,
                                     mean_logfc=mean_logfc(x, y),
                                     median_logfc = median_logfc(x, y),
                                     generalized_logfc = generalized_logfc(x, y))
  if(clr_transform){
    clr_x<-compositions::clr(x)
    func_by_group_list<-c("clr_mean", "clr_sd", "clr_median")
    tmp<-apply(clr_x, 2, function(x) tapply(x, y, function(x) c(mean(x), stats::sd(x), stats::median(x))));
    clr_desc_stats_by_group_df<-t(sapply(tmp, unlist));
    colnames(clr_desc_stats_by_group_df)<-unlist(lapply(levels(y), function(x) paste(func_by_group_list, x, sep="__")))
    desc_stats_by_group_df<-data.frame(desc_stats_by_group_df, clr_desc_stats_by_group_df)
  }
  desc_stats_by_group_df
}

#' @title Enr_by_q_cutoff
#' @description The significance of all features by q values.
#' @param test_output The output from function \code{mttest}.
#' @param desc_stats_by_group_df The output from function \code{desc_stats_by_group}.
#' @param positive_class A string indicating the specified class in the factor y.
#' @param q_cutoff A number indicating the cutoff of q values after fdr correction.
#' @return A data.frame with the significance (\code{IfSig}), and the enrichment (\code{Enr}: "Neutral",
#' "<positive_class>_depleted" or "<positive_class>_enriched") of each feature.
#' @examples
#' y <-factor(c(rep("A", 30), rep("B", 30)))
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' test_output<-mttest(x, y, p.adj.method="bonferroni", paired=FALSE)
#' desc_stats_by_group_df<-desc_stats_by_group(x, y, clr_transform=FALSE)
#' Enr_by_q_cutoff(test_output, desc_stats_by_group_df, q_cutoff=0.05, "A")
#' @export
Enr_by_q_cutoff<-function(test_output, desc_stats_by_group_df, positive_class=NA, q_cutoff=0.05){
  enriched<-paste(positive_class, "enriched", sep="_")
  depleted<-paste(positive_class, "depleted", sep="_")
  IfSig<-factor(ifelse(test_output[, "non.param.test_p.adj"]< q_cutoff, "Sig", "NotSig"), levels=c("NotSig", "Sig"))
  Enr0<-factor(ifelse(desc_stats_by_group_df$mean_logfc>0, enriched, depleted), levels=c(depleted, enriched))
  IfSigEnr<-interaction(IfSig, Enr0)
  Enr<-factor(ifelse(IfSig=="Sig", as.character(Enr0), "Neutral"), levels=c("Neutral", depleted, enriched))
  data.frame(IfSig,IfSigEnr,Enr)
}
