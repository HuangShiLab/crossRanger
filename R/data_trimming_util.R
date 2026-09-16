#' @importFrom graphics hist

#' @title check_metadata
#' @description Summarize the completeness and the data type (numeric or categorical) of each metadata variable.
#' @param metadata A dataframe with > two columns corresponds to samples (rownames) in the biological data.
#' @param more_missing_values A optional string(s) can be added to define the missing values.
#' @param unique_rate_thres The minimum ratio of unique values to non-missing values for an integer variable to be treated as numeric rather than categorical.
#' @return A data.frame summarizing each of metadata variables.
#' @examples
#' set.seed(123)
#' a<-factor(c(rep("A", 29), NA, rep("B", 29), NA))
#' b<-factor(c(rep("A", 27), NA, "Not applicable", "Missing:not collected", rep("B", 28), NA, NA))
#' c<-factor(c(rep("A", 20), rep("B", 20), rep("C", 20)))
#' d<-c(sample(1:18), NA, "Not applicable", sample(5:44))
#' e<-c(rnorm(40), NA, "Not applicable", rnorm(18, 4))
#' e0<-c(NA, "Not applicable", rnorm(58, 4))
#' f<-rep("C", 60)
#' g<-rep(4, 60)
#' h<-c(rep("C", 59), "B")
#' metadata<-data.frame(a, b, c, d, e, e0, f, g, h)
#' check_metadata(metadata)
#' @author Shi Huang
#' @export
check_metadata <- function(metadata, more_missing_values=NULL, unique_rate_thres=0.2){
  missing_values<-c("not applicable", "Not applicable", "Missing:not collected",
                    "Not provided", "missing: not provided", "unknown",
                    "not provided", "not collected", "NA", NA, "", more_missing_values)

  check_integers <- function(vector){
    checker <- grepl("^[0-9]+$", as.character(vector), perl = T)
    total_len <- length(vector)
    n_integers <- sum(checker)
    n_non_integers <- total_len - n_integers
    if_all_integers <- all(checker)
    out <- list(if_all_integers=if_all_integers,
                n_integers=n_integers,
                n_non_integers=n_non_integers)
    return(out)
  }

  check_characters <- function(vector){
    checker <- grepl("[A-Za-z-$_]+", as.character(vector), perl = T) # [A-Za-z[:punct:]]+
    total_len <- length(vector)
    n_characters <- sum(checker)
    n_non_characters <- total_len - n_characters
    if_character_existed <- any(checker) ## if any TRUE in the checker

    out <- list(if_character_existed=if_character_existed,
                n_characters=n_characters,
                n_non_characters=n_non_characters)
    return(out)
  }

  check_numeric <- function(vector){
    checker <- suppressWarnings(as.numeric(as.character(vector)) %% 1 != 0)
    total_len <- length(vector)
    n_numeric <- sum(checker[!is.na(checker)]) # if any NA values it means a string appear among the numeric values.
    n_non_numeric <- total_len - n_numeric
    if_numeric_existed <- any(checker[!is.na(checker)])
    out <- list(if_numeric_existed=if_numeric_existed,
                n_numeric=n_numeric,
                n_non_numeric=n_non_numeric)
    return(out)
  }

  n_unique_values <- sapply(metadata, function(x) nlevels(factor(x)))
  unique_values <- sapply(metadata, function(x) paste0(levels(factor(x)), collapse="|"))
  n_missing_values <- sapply(metadata, function(x) sum(x %in% missing_values))
  n_real_values <- sapply(metadata, function(x) sum(!x %in% missing_values))
  completeness <- n_real_values/nrow(metadata)
  unique_rate <- ifelse(n_real_values==0, 0, n_unique_values/n_real_values)
  if_character_existed <- sapply(metadata, function(x) check_characters(x[!x %in% missing_values])[[1]])
  n_characters <- sapply(metadata, function(x) check_characters(x[!x %in% missing_values])[[2]])
  if_all_integers <- sapply(metadata, function(x) check_integers(x[!x %in% missing_values])[[1]])
  n_integers <- sapply(metadata, function(x) check_integers(x[!x %in% missing_values])[[2]])
  if_numeric_existed <- sapply(metadata, function(x) check_numeric(x[!x %in% missing_values])[[1]])
  n_numeric <- sapply(metadata, function(x) check_numeric(x[!x %in% missing_values])[[2]])
  all_values_identical <- sapply(metadata, function(x) length(unique(x))==1)
  numeric_var <- if_numeric_existed | (if_all_integers & unique_rate >= unique_rate_thres)
  categorical_var <- if_character_existed | (if_all_integers & unique_rate < unique_rate_thres)
  metadata_summ<-data.frame(metadata=colnames(metadata),
                            all_values_identical,
                            total_n=nrow(metadata),
                            n_missing_values,
                            n_real_values,
                            completeness,
                            n_unique_values,
                            unique_values,
                            unique_rate,
                            if_character_existed,
                            n_characters,
                            if_all_integers,
                            n_integers,
                            if_numeric_existed,
                            n_numeric,
                            numeric_var,
                            categorical_var)

  metadata_summ
}

normalize_NA_in_metadata<-function(md){
  na_strings <- c("not provided", "Not provided", "Not Provided",
                  "not applicable", "Not applicable",
                  "Missing:not collected",
                  "NA", "na", "Na",
                  "none", "None", "NONE")
  md[] <- lapply(md, function(x) {
    x[x %in% na_strings] <- NA
    if(is.factor(x)) x <- droplevels(x)
    x
  })
  md
}

discard_uninfo_columns_in_metadata<-function(md){
  noninfo_idx<-which(sapply(md, function(x) length(unique(x))==1))
  if(length(noninfo_idx) > 0) md<-md[, -noninfo_idx, drop=FALSE]
  md
}

trim_metadata <- function(md, completeness_threshold=0.5, filter_cols_by_type="numeric"){
  if(!is.element(filter_cols_by_type, c("numeric", "categorical", NA)))
    stop("Only 'numeric', 'categorical', or NA are allowed for 'filter_cols_by_type'!")
  md<-normalize_NA_in_metadata(md)
  md_summ<-check_metadata(metadata = md)
  keep <- md_summ$n_real_values!=md_summ$n_unique_values &
    md_summ$completeness>completeness_threshold &
    md_summ$n_unique_values>1
  if(is.na(filter_cols_by_type)){
    cat("No column was filtered by data type!\n")
  }else if(filter_cols_by_type=="numeric"){
    keep <- keep & md_summ$numeric_var
  }else{
    keep <- keep & md_summ$categorical_var
  }
  vars<-md_summ$metadata[keep]
  return(md[, vars, drop=FALSE])
}

chop_seq_to_x_nt<-function(df, start=1, nt=100){
  cat("The length of feature sequences all equal to ", nt, " before :", all(nchar(colnames(df))==nt), "\n")
  if(!all(nchar(colnames(df))==nt)){
    if(start==1){end <- nt}else{end <- start + nt -1 }
    colnames(df)<-substr(colnames(df), start, end)
  }
  cat("The length of feature sequences all equal to ", nt, " after :", all(nchar(colnames(df))==nt), "\n")
  df
}


filter_features_allzero<-function(data, samples=TRUE, features=TRUE){
  if(samples & features){
    result<-data[which(rowSums(data)!=0), , drop=FALSE]
    result<-result[, which(colSums(result)!=0), drop=FALSE]
  }else if(samples & !features){
    result<-data[which(rowSums(data)!=0), , drop=FALSE]
  }else if(!samples & features){
    result<-data[, which(colSums(data)!=0), drop=FALSE]
  }else{
    stop("Nothing has been done!")
  }
  result
}

filter_features_by_threshold <- function(data, threshold) {
  if (threshold < 0 || threshold > 1) {
    stop("Threshold should be a value between 0 and 1.")
  }
  min_nonzero_count <- (1-threshold) * nrow(data)
  result <- data[, which(colSums(data != 0) >= min_nonzero_count), drop=FALSE]
  return(result)
}

filter_features_by_prev <- function(data, prev=0.001){
  data<-data[, which(colSums(data!=0) > prev * nrow(data)), drop=FALSE]
  data
}

filter_features_by_abundance <- function(data, mean_abd_cutoff=0.001, plot=TRUE){
  if(plot) hist(colMeans(data))
  data<-data[, which(colMeans(data) > mean_abd_cutoff), drop=FALSE]
  data
}

remove_noisy_feature_by_count_cutoff <- function(data, count_cutoff=1000){
  if(is.data.frame(data)) data<-data.matrix(data)
  cat("The total number of zeros in the table:", sum(data==0), "\n")
  data[data<=count_cutoff]<-0
  cat("The total number of zeros in the filtered table:", sum(data==0))
  data<-data.frame(data, check.names = FALSE)
  data
}

filter_features_with_NA <- function(data, min_prev=0){
  idx <- which(colSums(is.na(data)) <= min_prev * nrow(data))
  data<-data[, idx, drop=FALSE]
  data
}

filter_samples_with_NA <- function(data, min_prev=0){
  idx <- which(rowSums(is.na(data)) <= min_prev * ncol(data))
  data<-data[idx, , drop=FALSE]
  data
}


filter_samples_by_NA_in_y <- function(data, y){
  y_k<-y[which(!is.na(y))]
  data_k<-data[which(!is.na(y)), , drop=FALSE]
  result<-list()
  result$data_k<-data_k
  result$y_k<-y_k
  result
}

#' @title filter_samples_by_groups_in_target_field_of_metadata
#' @description Keep (or remove) the samples belonging to the specified groups of a metadata variable.
#' @param data A data.frame.
#' @param metadata A data.frame including multiple metadata variables corresponds to samples in the data.
#' @param target_field A character string indicating a target variable in the metadata.
#' @param groups A character string(s) indicating one or multiple groups in the target metadata variable specified.
#' @param negate A bool value indicates if samples in the specified groups should be removed (TRUE) or kept (FALSE).
#' @param ids_col The column of metadata holding sample IDs. If NA, the rownames of metadata are used.
#' @return A list with the filtered \code{data} and \code{metadata}.
#' @examples
#' set.seed(123)
#' data <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' z<-factor(c(rep("E", 20), rep("C", 20), rep("D", 20)))
#' y<-factor(c(rep("A", 30), rep("B", 30)))
#' y0<-factor(c(rep("A", 5), rep("B", 55)))
#' metadata <- data.frame(z, y, y0)
#' filter_samples_by_groups_in_target_field_of_metadata(data, metadata, target_field="z",
#'                                                      groups=c("D", "A"), negate=TRUE)
#' @export
filter_samples_by_groups_in_target_field_of_metadata <- function(data, metadata, target_field, groups, negate=FALSE, ids_col=NA){
  if(is.na(ids_col)){
    SampleIDs<-rownames(metadata)
  }else{
    SampleIDs<-metadata[, ids_col]
  }
  target_field_levels <- levels(factor(metadata[, target_field]))
  if(!identical(rownames(data), SampleIDs)) stop("The sample IDs should be idenical in feature table and metadata!")
  if_groups_exist<-target_field_levels %in% groups
  if(all(if_groups_exist==FALSE)) stop("All specified groups do not exist in the target field!")
  if(any(if_groups_exist==FALSE)) cat("Only specified groups:",  target_field_levels[if_groups_exist], "existed in the target field!\n")
  groups<-target_field_levels[if_groups_exist]
  if(negate){
    idx<-which(!metadata[, target_field] %in% groups)
  }else{
    idx<-which(metadata[, target_field] %in% groups)
  }
  metadata_k<-metadata[idx, ]
  metadata_k[, target_field] <- factor(metadata_k[, target_field])
  data_k<-data[idx, ]
  cat("The number of kept samples (only samples related to ", groups," in ",target_field,"): ", nrow(metadata_k) ,"\n")
  result<-list()
  result$data<-data_k
  result$metadata<-metadata_k
  result
}

filter_samples_by_NA_in_target_field_of_metadata <- function(data, metadata, target_field, ids_col=NA){
  if(is.na(ids_col)){
    SampleIDs<-rownames(metadata)
  }else{
    SampleIDs<-metadata[, ids_col]
  }
  if(!identical(rownames(data), SampleIDs)) stop("The sample IDs should be idenical in feature table and metadata!")
  NAN_values<-c("not applicable", "Not applicable", "Missing:not collected",
                "Not provided", "missing: not provided",
                "not provided", "not collected", "NA", NA, "")
  idx<-which(!metadata[, target_field] %in% NAN_values)
  metadata_k<-metadata[idx, ]
  data_k<-data[idx, ]
  cat("The number of kept samples (after filtering out samples with NA values in ",target_field,"): ", nrow(metadata_k) ,"\n")
  result<-list()
  result$data<-data_k
  result$metadata<-metadata_k
  result
}

split_dm_by_metadata<-function(dm, metadata, split_factor){
  if(!is.element(split_factor, colnames(metadata)))
    stop("The split_factor you specified should be one of the column names of input metadata.")
  if(inherits(dm, "dist")) dm <- data.matrix(dm)
  f<-metadata[, split_factor]
  sub_dm_name_list<-split(1:nrow(dm), f, drop=TRUE)
  sub_dm_list<-
    lapply(sub_dm_name_list, function(x){
      list(dm[x, x], metadata[x, ])
    })

  sub_dm_list

}

filter_dm_by_NA_in_target_field_of_metadata <- function(dm, metadata, target_field, ids_col=NA){
  if(is.na(ids_col)){
    SampleIDs<-rownames(metadata)
  }else{
    SampleIDs<-metadata[, ids_col]
  }
  if(!identical(rownames(dm), SampleIDs)) stop("The sample IDs should be idenical in feature table and metadata!")
  NAN_values<-c("not applicable", "Not applicable", "Missing:not collected",
                "Not provided", "missing: not provided", "unknown", "Unknown",
                "not provided", "not collected", "NA", NA, "")
  idx<-which(!metadata[, target_field] %in% NAN_values)
  metadata_k<-metadata[idx, ]
  dm_k<-dm[idx, idx]
  cat("The number of kept samples (after filtering out samples with NA values in ",target_field,"): ", nrow(metadata_k) ,"\n")
  result<-list()
  result$dm<-dm_k
  result$metadata<-metadata_k
  result
}

filter_samples_by_sample_ids_in_metadata <- function(data, metadata, ids_col=NA){
  if(is.na(ids_col)){
    shared_ids<-intersect(rownames(data), rownames(metadata))
    metadata_idx<-which(rownames(metadata) %in% shared_ids)
  }else{
    shared_ids<-intersect(rownames(data), metadata[, ids_col])
    metadata_idx<-which(metadata[, ids_col] %in% shared_ids)
  }
  data_matched<-data[shared_ids, ]
  data_matched<-data_matched[order(rownames(data_matched)),]
  cat("The number of samples in feature table (after filtering out samples with no metadata): ",
      nrow(data_matched) ,"\n")
  metadata_matched<-metadata[metadata_idx, ]

  cat("The number of samples metadata (after filtering out samples with no metadata): ",
      nrow(metadata_matched) ,"\n")

  if(is.na(ids_col)){
    metadata_matched<-metadata_matched[order(rownames(metadata_matched)),]
    cat("The sample IDs are idenical in feature table and metadata: ",
        identical(rownames(data_matched), rownames(metadata_matched)), "\n")
  }else{
    metadata_matched<-metadata_matched[order(as.character(metadata_matched[, ids_col])),]
    cat("The sample IDs are idenical in feature table and metadata: ",
        identical(rownames(data_matched), as.character(metadata_matched[, ids_col])), "\n")
  }


  result<-list()
  result$data<-data_matched
  result$metadata<-metadata_matched
  return(result)
}



filter_samples_by_seq_depth<-function(data, metadata, cutoff=1000){
  seq_dep_idx<-which(rowSums(data) > cutoff)
  cat("The number of kept samples with more than ",cutoff," reads: ", length(seq_dep_idx), "\n")
  if(length(seq_dep_idx)>0){
    metadata_k<-metadata[seq_dep_idx, ]
    data_k<-data[seq_dep_idx, ]
  }else{
    stop("The read count of all samples is less than sequencing depth threshold!")
  }
  cat("The sample IDs are idenical in feature table and metadata: ", identical(rownames(data_k), rownames(metadata_k)), "\n")
  result<-list()
  result$data<-data_k
  result$metadata<-metadata_k
  return(result)
}

convert_y_to_numeric<-function(y, reg=TRUE){
  if(reg & !is.numeric(y)){
    y=as.numeric(as.character(y))
  }else if(!reg){
    y=factor(y)
  }else{
    y=y
  }
  y
}

keep_shared_features<-function(train_x, test_x){
  common_idx<-intersect(colnames(train_x), colnames(test_x))
  train_x_shared<-train_x[, common_idx]
  test_x_shared<-test_x[, common_idx]
  cat("Number of features kept:", length(common_idx), "\n")
  cat("The proportion of commonly shared features in train and test dataset respectively: \n")
  cat("Train data: ", length(common_idx)/ncol(train_x), "\n")
  cat("Test data: ", length(common_idx)/ncol(test_x), "\n")
  result<-list()
  result$train_x_shared<-train_x_shared
  result$test_x_shared<-test_x_shared
  result
}

#' @title harmonize_features
#' @description Feature harmonization across multiple datasets: restrict each feature table
#' (samples in rows, features in columns) to the features shared by all datasets, in a common column order.
#' Features unique to individual datasets are excluded because they cannot contribute to cross-dataset prediction.
#' Pre-harmonized tables (e.g., after aggregation to a shared taxonomic level) can be supplied instead.
#' @param df_list A (named) list of data.frames or matrices with feature IDs as column names.
#' @param verbose A boolean value indicating if the number of shared features and
#' the proportion of features kept in each dataset are printed.
#' @return A list of data.frames with identical columns.
#' @examples
#' df_list <- list(A=data.frame(t1=1:3, t2=4:6, t3=7:9),
#'                 B=data.frame(t3=1:4, t1=2:5, t4=3:6))
#' harmonize_features(df_list)
#' @export
harmonize_features <- function(df_list, verbose=TRUE){
  if(!is.list(df_list) || is.data.frame(df_list) || length(df_list) < 2)
    stop("df_list should be a list of at least two feature tables.")
  feature_sets <- lapply(df_list, colnames)
  if(any(vapply(feature_sets, is.null, logical(1))))
    stop("All feature tables should have feature IDs as column names.")
  shared <- Reduce(intersect, feature_sets)
  if(length(shared)==0) stop("No feature is shared by all datasets.")
  if(verbose){
    cat("Number of features shared by all datasets:", length(shared), "\n")
    cat("The proportion of features kept in each dataset:\n")
    print(round(vapply(feature_sets, function(f) length(shared)/length(f), numeric(1)), 3))
  }
  lapply(df_list, function(d) data.frame(d[, shared, drop=FALSE], check.names=FALSE))
}


#' @noRd
add_ann<-function(tab, fmetadata, tab_id_col=1, fmetadata_id_col=1){
  fmetadata[, fmetadata_id_col]<-as.character(fmetadata[, fmetadata_id_col])
  tab[, tab_id_col]<-as.character(tab[, tab_id_col])
  matched_idx<-which(fmetadata[, fmetadata_id_col] %in% tab[, tab_id_col])
  uniq_features_len<-length(unique(tab[, tab_id_col]))
  if(uniq_features_len>length(matched_idx)){
    warning("# of features has no matching IDs in the taxonomy file: ", uniq_features_len-length(matched_idx), "\n")
  }
  fmetadata_matched<-fmetadata[matched_idx,]
  out<-merge(tab, fmetadata_matched, by.x=tab_id_col, by.y=fmetadata_id_col)
  out
}

rbind_na<-function(l){
  max_len<-max(unlist(lapply(l, length)))
  c_l<-lapply(l, function(x) {c(x, rep(NA, max_len - length(x)))})
  do.call(rbind, c_l)
}
expand_Taxon<-function(df, Taxon){
  taxa_df <- rbind_na(strsplit(as.character(df[, Taxon]), '; '))
  colnames(taxa_df) <- c("kingdom","phylum","class","order","family","genus","species") #"kingdom",
  data.frame(df, taxa_df)
}
