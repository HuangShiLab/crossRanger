#' Simulated multi-cohort case-control microbiome data
#'
#' A simulated relative abundance table of 60 genera in 120 samples from three cohorts
#' (40, 36 and 44 samples) with two samples (visits) per subject. In all cohorts, cases are enriched
#' in g__Genus01-03 and depleted in g__Genus04-05 (universal markers); cases are also enriched in one
#' cohort-specific genus (g__Genus06, g__Genus07 and g__Genus08 in Cohort_A, Cohort_B and Cohort_C).
#' The cohorts differ in their baseline compositions, and g__Genus10-12 increase with age.
#' @format A list with two data.frames sharing the sample IDs as rownames:
#' \describe{
#'   \item{abundance}{A data.frame of relative abundances: 120 samples by 60 genera.}
#'   \item{metadata}{A data.frame with the cohort, disease_status (Case/Control), subject_id, visit and age of each sample.}
#' }
#' @source Simulated by \code{data-raw/simulate_crossRanger_sim.R} in the source repository.
#' @examples
#' data(crossRanger_sim)
#' table(crossRanger_sim$metadata$cohort, crossRanger_sim$metadata$disease_status)
"crossRanger_sim"
