# Simulate the example dataset `crossRanger_sim`: a multi-cohort case-control microbiome study
# with two samples (visits) per subject, universal and cohort-specific disease markers,
# cohort-specific baseline compositions (batch effects) and an age-associated signal.
# Run from the package root: Rscript data-raw/simulate_crossRanger_sim.R
set.seed(2024)
n_taxa <- 60
taxa <- sprintf("g__Genus%02d", seq_len(n_taxa))
cohorts <- c(Cohort_A=40, Cohort_B=36, Cohort_C=44)
universal_up <- 1:3
universal_down <- 4:5
cohort_specific <- c(Cohort_A=6, Cohort_B=7, Cohort_C=8)
age_taxa <- 10:12
base <- as.numeric(gtools::rdirichlet(1, rep(1, n_taxa)))
base[c(universal_up, universal_down, cohort_specific, age_taxa)] <- 1/n_taxa

abundance_list <- list()
metadata_list <- list()
for(cohort in names(cohorts)){
  n_subjects <- cohorts[[cohort]]/2
  cohort_base <- base * exp(rnorm(n_taxa, 0, 0.5))
  for(s in seq_len(n_subjects)){
    subject_id <- sprintf("%s_S%02d", sub("Cohort_", "", cohort), s)
    status <- if(s %% 2 == 1) "Case" else "Control"
    age <- round(runif(1, 20, 70))
    subject_prob <- cohort_base * exp(rnorm(n_taxa, 0, 0.3))
    for(visit in 1:2){
      prob <- subject_prob * exp(rnorm(n_taxa, 0, 0.2))
      if(status=="Case"){
        prob[universal_up] <- prob[universal_up] * 3
        prob[universal_down] <- prob[universal_down] / 3
        prob[cohort_specific[[cohort]]] <- prob[cohort_specific[[cohort]]] * 4
      }
      prob[age_taxa] <- prob[age_taxa] * exp((age + visit - 46)/20)
      counts <- as.numeric(rmultinom(1, 5000, prob/sum(prob)))
      sample_id <- paste0(subject_id, "_V", visit)
      abundance_list[[sample_id]] <- counts/sum(counts)
      metadata_list[[sample_id]] <- data.frame(cohort=cohort, disease_status=status, subject_id=subject_id,
                                               visit=visit, age=age + visit - 1)
    }
  }
}
abundance <- data.frame(do.call(rbind, abundance_list), check.names=FALSE)
colnames(abundance) <- taxa
metadata <- do.call(rbind, metadata_list)
rownames(metadata) <- rownames(abundance)
metadata$cohort <- factor(metadata$cohort)
metadata$disease_status <- factor(metadata$disease_status, levels=c("Case", "Control"))
metadata$subject_id <- factor(metadata$subject_id)
crossRanger_sim <- list(abundance=abundance, metadata=metadata)
save(crossRanger_sim, file="data/crossRanger_sim.rda", compress="xz")
