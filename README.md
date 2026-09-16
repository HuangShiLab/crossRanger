# crossRanger
*This package provide functionalities to perform random forests (RF) meta-analysis of microbiome data*

***
## Introduction
crossRanger: this package enables random forests classification or regression on a single target variable for the multiple datasets and allows users to apply the prebuilt random forests model from one dataset to predict the target variable in other datasets.

[ranger](https://github.com/imbs-hl/ranger) is a fast implementation of random forests (Breiman 2001) or recursive partitioning, particularly suited for high dimensional data. Classification, regression, and survival forests are supported. To meet the meta-analysis requirements we further use `ranger` as a core RF implementation function to develop this package for microbiome-centric analyses. 

This R package basically provides a variety of functions for the RF analyses  within a single or across multiple microbiome datasets. Specifically, this package allows 
* RF classification or regression on a single target variable in a single microbiome dataset with out-of-bag estimation, stratified k-fold, group k-fold or leave-one-group-out cross-validation
* RF classification or regression on a single target variable (such as diseased status or age) stratified by another covariate (such as sex, body sites or timepoints) for a single microbiome dataset
* RF classification or regression on a single covariate for multiple microbiome datasets (such as microbiome studies focused on the same phenotype), with feature harmonization across datasets
* The prediction performance comparisons across multiple RF models (including AUROC, AUPRC, and accuracy etc. for classification models; MAE, MSE, RMSE, MAPE, R squared etc. for regression models)
* the application of multiple RF models to each other and output the prediction performance (such as application of an RF regression model on aging on the female to male cohort and vice versa), and leave-one-dataset-out validation
* recursive feature elimination (RFE) for biomarker discovery
* to output the univariate associations in the microbiome with a single covariate (e.g Wilcoxon rank-sum test on the clr-transformed relative abundance of each microbial feature in microbiome data), and to intersect them with the RF-selected features
* to visualize the results above, produce publication-quality figures and render R Markdown reports

## Authors ##
Shi Huang, UC San Diego 

## Installation ##

The development version is maintained on GitHub and can be downloaded as follows:
``` r 
## install.packages('devtools') # if devtools not installed
devtools::install_github('HuangShiLab/crossRanger', build_vignettes = TRUE)
```

It can also be installed from a source tarball on Windows, macOS and Linux:
``` r
install.packages("crossRanger_0.2.0.tar.gz", repos = NULL, type = "source")
```

## Parallel computation ##

The functions that build models for multiple datasets (`rf_clf.by_datasets`, `rf_reg.by_datasets`,
`rf_clf.comps`, `rf_clf.lodo` and `rf_reg.lodo`) run in parallel through `foreach` and `doParallel`,
which forks the R session on macOS/Linux and starts a PSOCK cluster on Windows. The number of cores
is set with `n_cores` (by default, all but four cores).

If a `foreach` backend with more than one worker is registered before the call, crossRanger uses it
instead of registering its own, so the way of parallelization can be chosen freely:
``` r
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
res <- rf_clf.by_datasets(df, metadata, s_category = "cohort", c_category = "disease_status")
parallel::stopCluster(cl)
```

## Examples ##
* A simulated multi-cohort case-control dataset (two samples per subject)
``` r
library(crossRanger)
data(crossRanger_sim)
df <- crossRanger_sim$abundance
metadata <- crossRanger_sim$metadata
```
* RF modeling on one responsive variable
``` r
rf.out.of.bag(x=df, y=metadata$disease_status, imp_pvalues=FALSE)
```
* RF modeling on one responsive variable stratified by another covariate
(here group k-fold CV keeps the samples of a subject in the same fold)
``` r
res_list <- rf_clf.by_datasets(df, metadata, s_category='cohort', c_category='disease_status',
                               positive_class="Case", nfolds=5, cv_type="group", g_category="subject_id")
res_list$perf_summ
plot_clf_ROC(res_list)
```
* Apply RF models to other datasets, and leave-one-dataset-out validation
``` r 
cross_rf <- rf_clf.cross_appl(rf_model_list=res_list$rf_model_list, 
                              x_list=res_list$x_list, 
                              y_list=res_list$y_list, positive_class="Case")
cross_rf$perf_summ   # Train_data, Test_data, Validation_type, Accuracy, AUROC, AUPRC, Kappa, Sensitivity, ...
plot_cross_appl(cross_rf, metric="AUROC")
lodo <- rf_clf.lodo(df, metadata, s_category='cohort', c_category='disease_status', positive_class="Case")
lodo$perf_summ
```
* Biomarker discovery
``` r
x_A <- df[metadata$cohort=="Cohort_A", ]
y_A <- metadata$disease_status[metadata$cohort=="Cohort_A"]
rfe <- rf_clf.rfe(x_A, y_A, nfolds=5, metric="AUROC", tolerance=0.01)
plot_clf_feature_selection(rfe)$plot
bg_test <- BetweenGroup.test(x_A, y_A, clr_transform=TRUE, p.adj.method="fdr", q_cutoff=0.05)
id_robust_markers(rfe, bg_test)
plot_logfc_heatmap(res_list)$plot
corr_datasets_by_imps(res_list, top_n=5)$feature_consistency
```
* Reports
``` r
crossRanger_report(Within_study=res_list, Cross_application=cross_rf, LODO=lodo, RFE=rfe,
                   output_format="html_document", outdir="crossRanger_results")
```

## References ##
* Breiman, L. (2001). Random forests. Mach Learn, 45:5-32. https://doi.org/10.1023/A:1010933404324.
* Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. https://doi.org/10.18637/jss.v077.i01.

## License ##
All source code freely availale under [GPL-3 License](https://www.gnu.org/licenses/gpl-3.0.en.html). 

## Documentation ##
Each exported function in the package has been documented and we have also written an introductory vignette that is accessible by calling 
``` r
vignette('crossRanger--intro', package='crossRanger')
```

## Bugs/Feature requests ##
I appreciate bug reports and feature requests. Please post to the github issue tracker [here](https://github.com/HuangShiLab/crossRanger/issues). 

## Acknowledgements

 This work is supported by IBM Research AI through the AI Horizons Network. For
 more information visit the [IBM AI Horizons Network website](https://www.research.ibm.com/artificial-intelligence/horizons-network/).
