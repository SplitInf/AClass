
# AClass

`AClass` is an R package designed to perform tumor subgroup classification using transcriptomic data, particularly from small cohorts. The method was originally developed for the classification of Atypical Teratoid Rhabdoid Tumors (ATRT) using probe-based NanoString nCounter data, but the core framework is compatible with other gene expression platforms.

The motivation behind AClass is to enable robust subgroup classification in clinical settings when working with limited sample size, degraded/archival RNA, or cost- and time-sensitive workflows. To address this, AClass implements an ensemble classification strategy over a minimal biomarker panel.


#### Method Summary
The AClass workflow consists of the following components:
- Subgroup signature gene list: A minimal gene signature list, optimized for classification performance. The default panel, [CodeSet30](https://github.com/SplitInf/AClass/blob/main/inst/extdata/probes_list_hgnc.txt)
, was derived from subgroup-specific overexpressed genes in ATRT and is included in the package.
- Model selection: Multiple machine learning algorithms (e.g., rf, pam, glmnet, nb, knn) are evaluated across a range of gene subset sizes (e.g., 20–30 genes). As described in the publication, the top 5 algorithms over 20-30 genes were used to create a total of 55 classification models. Pre-trained models are included with the package.
- Ensemble prediction: Classification is performed using a hard-voting strategy across all selected models. The final class assignment is based on majority vote, and a prediction score is computed by averaging probabilities across models that agree on the predicted class.

#### Input Requirements
- Expression data matrix (NanoString or other transcriptomic platform)
- Training labels (required only for retraining)
- Pre-ranked probe list (required only for retraining)

#### Output
- Predicted class labels and prediction scores
- Classification report summarizing model agreement and QC
- Visualization tools including multidimensional scaling (MDS) plots


## Installation

``` r
#### Using devtools in R:
Using the R package `devtools`, run
`devtools::install_github('https://github.com/SplitInf/AClass')`

#### From source:
Clone the repository: ` https://github.com/SplitInf/AClass.git`
Open R in the directory you cloned the package in and run `install.packages('AClass', repos = NULL)`

```

## Tutorial

For detailed usage instructions, see the AClass [tutorial](http://htmlpreview.github.io/?https://github.com/SplitInf/AClass/blob/main/doc/tutorial.html).



## Publication

If you have used AClass or CodeSet30 for your work, please cite [publication](https://doi.org/10.1093/noajnl/vdae004).
Ben Ho , Anthony Arnoldo , Yvonne Zhong , Mei Lu , Jonathon Torchia , Fupan Yao , Cynthia Hawkins , Annie Huang. **Rapid, economical diagnostic classification of ATRT molecular subgroup using NanoString nCounter platform**. Neuro-Oncology Advances (2024), https://doi.org/10.1093/noajnl/vdae004 .
