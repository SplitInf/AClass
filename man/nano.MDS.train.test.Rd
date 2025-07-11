% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AClass_230815_2.R
\name{nano.MDS.train.test}
\alias{nano.MDS.train.test}
\title{Plot both training (train.data.main) and testing (norm.t) data by merging training and test results.}
\usage{
nano.MDS.train.test(
  prefix,
  train.data,
  test.data,
  colour_code,
  plot_type = c("plot", "ggplot", "ggplot_label", "ggplot_label_test", "plotly"),
  train_ellipse = FALSE,
  memberships = NULL,
  gene_list = NULL,
  omit_sample = NULL,
  prob = NULL
)
}
\arguments{
\item{prefix}{Character string used in plot titles.}

\item{train.data}{AClass object containing normalized training data (must include train.data.main).}

\item{test.data}{AClass object containing normalized test data (must include norm.t and test_results).}

\item{colour_code}{Data frame specifying group-to-colour mappings. Typically set by nano.set.colour().}

\item{plot_type}{Type of plot to generate. Options: "plot", "ggplot", "ggplot_label", "ggplot_label_test", "plotly".}

\item{train_ellipse}{Logical. If TRUE, adds confidence ellipse around training samples. Default is FALSE.}

\item{memberships}{Optional data frame of sample names (as rownames) and subgroup labels. Overrides default assignments from training and test data.}

\item{gene_list}{Optional character vector of genes to include in the plot. Default is NULL (use all).}

\item{omit_sample}{Optional character vector of sample names to omit from the plot. Default is NULL.}

\item{prob}{Optional column name in test_results to use for coloring test samples by prediction score. Default is NULL.}
}
\value{
A plot is printed to the active device.
}
\description{
Requires nano.set.colour() to be run first to determine colour_code. Requires nano.test() to obtain testing results.
Memberships and gene_list are optional.
}
\examples{
nano.MDS.train.test(prefix = "demo", train.data = train, test.data = test, colour_code = train$colour_code)
}
