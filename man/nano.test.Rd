% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AClass_230815_2.R
\name{nano.test}
\alias{nano.test}
\title{Test classification model on NanoString or other transcriptomic expression data}
\usage{
nano.test(
  prefix,
  training_model_obj,
  data,
  alg_list = NULL,
  min_test_features = NULL,
  max_test_features = NULL,
  out_path = NULL
)
}
\arguments{
\item{prefix}{Output file prefix (e.g., \code{paste(project, format(Sys.time(), "\%Y-\%m-\%d_\%H\%M"), sep = "_")}).}

\item{training_model_obj}{Object returned from \code{nano.train()}, must include \verb{$training_model_list}.}

\item{data}{AClass object containing test data in \code{norm.t}. If known class labels (\code{Group}) are present, a confusion matrix is computed to evaluate performance.}

\item{alg_list}{Optional. Subset of algorithms to test. If not provided, all algorithms in the training object are used.}

\item{min_test_features}{Minimum number of features to test. If not specified, inferred from training model.}

\item{max_test_features}{Maximum number of features to test. If not specified, inferred from training model.}

\item{out_path}{Output path. If not provided, it will be inferred from \code{run_info}.}
}
\value{
AClass object with \verb{$test_results_full} and updated \code{run_info}. If \code{Group} is present, also outputs confusion matrix results to disk.
}
\description{
Applies trained models from \code{nano.train()} to test data for classification. If known class labels (\code{Group}) are present, a confusion matrix is computed to evaluate performance. Prediction results are saved to disk and appended to the AClass object for downstream analysis.
}
\examples{
test.obj <- nano.test(prefix = "demo", training_model_obj = trained_models, data = test.obj)
}
