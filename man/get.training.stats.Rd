% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AClass_230815_2.R
\name{get.training.stats}
\alias{get.training.stats}
\title{Summarize training configuration}
\usage{
get.training.stats(training_model_obj)
}
\arguments{
\item{training_model_obj}{\itemize{
\item training model object where training performance is pulled out from
}}
}
\description{
Extracts and summarizes key model training settings from a \code{nano.train()} output AClass object, including algorithm type, cross-validation method, number of folds, number of repeats, number of trained models, and feature range tested. Used internally in \code{nano.train.report()}.
}
