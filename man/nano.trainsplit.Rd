% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AClass_230815_2.R
\name{nano.trainsplit}
\alias{nano.trainsplit}
\title{Split data into training and testing dataset}
\usage{
nano.trainsplit(data, training_memberships_path, N.train.per, seed = NULL)
}
\arguments{
\item{data}{AClass object with \code{$norm.t}.}

\item{training_memberships_path}{Tab-delimited file with sample IDs (column 1) and group labels (column 2)}

\item{N.train.per}{Proportion of data to use for training (e.g., 0.8 for 80\%).}

\item{seed}{Set to ensure reproducibility. Note it is needed to wrap function in loop}
}
\value{
AClass object with added train.data.main and train.data.validate components
}
\description{
Split normalized training data (NanoString or transcriptomic data) and assign 'Group' label to output
}
\examples{
data.obj <- nano.trainsplit(data.obj, training_memberships_path, 0.8)
}
