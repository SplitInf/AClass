% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AClass_230815_2.R
\name{nano.trainsplit}
\alias{nano.trainsplit}
\title{Split training data and assign memberships to training data (or testing). Returns index for splitting data.}
\usage{
nano.trainsplit(data, training_memberships_path, N.train.per, seed = NULL)
}
\arguments{
\item{training_memberships_path}{path to file that contains sample name in column 1 and subgroup in column 2}

\item{N.train.per}{\subsection{percentage of training samples use for training, the remaining would be used for validation e.g. 0.8}{
}}

\item{seed}{integer to set seed when needing to wrap function in loop}

\item{train.norm}{Normalized dataframe.}

\item{have_test_label}{If ground truth labels are avaliable for test data}
}
\value{
train.idx
}
\description{
Split training data and assign memberships to training data (or testing). Returns index for splitting data.
}
