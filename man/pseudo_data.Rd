% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HeckmanSelect.R
\name{pseudo_data}
\alias{pseudo_data}
\title{Title Pseudo data generation for the Least Square Approximation Method}
\usage{
pseudo_data(selection, outcome, data, Model = c("Normal", "AMH"))
}
\arguments{
\item{selection}{selection equation}

\item{outcome}{outcome equation}

\item{data}{data data matrix containing both the outcome and selection variables}

\item{Model}{can either be Normal error of AMH (Ali-Mikhail-Haq) copula function}
}
\value{
xstar: a square matrix containing all the predictors along with
the correlation and the intercepts

ystar: corresponding pseudo response
}
\description{
Title Pseudo data generation for the Least Square Approximation Method
}
\examples{
pseudo_data(selection, outcome,data=data, Model="AMH")
}
