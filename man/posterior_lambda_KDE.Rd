% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nhpp_functions.R
\name{posterior_lambda_KDE}
\alias{posterior_lambda_KDE}
\title{Get analytical posterior}
\usage{
posterior_lambda_KDE(tree.list)
}
\arguments{
\item{tree.list}{multiSimmap object.}
}
\value{
A list with mean ($Mean), standard deviation ($SD), and 95HPD interval ($Q_2.5 and $Q_97.5) calculated for the posterior distribution.
}
\description{
Calculates the required statitics for the posterior distribution of number of state changes across all branches of all trees.
}
\examples{
tree.list <- Trees.focal.BRs[[1]]
posterior_lambda_KDE(tree.list)

}
\author{
Sergei Tarasov
}
