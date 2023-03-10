% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nhpp_functions.R
\name{normalize_KDE}
\alias{normalize_KDE}
\title{Normalize loess smoothing}
\usage{
normalize_KDE(tree.discr, Maps.mean.loess)
}
\arguments{
\item{tree.discr}{simmap or phylo object. A discretized tree using the 'discr_Simmap' function.}

\item{Maps.mean.loess}{list. A list with the loess smoothing calculated for each edge using the 'loess_smoothing_KDE' function.}
}
\value{
A list with the normalized loess smoothing calculated for each edge.
}
\description{
Normalizes the loess smoothing for the Markov KDE.
}
\examples{
edge.KDE$Maps.mean.loess.norm <- normalize_KDE(tree.discr, edge.KDE$Maps.mean.loess)

}
\author{
Sergei Tarasov
}
