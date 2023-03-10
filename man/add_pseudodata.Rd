% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nhpp_functions.R
\name{add_pseudodata}
\alias{add_pseudodata}
\title{Add pseudodata}
\usage{
add_pseudodata(Edge.groups, Pseudo.data, Path.data)
}
\arguments{
\item{Edge.groups}{list. A list with groups of edge IDs.}

\item{Pseudo.data}{numeric. A vector with values of pdeusodata.}

\item{Path.data}{numeric. A list of path data obtained from the 'make_data_NHPP_KDE_Markov_kernel' function.}
}
\value{
A list of path data with the pseudodata added.
}
\description{
Adds a vector of pseudodata to the path data obtained from the 'make_data_NHPP_KDE_Markov_kernel' function.
}
\examples{
add_pseudodata(Edge.groups, Pseudo.data, Path.data)

}
\author{
Sergei Tarasov
}
