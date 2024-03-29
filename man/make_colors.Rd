% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphics_functions.R
\name{make_colors}
\alias{make_colors}
\title{Make color palette for image plotting}
\usage{
make_colors(Stat, palette)
}
\arguments{
\item{Stat}{numeric. A named vector where values are the statistics, and names are ontology term labels.}

\item{palette}{A character vector or function defining a color palette.}
}
\value{
A character vector where elements are color IDs and names are the input ontology term labels.
}
\description{
Produces a color scale for a given statistic of evolutionary rate.
}
\examples{
stat <- setNames(runif(5, 0.1, 10), c("cranium", "fore_wing", "hind_wing", "pronotum", "propectus") )
hm.palette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"), space = "Lab")
cols.maps <- make_colors(stat, palette = hm.palette(100))
cols.maps

}
\author{
Sergei Tarasov
}
