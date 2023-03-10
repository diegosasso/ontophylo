% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paramo_functions.R
\name{get_descendants_chars}
\alias{get_descendants_chars}
\title{Get characters that are the descendants of a selected ontology term}
\usage{
get_descendants_chars(ontology, annotations = "auto", terms, ...)
}
\arguments{
\item{ontology}{ontology_index object.}

\item{annotations}{character. Sets which annotations to use: "auto" means automatic annotations, "manual" means manual annotations.
Alternatively, any other list of element containing annotations can be specified.}

\item{terms}{character. IDs of ontology terms for which descendants are queried.}

\item{...}{other parameters for ontologyIndex::get_descendants() function.}
}
\value{
The vector of character IDs.
}
\description{
Returns all characters located (associated) with a given ontology term.
}
\examples{
ontology <- HAO
ontology$terms_selected_id <- list(`CHAR:1` = c("HAO:0000653"), `CHAR:2` = c("HAO:0000653"))
get_descendants_chars(ontology, annotations = "manual", "HAO:0000653")

}
\author{
Sergei Tarasov
}
