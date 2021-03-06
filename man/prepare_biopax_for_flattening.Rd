% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_biopax_for_flattening.R
\name{prepare_biopax_for_flattening}
\alias{prepare_biopax_for_flattening}
\title{Prepare BioPAX for Flattening}
\usage{
prepare_biopax_for_flattening(biopax)
}
\arguments{
\item{biopax}{BioPAX object.}
}
\description{
Prepares a BioPAX object for correct flattening:  \cr
1. Applies \code{clean_biopax} function;  \cr
2. For instances of class \code{Complex}, replaces \code{memberPhysicalEntity} with \code{component};  \cr
3. Removes all \code{memberPhysicalEntity}-related instances by replacing 
references to instances with \code{memberPhysicalEntity} property with references to the instances they reference.
}
\author{
Ivan Grishagin
}
