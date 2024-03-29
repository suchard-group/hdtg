% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSampleFunctions.R
\name{getInitialPosition}
\alias{getInitialPosition}
\title{Get an eligible initial value for a MTN with given mean and truncations}
\usage{
getInitialPosition(mean, lowerBounds, upperBounds)
}
\arguments{
\item{mean}{a d-dimensional mean vector.}

\item{lowerBounds}{a d-dimensional vector specifying the lower bounds.}

\item{upperBounds}{a d-dimensional vector specifying the lower bounds.}
}
\value{
an eligible d-dimensional initial vector.
}
\description{
For a given MTN the function returns an initial vector whose elements are one of:
(1) middle point of the truncation interval if both lower and upper bounds are
finite (2) lower (upper) bound +0.1 (-0.1) if only the lower (upper) bound is finite
(3) the corresponding mean value if lower bound = \code{-Inf} are upper bound = \code{Inf}.
}
