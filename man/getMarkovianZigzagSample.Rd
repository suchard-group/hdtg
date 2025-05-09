% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSampleFunctions.R
\name{getMarkovianZigzagSample}
\alias{getMarkovianZigzagSample}
\title{Draw one Markovian zigzag sample}
\usage{
getMarkovianZigzagSample(position, velocity = NULL, engine, travelTime)
}
\arguments{
\item{position}{a d-dimensional position vector.}

\item{velocity}{optional d-dimensional velocity vector. If NULL, it will be generated within the function.}

\item{engine}{an object representing the Markovian zigzag engine, typically containing settings and state required for the simulation.}

\item{travelTime}{the duration for which the dynamics are simulated.}
}
\value{
A list containing the position and velocity after simulating the dynamics.
}
\description{
Simulate the Markovian zigzag dynamics for a given position over a specified travel time.
}
