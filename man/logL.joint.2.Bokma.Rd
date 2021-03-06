% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logL.joint.2.Bokma.R
\name{logL.joint.2.Bokma}
\alias{logL.joint.2.Bokma}
\title{Log-likelihoods for evolutionary models}
\usage{
logL.joint.2.Bokma(p, y, gg)
}
\arguments{
\item{p}{parameters of the model to be optimized}

\item{y}{a paleoTS object}

\item{gg}{numeric vector indicating membership of each sample in a segment}
}
\value{
The log-likelihood of the parameter estimates, given the data.
}
\description{
Returns log-likelihood for a the Bokma model with two punctuations.
}
\details{
In general, users will not be access these functions directly, but instead use the optimization functions, which use these functions to find the best-supported parameter values.
}
\author{
Kjetil Lysne Voje
}
