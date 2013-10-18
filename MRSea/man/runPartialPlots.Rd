\name{runPartialPlots}
\alias{runPartialPlots}
\title{Plot partial plots for each of the variables listed in \code{factorlist} or \code{varlist}.}
\usage{
  runPartialPlots(model, data, factorlist = NULL,
    varlist = NULL, showKnots = FALSE, save = FALSE)
}
\arguments{
  \item{model}{Fitted model object (glm or gam)}

  \item{data}{Data frame of data information used to fit
  \code{model}}

  \item{factorlist}{(\code{default=NULL}). Vector or names
  of factor variables}

  \item{varlist}{(\code{default=NULL}). Vector of names of
  continuous variables}

  \item{showKnots}{(\code{default=FALSE}). Logical stating
  whether knot locations should be plotted.}

  \item{save}{(\code{default=FALSE}). Logical stating
  whether plot should be saved into working directory.}
}
\value{
  Partial plots, one for each covariate in
  \code{factorlist} and \code{varlist}
}
\description{
  Plot partial plots for each of the variables listed in
  \code{factorlist} or \code{varlist}.
}
\examples{
#' # load data
data(ns.data.re)

model<-glm(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),
           family='quasipoisson', data=ns.data.re)

runPartialPlots(model, ns.data.re, factorlist=c('floodebb', 'impact'),
                varlist=c('observationhour'))
}

