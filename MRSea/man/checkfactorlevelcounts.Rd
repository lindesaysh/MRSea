\name{checkfactorlevelcounts}
\alias{checkfactorlevelcounts}
\title{Factor level response check}
\usage{
  checkfactorlevelcounts(factorlist, data, response)
}
\arguments{
  \item{factorlist}{Vector of factor variables specified in
  \code{model}.  Specified so that a check can be made that
  there are non-zero counts in all levels of each factor.}

  \item{data}{Data frame containing columns of covariates
  listed in \code{factorlist}.  Column names must match
  with names in \code{factorlist}}

  \item{response}{A vector of response values}
}
\description{
  This function checks that there are some non-zero counts
  in each level of each factor variable for consideration
  in a model
}
\examples{
# load data
data(ns.data.re)

checkfactorlevelcounts(factorlist=c('floodebb', 'impact'), ns.data.re,
     ns.data.re$birds)
}

