% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runsTest.R
\name{runsTest}
\alias{runsTest}
\title{Runs Test for Randomness}
\usage{
runsTest(
  y,
  plot.it = FALSE,
  alternative = c("two.sided", "positive.correlated", "negative.correlated"),
  emp.distribution = NULL
)
}
\arguments{
\item{y}{a numeric vector of data values}

\item{plot.it}{logical flag.  If 'TRUE' then the graph will be plotted. If 'FALSE', then it is not plotted.}

\item{alternative}{a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "negative.correlated" or "positive.correlated"}

\item{emp.distribution}{vector containing the empirical distribution of test statistics under the nyll hypothesis. Generated using \code{\link{getEmpDistribution}}.}
}
\value{
A list with the following components.

\code{statistic}	the value of the standardized Runs statistic.
\code{p.value} the p-value for the test.
\code{data.name} a character string giving the names of the data.
\code{alternative} a character string describing the alternative hypothesis.
}
\description{
This function performs the runs test for randomness. Users can choose whether to plot the correlation graph or not, and whether to test against two-sided, negative or positive correlation. NAs from the data are omitted. An empirical distribution may be used for the distribution of the test statistic under the null hypothesis of independence.
}
\details{
On the graph observations which are less than the sample median are represented by letter "A" in red color, and observations which are greater or equal to the sample median are represented by letter "B" in blue color.
}
\examples{

data(ns.data.re)

model<-gamMRSea(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),  
              family='poisson', data=ns.data.re)

runsTest(residuals(model))

## Empirical distribution:

simData<-generateNoise(n=500, response=fitted(model), family='poisson')

empdist<-getEmpDistribution(500, simData, model, data=ns.data.re, plot=FALSE, 
                            returnDist=TRUE,dots=FALSE)
runsTest(residuals(model), emp.distribution=empdist)

}
\references{
Mendenhall, W (1982), Statistics for Management and Economics, 4th Ed., 801-807, Duxbury Press, Boston.

J.L. Gastwirth; Y.R. Gel, W. L. Hui, V. Lyubchich, W. Miao and K. Noguchi (2015). lawstat: Tools for Biostatistics, Public Policy, and Law. R package version 3.0
}
