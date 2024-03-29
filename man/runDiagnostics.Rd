% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runDiagnostics.R
\name{runDiagnostics}
\alias{runDiagnostics}
\title{functions to create observed vs fitted and fitted vs scaled pearsons residual plots}
\usage{
runDiagnostics(model, plotting = "b", save = FALSE, label = NULL)
}
\arguments{
\item{model}{Fitted model object (glm or gam)}

\item{plotting}{Plotting options (\code{default='b'}). \code{b}: returns both plots, \code{f}: returns observed vs fitted only and  \code{r}: returns scale pearsons residual plot only.}

\item{save}{(\code{default=FALSE}). Logical stating whether plot should be saved into working directory. See \code{label} to change directory.}

\item{label}{Character string indicating an label to be added to the plot when using \code{save = TRUE}. Can also include a pathway to a directory of choice.}
}
\value{
Two plots:
\item{Observed vs Fitted}{Plot of observed vs fitted with concordence correlation and marginal R-squared printed in the plot title.}
\item{Fitted vs scaled Pearsons residuals}{The red line is a locally weighted least squares regression line of all of the residuals.}
}
\description{
functions to create observed vs fitted and fitted vs scaled pearsons residual plots
}
\examples{
# load data
data(ns.data.re)

model<-glm(birds ~ observationhour + as.factor(floodebb) + as.factor(impact), 
           family='quasipoisson', data=ns.data.re)

runDiagnostics(model)

}
