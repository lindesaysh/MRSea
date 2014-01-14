\name{getPvalues}
\alias{getPvalues}
\title{Calculate marginal p-values from a \code{model}.}
\usage{
getPvalues(model, varlist = NULL, factorlist = NULL)
}
\arguments{
  \item{model}{Fitted model object (gee)}

  \item{varlist}{(default =\code{NULL}). Vector of
  covariate names (continous covariates only) used to make
  the output table names shorter.  Useful if spline
  parameters are specified in the model.}

  \item{factorlist}{(default =\code{NULL}). Vector of
  covariate names (factor covariates only) used to make the
  output table names shorter. Useful if spline parameters
  are specified in the model.}
}
\value{
Print out table of each variable and its associated
marginal p-value.
}
\description{
An ANOVA is fitted repeatedly with each covariate being the
last so that the output is marginal.  \code{varlist} and
\code{factorlist} are optional and shorten the variable
names in the output.
}
\examples{
# load data
data(ns.data.re)

# make blocking structure
ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear,
                    ns.data.re$DayOfMonth, sep='')
ns.data.re$blockid<-as.factor(ns.data.re$blockid)

initialModel<- geeglm(birds ~ as.factor(floodebb) + as.factor(impact) + observationhour + x.pos +
              y.pos + offset(log(area)), family='poisson',data=ns.data.re, id=blockid)

getPvalues(initialModel, varlist=c('observationhour', 'x.pos', 'y.pos'),
            factorlist=c('floodebb', 'impact'))

getPvalues(initialModel)
}

