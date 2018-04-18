#' Function to generate noisy data
#'
#' The function generates a random sample from poisson, overdispersed poisson, binomial and zero inflated binomial samples.
#'
#' @param n number of simulations to generate
#' @param response vector of 'true' means to genereate from
#' @param family one of \code{poisson}, \code{binomial} or \code{zibinomial}
#' @param ... Other parameters required for the family specified
#'
#' @details
#' An additional parameter for the Poisson distribution is the dispersion parameter, specified by d=
#' The additional parameters for the Binomial distribution can be found in \link{rbinom}
#' The zibinomial family requires the \code{VGAM} library to generate zero inflated binomial data. Additional parameters can be found in the help for \link{rzibinom}.
#'
#' @examples
#'
#' data(ns.data.re)
#' 
#' model<-gamMRSea(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),  
#'               family='poisson', data=ns.data.re)
#'
#' simData<-generateNoise(n=500, response=fitted(model), family='poisson')
#'
#' @author LAS Scott-Hayward, University of St Andrews
#' @export
#'
generateNoise<-function(n, response, family, gamma.variance...){
  simData<-matrix(NA, nrow=length(response), ncol=n)
  for(i in 1:n){
    if(family=='poisson'){
      simData[,i]<-rpois.od(length(response), response, ...)
    }
    if(family=='binomial'){
      simData[,i]<-rbinom(n=length(response), prob=response, ...)
    }
    if(family=='zibinomial'){
      simData[,i]<-rzibinom(n=length(response), prob=response, ...)
    }
    if(family=='gaussian'){
      simData[,i]<-rnorm(n=length(response), mean = response, ...)
    }
    if(family=='gamma'){
      mode=response
      sd=sqrt(gamma.variance)
      ra1=(mode + sqrt(mode^2 + 4*sd^2))/(2*sd^2)
      shape1=1+mode*ra1
      simData[,i]<-rgamma(n = length(response), rate=ra1, shape=shape1)
    }
  }
  return(simData)
}


rpois.od<-function(n, lambda, d=1){
  if(d[1]==1)
    rpois(n, lambda)
  else
    rnbinom(n, size=(lambda/(d-1)), mu=lambda)
}


#' qzibinom function from the VGAM package
#'
#'
#' @author VGAM package
#'
#' @export
#' 
qzibinom<-function (p, size, prob, pstr0 = 0)
{
  LLL <- max(length(p), length(size), length(prob), length(pstr0))
  p <- rep_len(p, LLL)
  size <- rep_len(size, LLL)
  prob <- rep_len(prob, LLL)
  pstr0 <- rep_len(pstr0, LLL)
  ans <- p
  ans[p <= pstr0] <- 0
  ans[p > pstr0] <- qbinom((p[p > pstr0] - pstr0[p > pstr0])/(1 -
                                                                pstr0[p > pstr0]), size[p > pstr0], prob[p > pstr0])
  prob0 <- (1 - prob)^size
  deflat.limit <- -prob0/(1 - prob0)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 < 0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[p[ind0] <= pobs0] <- 0
    pindex <- (1:LLL)[ind0 & (p > pobs0)]
    Pobs0 <- pstr0[pindex] + (1 - pstr0[pindex]) * prob0[pindex]
    ans[pindex] <- qposbinom((p[pindex] - Pobs0)/(1 - Pobs0),
                             size = size[pindex], prob = prob[pindex])
  }
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN
  ans
}


rzibinom<-function (n, size, prob, pstr0 = 0)
{
  qzibinom(runif(n), size, prob, pstr0 = pstr0)
}


qposbinom<-function (p, size, prob){
  ans <- qbinom(pbinom(0, size, prob, lower.tail = FALSE) *
                  p + dbinom(0, size, prob), size, prob)
  ans[p == 1] <- size[p == 1]
  ans[p == 0] <- 1
  ans[prob == 0] <- NaN
  ans[p < 0] <- NaN
  ans[1 < p] <- NaN
  ans
}

