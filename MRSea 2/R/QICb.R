#'
#' Function to calculate QICb 
#'
#' @author Lindesay Scott-Hayward, University of St Andrews
#' 
#' @export
#'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~ QICb ~~~~~~~~~~~~~~~~~~~~~~
# function to calculate a quasi version of a bic
# data inputs:
# data: vector of data
# fits: vector of fitted values from a model
# k: number of parameters estimated (length(coefficients))
# n: number of data points

QICb<-function(model){
  data<-model$data$response
  fits<-fitted(model)
  k<-length(model$coeff)
  n<-length(data)
  
  if(model$family$family=='poisson' | model$family$family=='quasipoisson'){
    ql<-sum(data*log(fits) - fits)  
  }
  if(model$family$family=='binomial' | model$family$family=='quasibinomial'){
    ql<-sum(data*(log(fits/(1-fits))) + log(1-fits))  
  }
  
  -2*ql + k*log(n)
}