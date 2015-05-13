#'
#' Function to calculate AICh (Hardin and Hilbe 2013) 
#'
#'
#' @author Lindesay Scott-Hayward, University of St Andrews
#' 
#' @export
#'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~ AICh ~~~~~~~~~~~~~~~~~~~~~~
# function to calculate a quasi version of a bic
# data inputs:
# data: vector of data
# k: number of location and scale parameters estimated
# p: number of model coefficients
# n: number of data points

AICh<-function(model){
  data<-model$data$response
  fits<-fitted(model)
  Ll<-as.numeric(logLik(model))
  p<-length(model$coeff)
  n<-length(data)
  if(model$family$family=='gaussian' | model$family$family=='quasipoisson' | model$family$family=='quasibinomial'){
    k<-2
    if(model$family$family=='quasipoisson'){
      #Ll<-update(model, .~., family=poisson)$aic
      Ll<-sum(data*log(fits) - fits)  
    }
    if(model$family$family=='quasibinomial'){
      #Ll<-update(model, .~., family=binomial)$aic
      Ll<-sum(data*(log(fits/(1-fits))) + log(1-fits))
    }
  }else{
    k<-1  
  }
  
  top<-4*(p^2 - p*k - 2*p)*(p+k+1)*(p+k+2)
  -2*Ll + (top/(n-p-k-2))
}