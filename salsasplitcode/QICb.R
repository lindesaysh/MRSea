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

QICb<-function(data, fits, k, n){
  ql<-sum(data*log(fits) - fits)
  -2*ql + k*log(n)
}
