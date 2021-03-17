#' dispersion parameter
#' 
#' This function calculates the dispersion parameter for Normal, Binomial, Poisson and Gamma distributions
#' 
#' @param model gamMRSea or glm model object
#' 
#' @details
#' some details
#' 
#' @return 
#' single number of dispersion parameter estimation
#' 
#' @author LAS Scott-Hayward
#' 
#' @export
#' 
#' 
getDispersion<-function(model){
  
    n<-model$df.null + 1
    k<-model$rank
    fits<- fitted(model)
    raw<- model$y
    sqresid<-(raw-fits)^2
    
    if(model$family[[1]]=='gaussian'){
      phi = sum(sqresid) * 1/(n-k)
    }
    if(model$family[[1]]=='poisson' | model$family[[1]]=='quasipoisson'){
      phi = sum(sqresid/fits) * 1/(n-k)
    }
    # if(model$family[[1]]=='binomial' | model$family[[1]]=='quasibinomial'){
    #   phi = sum(sqresid) * (1/n-k)
    # }
    if(model$family[[1]]=='quasibinomial'){
      phi = sum(sqresid) * (1/n-k)
    }
    if(model$family[[1]]=='binomial'){
       phi = 1
    }
    if(model$family[[1]]=='Gamma'){
      phi = sum(sqresid/fits^2) * 1/(n-k)
      
    }
  return(phi)
}
