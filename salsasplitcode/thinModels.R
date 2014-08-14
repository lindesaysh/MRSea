#' function to thin the number of models
#' 
#' @author Cameron Walker, Department of Enginering Science, University of Auckland.
#' 
#' @export
#' 

thinModels<- function(models){
  
  bic = vector(length=length(models))
  
  for(i in 1:length(models)){
    bic[i] = models[[i]][[4]]
  }
  
  id = which(bic<=(min(bic)+10))
  tempModels = models[id]
  
  return(models = tempModels)
}


