#' function to thin the number of models
#' 
#' @param models list of models fitted
#' 
#' @author Cameron Walker, Department of Engineering Science, University of Auckland.
#' 
#' @export
#' 

thinModels<- function(models){
  
  if(length(models)>0){
    bic = vector(length=length(models))
    
    for(i in 1:length(models)){
      bic[i] = models[[i]][[4]]
    }
    
    id = which(bic<=(min(bic)+10))
    tempModels = models[id]  
  }else{
    tempModels = models
  }
  return(models = tempModels)
}


