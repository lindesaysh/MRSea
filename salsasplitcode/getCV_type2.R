#' Function to calculate a Cross validation score
#' 
#' 
#' @export
#' 


getCV_type2<- function(folds, baseModel){
  data<- baseModel$data
  data$response<- baseModel$y
  if(baseModel$family[[1]]=="poisson"){data$response<- round(baseModel$y)}
  data$density <- data$response/exp(baseModel$offset)
  cvscore <- vector(length=folds)
  id_cv<- getCVids(data, folds)
  
  for(i in 1:folds){
    #cat('fold: ', i, '\n')
    tempid<- which(id_cv!=i)
    data2<- data.frame(response=baseModel$y[tempid], model.matrix(baseModel)[tempid,2:length(coefficients(baseModel))], offset=baseModel$offset[tempid])
    names(data2)<- c("response", paste("V", 1:(length(coefficients(baseModel))-1), sep=""), "offset")
    textForEval<- paste("tempCVFit<-glm(response ~ ", paste("V", 1:(length(coefficients(baseModel))-1), sep="", collapse="+"), ", family=family(baseModel), data=data2, offset = offset)")
    eval(parse(text=textForEval))
    # make predictions to new model  
    idpred<- which(id_cv==i)
    newdata<- data.frame(model.matrix(baseModel)[idpred,2:length(coefficients(baseModel))], offset = rep(1, length(idpred)))
    names(newdata)<- c(paste("V", 1:(length(coefficients(baseModel))-1), sep=""), "offset")
    # make predictions to values ==i
    preds<- predict(tempCVFit, newdata, type='response')
    # find cv score between the data - turned into a density using offset and the predictions (calc as a density)
    cvscore[i] <- sum((data$density[idpred] - preds)^2)/length(preds)
  }
  return(cvscore)
}