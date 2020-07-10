# Assumes a vglmMRSea object calculates accuracy from model

mn.accuracy <- function(modelobject){
  response <- modelobject@y
  predictions <- modelobject@fitted.values
  resp_class <- apply(response,1,which.max)
  pred_class <- apply(predictions,1,which.max)
  correct_preds <- pred_class == resp_class
  accuracy <- sum(correct_preds) / length(correct_preds)
  cost <- 1- accuracy
  return(cost)
}
