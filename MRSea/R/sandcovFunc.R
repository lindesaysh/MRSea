sandcov<-function (model, id=NULL){
  if(is.null(id)){
    id<-1:length(model$y)
  }
  model$weights <- pmax(1e-08, model$weights)
  l <- length(coef(model))
  vv <- summary.glm(model)$cov.unscaled
  mm <- model.matrix(model)
  ii <- match(colnames(vv), colnames(mm))
  infl <- (residuals(model, "response") * model$prior.weights * 
             mm[, ii]) %*% vv
  jj <- match(colnames(mm), colnames(vv))
  infl <- infl[, jj]
  colnames(infl) <- colnames(mm)
  infl <- rowsum(infl, id)
  t(infl) %*% infl
}