"get.measure" <- function(fitnessMeasure,measures,out.lm, initDisp){
  
  attributes(out.lm$formula)$.Environment<-environment()
  data<-out.lm$data
  if(length(fitnessMeasure)>1){
    extra<-as.numeric(fitnessMeasure[2:length(fitnessMeasure)])
    fitnessMeasure<-as.character(fitnessMeasure[1])
    require(MuMIn)
  }
  
  #print("Getting measure...")
  
  tempMeasure <- measures[1]
  if(fitnessMeasure=="AIC"){       
    fitStat <- AIC(out.lm)}
  
  if(fitnessMeasure=="AICc"){       
    fitStat <- AICc(out.lm)}
  
  if(fitnessMeasure=="BIC"){       
    fitStat <- BIC(out.lm)}
  
  if(fitnessMeasure=="QAIC"){   
    if(out.lm$family[1]=='quasipoisson'){
      fitStat <- QAIC(update(out.lm, round(response) ~ ., family=poisson), chat = extra)
    }else{
      fitStat <- QAIC(out.lm, chat = extra)
    }
    if(out.lm$family[1]=='quasibinomial'){
      fitStat <- QAIC(update(out.lm,.~. ,family=binomial), chat = extra)
    }else{
      fitStat <- QAIC(out.lm, chat = extra)
    }
  }
  
  if(fitnessMeasure=="QBIC"){   
    if(out.lm$family[1]=='quasipoisson'){
      fitStat <- QAIC(update(out.lm,  round(response) ~ ., family=poisson), chat = extra, k=log(nrow(out.lm$data)))
    }else{
      fitStat <- QAIC(out.lm, chat = extra, k=log(nrow(out.lm$data)))
    }
    if(out.lm$family[1]=='quasibinomial'){
      fitStat <- QAIC(update(out.lm, .~.,family=binomial), chat = extra,k=log(nrow(out.lm$data)))
    }else{
      fitStat <- QAIC(out.lm, chat = extra, k=log(nrow(out.lm$data)))
    }
  }
  
  if (fitnessMeasure == "newCrit") {
    fitStat <- mean((residuals(out.lm)/(1-influence(out.lm)$h))**2)
    
  }
  
  
  if(fitnessMeasure=="QAICc"){       
    if(out.lm$family[1]=='quasipoisson'){
      fitStat <- QAICc(update(out.lm,  round(response) ~ ., family=poisson), chat = extra)
    }else{
      fitStat <- QAICc(out.lm, chat = extra)
    }
    if(out.lm$family[1]=='quasibinomial'){
      fitStat <- QAICc(update(out.lm, .~.,family=binomial), chat = extra)
    }else{
      fitStat <- QAICc(out.lm, chat = extra)
    }
  }
  
  # calculates a QIC but using a bic penalty
  if(fitnessMeasure=="QICb"){       
    fitStat <- QICb(out.lm)
  }
  
  # Hardin and Hilbe AIC statistic.  See Hilbe 2014 modelling count data book
  if(fitnessMeasure=="AICh"){
    fitStat<-AICh(out.lm)
  }
  
  # calculates CV using lsh code below - includes offset
  if(fitnessMeasure=="CV.offset"){
    fitStat<- mean(getCV_type2(folds=5, out.lm))}
  # calculates CV using cv.glm function - no offset included
  if(fitnessMeasure=="cv.glm"){
    
    if(dim(model.matrix(out.lm))[2]==1){
      data2<- data.frame(response=out.lm$y)
      textForEval<- "tempCVFit<-glm(response~1, data=data2, family=family(out.lm))" 
    }
    if(dim(model.matrix(out.lm))[2]>1){
      data2<- data.frame(response=out.lm$y, model.matrix(out.lm)[,2:length(coefficients(out.lm))])
      names(data2)<- c("response", paste("V", 1:(length(coefficients(out.lm))-1), sep=""))
      textForEval<- paste("tempCVFit<-glm( response ~ ", paste("V", 1:(length(coefficients(out.lm))-1), sep="", collapse="+"), ", family=family(out.lm), data=data2)")
    }
    eval(parse(text=textForEval))  
    require(boot)
    fitStat<-cv.glm(data2,tempCVFit, K=5)$delta[2]
  }
  
  # calculates cv with blocking structure
  if(fitnessMeasure=="CV"){
    #out.lm<<- out.lm
    store<- matrix(0, nrow=nfolds, ncol=1)
    for(f in 1:nfolds){
      #explanatory<- out.lm$data
      subdata<<- combinedData[combinedData$fold!=f,]
      foldedFit<- glm(out.lm$formula, family=out.lm$family, data = subdata, offset = log(subdata$area))
      preds<- exp(model.matrix(out.lm)[combinedData$fold==f,]%*%coef(foldedFit))*exp(out.lm$offset)[combinedData$fold==f]
      store[f]<- mean((out.lm$data$response[combinedData$fold==f]-preds)**2)
    }
    fitStat<- mean(store)
    #fitStat<- outer.func(out.lm)
  }
  
  if(fitnessMeasure=="CV2"){
    fitStat<-getCV_CReSS(data, out.lm)
  }
  
  #cat("Evaluating new fit: ", fitStat, "\n")
  if(is.na(fitStat)){
    fitStat<- tempMeasure + 1000
    cat("Change Fit due to NA: ", fitStat, "\n")
  }
  if(getDispersion(out.lm)>initDisp){
    if(tempMeasure=='NA' | is.null(tempMeasure) | is.na(tempMeasure)){
      fitStat <- fitStat + 10000
    }else{
      fitStat<- tempMeasure + 1000  
    }
    #cat("Change Fit due to large dispersion: ",getDispersion(out.lm), ', init: ', initDisp, "\n")
  }
  
  list(tempMeasure=tempMeasure,fitStat=fitStat)
}