get.measure_2d<- function(fitnessMeasure,measures,out.lm, data, dists,aR,radii,radiusIndices, initDisp, seed.in){
  #print("ooooooooooooooooooooooooooooooooooooooo")
  #print("Getting measure...")
  #print("ooooooooooooooooooooooooooooooooooooooo")
  
  attributes(out.lm$formula)$.Environment<-environment()
  
  tempMeasure <- measures[1]
  if(fitnessMeasure=="AIC"){       
    fitStat <- AIC(out.lm)}
  
  if(fitnessMeasure=="AICc"){       
    fitStat <- AICc(out.lm)}
  
  if(fitnessMeasure=="BIC"){       
    fitStat <- BIC(out.lm)}
  
  if (fitnessMeasure == "newCrit") {
    fitStat <- mean((residuals(out.lm)/(1-influence(out.lm)$h))**2)
    
  }
  
  
  if(fitnessMeasure=="QAIC"){       
    
    if(out.lm$family[1]=="quasipoisson"){
      PoisMod<-update(out.lm, round(.)~., family=poisson)
      fitStat <- QAIC(PoisMod, chat = initDisp)}
    
    if(out.lm$family[1]=="quasibinomial"){
      BinMod<-update(out.lm, family=binomial)
      fitStat <- QAIC(BinMod, chat = initDisp)}
  }
  
  if(fitnessMeasure=="QAICc"){       
    if(out.lm$family[1]=="quasipoisson"){
      PoisMod<-update(out.lm, family=poisson)
      fitStat <- QAICc(PoisMod, chat = initDisp)}
    
    if(out.lm$family[1]=="quasibinomial"){
      BinMod<-update(out.lm, family=binomial)
      fitStat <- QAICc(BinMod, chat = initDisp)}
    
  }
  
  if(fitnessMeasure=="QBIC"){   
    if(out.lm$family[1]=='quasipoisson'){
      fitStat <- QAIC(update(out.lm,  round(response) ~ ., family=poisson), chat = initDisp, k=log(nrow(out.lm$data)))
    }
    if(out.lm$family[1]=='quasibinomial'){
      fitStat <- QAIC(update(out.lm, family=binomial), chat = initDisp,k=log(nrow(out.lm$data)))
    }
  }
  
  if(fitnessMeasure=="cv.offset"){
    #     if(dim(model.matrix(out.lm))[2]==1){
    #       data2<- data.frame(response=response)
    #       textForEval<- "tempCVFit<-glm(response~1, data=data2, family=family(out.lm))" 
    #     }
    #     if(dim(model.matrix(out.lm))[2]>1){
    #       data2<- data.frame(response=response, model.matrix(out.lm)[,2:length(coefficients(out.lm))], offset = exp(baseModel$offset))
    #       names(data2)<- c("response", paste("V", 1:(length(coefficients(out.lm))-1), sep=""), "offset")
    #       textForEval<- paste("tempCVFit<-glm(round(response) ~ ", paste("V", 1:(length(coefficients(out.lm))-1), sep="", collapse="+"), ", family=family(out.lm), data=data2, offset = log(offset))")
    #     }
    #     eval(parse(text=textForEval))  
    #     require(boot)
    #fitStat<-cv.glm(data2,tempCVFit, K=5)$delta[2]
    fitStat<- mean(getCV_type2(folds = 5, out.lm))
  }
  
  if(fitnessMeasure=="CV"){  
    fitStat <- getCV_CReSS_2D(data, out.lm, dists,aR,radii,radiusIndices)
  }
  
  if(fitnessMeasure=="cv.gamMRSea"){  
    set.seed(seed.in)
    fitStat <- cv.gamMRSea(data, out.lm, K=10)$delta[2]
  }
  
  if(fitnessMeasure=="PRESS"){  
    fitStat <- getPRESS_CReSS(data, out.lm)
  }
  
  
  # calculate a QIC with bayesian penalty
  if(fitnessMeasure=="QICb"){       
    fitStat <- QICb(out.lm)
  }
  
  # Hardin and Hilbe AIC statistic.  See Hilbe 2014 modelling count data book
  if(fitnessMeasure=="AICh"){
    fitStat<-AICh(out.lm)
  }
  
  # cat("Evaluating new fit: ", fitStat, "\n")
  if(is.na(fitStat)){
    fitStat<- tempMeasure + 10000000
    cat("Change Fit due to fitStat=NA: ", fitStat, "\n")
  }
  if(getDispersion(out.lm)>initDisp){
    fitStat<- tempMeasure + 10000000
    #cat("Change Fit due to large dispersion: ",getDispersion(out.lm), ', init: ', initDisp, "\n")
  }
  #if(length(which(is.na(out.lm$coefficients)))>0){
  #  fitStat<- tempMeasure + 1000
  #  cat("Change Fit due to NA coefficients: ", fitStat, "\n")
  #}
  
  #print("ooooooooooooooooooooooooooooooooooooooo")
  #print("Got measure...")
  #print("ooooooooooooooooooooooooooooooooooooooo")
  list(tempMeasure=tempMeasure,fitStat=fitStat)
}
