get.measure_2d<- function(fitnessMeasure,measures,out.lm, data, dists,aR,radii,radiusIndices, initDisp, cv.opts, printout){
  
  if (isS4(out.lm)) {
    attributes(out.lm@misc$formula)$.Environment<-environment()
  } else {
    attributes(out.lm$formula)$.Environment<-environment()  
  }
  
  tempMeasure <- measures[1]
  if(fitnessMeasure=="AIC"){
    fitStat <- AIC(out.lm)
  }
  
  if(fitnessMeasure=="AICc"){       
    fitStat <- AICc(out.lm)
  }
  
  if(fitnessMeasure=="BIC"){       
    fitStat <- BIC(out.lm)
  }
  
  if (fitnessMeasure == "newCrit") {
    if (isS4(out.lm)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      fitStat <- mean((residuals(out.lm)/(1-influence(out.lm)$h))**2)
    }
  }
  
  
  if(fitnessMeasure=="QAIC"){
    if (isS4(out.lm)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      if(out.lm$family[1]=="quasipoisson"){
        PoisMod<-update(out.lm, round(.)~., family=poisson)
        fitStat <- MuMIn::QAIC(PoisMod, chat = initDisp)
      }
      if(out.lm$family[1]=="quasibinomial"){
        BinMod<-update(out.lm, family=binomial)
        fitStat <- MuMIn::QAIC(BinMod, chat = initDisp)
      }
    }
  }
  
  if(fitnessMeasure=="QAICc"){
    if (isS4(out.lm)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      if(out.lm$family[1]=="quasipoisson"){
        PoisMod<-update(out.lm, family=poisson)
        fitStat <- MuMIn::QAICc(PoisMod, chat = initDisp)}
      
      if(out.lm$family[1]=="quasibinomial"){
        BinMod<-update(out.lm, family=binomial)
        fitStat <- MuMIn::QAICc(BinMod, chat = initDisp)}
    }
  }
  
  if(fitnessMeasure=="QBIC"){
    if (isS4(out.lm)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      if(out.lm$family[1]=='quasipoisson'){
        fitStat <- MuMIn::QAIC(update(out.lm,  round(response) ~ ., family=poisson), chat = initDisp, k=log(nrow(out.lm$data)))
      }
      if(out.lm$family[1]=='quasibinomial'){
        fitStat <- MuMIn::QAIC(update(out.lm, family=binomial), chat = initDisp,k=log(nrow(out.lm$data)))
      }
    }
  }
  
  if(fitnessMeasure=="cv.offset"){
    if (isS4(out.lm)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      fitStat<- mean(getCV_type2(folds = 5, out.lm))
    }
  }
  
  if(fitnessMeasure=="CV"){ 
    if (isS4(out.lm)) {
      fitStat <- getCV_CReSS_2D(data, out.lm, dists,aR,radii,radiusIndices)
    } else {
      fitStat <- getCV_CReSS_2D(data, out.lm, dists,aR,radii,radiusIndices)
    }
  }
  
  if(fitnessMeasure=="cv.gamMRSea"){
    if (isS4(out.lm)) {
      set.seed(cv.opts$cv.gamMRSea.seed)
      fitStat <- cv.gamMRSea(data, out.lm, K=cv.opts$K, cost=cv.opts$cost)$delta[2]
    } else {
     fitStat <- cv.gamMRSea(data, out.lm, K=cv.opts$K, cost=cv.opts$cost, s.eed = cv.opts$cv.gamMRSea.seed)$delta[2]
    }
  }
  
  if(fitnessMeasure=="PRESS"){  
    if (isS4(out.lm)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      fitStat <- getPRESS_CReSS(data, out.lm)
    }
  }
  
  
  # calculate a QIC with bayesian penalty
  if(fitnessMeasure=="QICb"){ 
    if (isS4(out.lm)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      fitStat <- QICb(out.lm)
    }
  }
  
  # Hardin and Hilbe AIC statistic.  See Hilbe 2014 modelling count data book
  if(fitnessMeasure=="AICh"){
    if (isS4(out.lm)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      fitStat<-AICh(out.lm)
    }
  }
  

  # calculate accuracy for vglm based multinomial
  if(fitnessMeasure=="mn.accuracy"){ 
    if (isS4(out.lm)) {
      fitStat <- mn.accuracy(out.lm)
    } else {
      stop('Fitness measure only supported for multinomial with vglm')
    }
  }

  if(fitnessMeasure=="AICtweedie"){
    fitStat<-tweedie::AICtweedie(out.lm)
  }
  
  if(fitnessMeasure=="BICtweedie"){
    fitStat<-tweedie::AICtweedie(out.lm, k=log(nrow(out.lm$data)))

  }
  
  # cat("Evaluating new fit: ", fitStat, "\n")
  if(is.na(fitStat)){
    fitStat<- tempMeasure + 10000000
    if(printout){
      cat("Change Fit due to fitStat=NA: ", fitStat, "\n")
    }
  }
  
  if(getDispersion(out.lm)>initDisp){
    fitStat<- fitStat + 10000000
  }
  list(tempMeasure=tempMeasure,fitStat=fitStat)
}
