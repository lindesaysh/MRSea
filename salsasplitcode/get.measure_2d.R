#' Function to get the statistic required to determine fit of the models
#' 
#' @author Cameron Walker, Department of Enginering Science, University of Auckland.
#' 
#' @export
#' 

###################################################################################

"get.measure_2d" <- function(fitnessMeasure,measures,out.lm, data, dists,aR,radii,radiusIndices){
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
  
  if(fitnessMeasure=="QAIC"){       
    
    if(baseModel$family[1]=="quasipoisson"){
      PoisMod<-update(out.lm, round(.)~., family=poisson)
      fitStat <- QAIC(PoisMod, chat = summary(out.lm)$dispersion)}
    
    if(baseModel$family[1]=="quasibinomial"){
      BinMod<-update(out.lm, family=binomial)
      fitStat <- QAIC(BinMod, chat = summary(out.lm)$dispersion)}
  }
  
  if(fitnessMeasure=="QAICc"){       
    if(baseModel$family[1]=="quasipoisson"){
      PoisMod<-update(out.lm, family=poisson)
      fitStat <- QAICc(PoisMod, chat = summary(out.lm)$dispersion)}
    
    if(baseModel$family[1]=="quasibinomial"){
      BinMod<-update(out.lm, family=binomial)
      fitStat <- QAICc(BinMod, chat = summary(out.lm)$dispersion)}
    
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
    
    nfolds<-length(unique(data$foldid))
    store<- matrix(0, nrow=nfolds, ncol=1)
    d2k<-dists
    for(f in 1:nfolds){
      dists<-d2k[data$foldid!=f,]
      foldedFit<- update(out.lm, .~., data=data[data$foldid!=f,])
      if(length(coef(foldedFit))==1){
        dists<-d2k[data$foldid==f,]
        predscv<- exp(as.matrix(model.matrix(out.lm)[data$foldid==f])%*%coef(foldedFit)) * exp(baseModel$offset)[data$foldid==f]
      }else{  
        dists<-d2k[data$foldid==f,]
        predscv<- exp(model.matrix(out.lm)[data$foldid==f,]%*%coef(foldedFit)) * exp(out.lm$offset)[data$foldid==f]
      }        
      store[f]<- mean((data$response[data$foldid==f]-predscv)**2)
    }
    
    fitStat<- mean(store)
  }
  
  # calculate a QIC with bayesian penalty
  if(fitnessMeasure=="QICb"){       
    fitStat <- QICb(out.lm$data$response, fitted(out.lm), length(out.lm$coeff), length(out.lm$data$response))
  }
  
  # cat("Evaluating new fit: ", fitStat, "\n")
  if(is.na(fitStat)){
    fitStat<- tempMeasure + 1000
    cat("Change Fit due to fitStat=NA: ", fitStat, "\n")
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
