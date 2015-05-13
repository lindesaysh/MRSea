initialise.measures_2d<- function(knotDist,maxIterations,gap,radii,dists,gridResp,explData,startKnots,xvals, yvals, explanatory, response, baseModel,radiusIndices, initialise, initialKnots,fitnessMeasure, interactionTerm, data, knot.seed, initDisp){
  
  attributes(baseModel$formula)$.Environment<-environment()
  baseModel<-update(baseModel, data=data)
  splineParams<-baseModel$splineParams
  
  
  print("******************************************************************************")
  print("Initialising...")
  print("******************************************************************************")
  #greedy pick from the 1st grid point
  fuse=0
  
  mapKnotPoint <- c()
  legPos <- c()
  numRand <- 0
  
  track <- cbind()
  # recall: gridResp has x coordinates of possible knot locations
  #         or NA if knot location is outside legal region
  oInd<- 1:length(gridResp)
  resp<- na.omit(gridResp)
  
  # determine  index of legal knotpoints in gridResp (and hence gridData)
  if (length(gridResp) > length(resp)) {
    mapInd<- oInd[-na.action(resp)]       
  } else {
    mapInd<- oInd
  }
  # make pointer to determine where in mapInd each legal knot is stored - else = 0
  # e.g. if 693 grid points but only 500 legal positions, if the 300th knot poistion is
  # the 200th legal position, then invInd[300] = 200
  invInd <- rep(0,length(gridResp))
  for (i in 1:length(mapInd)) {
    invInd[mapInd[i]] <- i
  }
  
  #if initialise is TRUE:
  if (initialise) {
    require(fields)
    numNeeded = startKnots
    numGot = 0
    # while ((numGot < numNeeded) && (fuse < maxIterations)) {
    #  fuse = fuse + 1
    #  legPos = mapInd
    #   posKnots = cbind()
    #  for (i in 1:numNeeded) {
    #     newKnot=legPos[floor(runif(1,1,length(legPos)+1))]
    #     posKnots = c(posKnots,newKnot)
    #     numGot=length(posKnots)
    #     if (length(legPos>1)) {
    #       entries=which(apply(as.matrix(knotDist[invInd[legPos],invInd[posKnots]]),1,min)>=gap)
    #       badEntries=which(apply(as.matrix(knotDist[invInd[legPos],invInd[newKnot]]),1,min)<gap)
    #     } else {
    #       if (any(knotDist[invInd[legPos],invInd[posKnots]]<gap)) {
    #         entries=c();badEntries=invInd[legPos]
    #       } else {
    #         badEntries=c();entries=invInd[legPos]
    #       }
    #     }     
    #     fixFault<-na.omit(match(mapInd[invInd[legPos][badEntries]],legPos))
    #     legPos<-legPos[-fixFault]
    #     if (length(legPos)==0) break
    #  }
    #}
    #if (numGot < numNeeded) print("WARNING: less knots fitted than desired")
    #knotPoint<- posKnots
    #print(knotPoint)
    #aR <- knotPoint
    #radiusIndices <-rep((1:length(radii))[ceiling(length(radii)/2)],length(aR))
    print("Space-filling knots....")
    set.seed(knot.seed)
    
    
    options('warn'=-1)
    
    if(nrow(explData)<1000){
      
      dupPoints <-paste(explData[,1], explData[,2], sep='E')
      
      if(length(which(duplicated(dupPoints)==T))>0){                                
        initialKnots<- cover.design((explData)[-which(duplicated(dupPoints)==T),], nd=numNeeded, nruns=1)$design
      }                      
      if(length(which(duplicated(dupPoints)==T))==0){                                
        initialKnots<- cover.design(explData, nd=numNeeded, nruns=1)$design
      }
      
      
    }else{
      SampledPoints<- sample(1:dim(explData)[1], min(1000, dim(explData)[1]))
      
      #space-fill data (subsample - see line above) to get knot locations
      # remove any duplicated points as doesnt work in cover design
      dupPoints <-paste(explData[SampledPoints,1], explData[SampledPoints,2], sep='E')
      
      if(length(which(duplicated(dupPoints)==T))>0){                                
        initialKnots<- cover.design((explData)[SampledPoints,][-which(duplicated(dupPoints)==T),], nd=numNeeded, nruns=1)$design
      }                      
      if(length(which(duplicated(dupPoints)==T))==0){                                
        initialKnots<- cover.design((explData)[SampledPoints,], nd=numNeeded, nruns=1)$design
      }
    }
    
    options('warn'=0)
    
    posKnots = cbind()
    legPos=mapInd
    for (i in 1:(dim(initialKnots)[1])) {
      new<-scale(explanatory[legPos,],center=c(initialKnots[i,1],initialKnots[i,2]))
      ####Pick nearest grid point that is also far enough away from another knot
      ind<-which.min(abs(new[,1])+abs(new[,2]))
      newKnot = legPos[ind]
      posKnots = c(posKnots,newKnot)
      if (length(legPos>1)) {
        entries=which(apply(as.matrix(knotDist[invInd[legPos],invInd[posKnots]]),1,min)>=gap)
        badEntries=which(apply(as.matrix(knotDist[invInd[legPos],invInd[newKnot]]),1,min)<gap)
      } else {
        if (any(knotDist[invInd[legPos],invInd[posKnots]]<gap)) {
          entries=c();badEntries=invInd[legPos]
        } else {
          badEntries=c();entries=invInd[legPos]
        }
      }     
      fixFault<-na.omit(match(mapInd[invInd[legPos][badEntries]],legPos))
      if(length(fixFault)>0){
        legPos<-legPos[-fixFault]
      }
      if (length(legPos)==0) break
    }
    
    
    
    numGot=length(posKnots)
    
    
    if (numGot < dim(initialKnots)[1]) {
      print("WARNING: less knots positioned than desired")
    }
    knotPoint<- posKnots
    # print(c('knots: ',knotPoint))
    aR <- knotPoint
    # print(c('actual knots: ',invInd[aR]))
    radiusIndices <-rep((1:length(radii))[ceiling(length(radii)/2)],length(aR))
    
  } else {
    posKnots = cbind()
    legPos=mapInd
    for (i in 1:(dim(initialKnots)[1])) {
      new<-scale(explanatory[legPos,],center=c(initialKnots[i,1],initialKnots[i,2]))
      ####Pick nearest grid point that is also far enough away from another knot
      ind<-which.min(abs(new[,1])+abs(new[,2]))
      newKnot = legPos[ind]
      posKnots = c(posKnots,newKnot)
      if (length(legPos>1)) {
        entries=which(apply(as.matrix(knotDist[invInd[legPos],invInd[posKnots]]),1,min)>=gap)
        badEntries=which(apply(as.matrix(knotDist[invInd[legPos],invInd[newKnot]]),1,min)<gap)
      } else {
        if (any(knotDist[invInd[legPos],invInd[posKnots]]<gap)) {
          entries=c();badEntries=invInd[legPos]
        } else {
          badEntries=c();entries=invInd[legPos]
        }
      }     
      fixFault<-na.omit(match(mapInd[invInd[legPos][badEntries]],legPos))
      if(length(fixFault)>0){
        legPos<-legPos[-fixFault]
      }
      if (length(legPos)==0) break
    }
    numGot=length(posKnots)
    if (numGot < length(initialKnots)) {
      print("WARNING: less knots positioned than desired")
      radiusIndices=radiusIndices[1:numGot]
    }
    knotPoint<- posKnots
    #print(c('knots: ',knotPoint))
    aR <- knotPoint
    # print(c('actual knots: ',invInd[aR]))
  } # end of ifelse statement related to initialise  = T/F
  
  
  print("Initialising model...")
  models = vector("list",0)
  
  if(fitnessMeasure=="AIC"){       
    fitStat <- AIC(baseModel)}
  
  if(fitnessMeasure=="AICc"){       
    fitStat <- AICc(baseModel)}
  
  if(fitnessMeasure=="BIC"){       
    fitStat <- BIC(baseModel)}
  
  if (fitnessMeasure == "newCrit") {
    fitStat <- mean((residuals(baseModel)/(1-influence(baseModel)$h))**2)
    
  }
  
  if(fitnessMeasure=="QAIC"){       
    if(baseModel$family[1]=="quasipoisson"){
      PoisMod<-update(baseModel, round(.)~., family=poisson)
      fitStat <- QAIC(PoisMod, chat = initDisp)}
    if(baseModel$family[1]=="quasibinomial"){
      BinMod<-update(baseModel, family=binomial)
      fitStat <- QAIC(BinMod, round(.)~., chat = initDisp)}
  }
  
  if(fitnessMeasure=="QAICc"){       
    if(baseModel$family[1]=="quasipoisson"){
      PoisMod<-update(baseModel, family=poisson)
      fitStat <- QAICc(PoisMod, chat = initDisp)}
    if(baseModel$family[1]=="quasibinomial"){
      BinMod<-update(baseModel, family=binomial)
      fitStat <- QAICc(BinMod, chat = initDisp)}
  }
  
  if(fitnessMeasure=="QBIC"){   
    if(baseModel$family[1]=='quasipoisson'){
      fitStat <- QAIC(update(baseModel,  round(response) ~ ., family=poisson), chat = initDisp, k=log(nrow(baseModel$data)))
    }
    if(baseModel$family[1]=='quasibinomial'){
      fitStat <- QAIC(update(baseModel, family=binomial), chat = initDisp,k=log(nrow(baseModel$data)))
    }
  }
  
  # fitStat <- get.measure_2d(fitnessMeasure, NULL, baseModel)
  if(fitnessMeasure=="CV.offset"){    
    #    if(dim(model.matrix(baseModel))[2]==1){
    #       data2<- data.frame(response=response)
    #       textForEval<- "tempCVFit<-glm(response~1, data=data2)" 
    #             
    #       }
    #     
    #     if(dim(model.matrix(baseModel))[2]>1){
    #       data2<- data.frame(response=response, model.matrix(baseModel)[,2:length(coefficients(baseModel))])
    #       
    #   names(data2)<- c("response", paste("V", 1:(length(coefficients(baseModel))-1), sep=""))
    #   textForEval<- paste("tempCVFit<-glm( response ~ ", paste("V", 1:(length(coefficients(baseModel))-1), sep="", collapse="+"), ", data=data2)")
    #}
    
    #   eval(parse(text=textForEval))  
    #   require(boot)
    #   fitStat<-cv.glm(data2,tempCVFit, K=5)$delta[2]
    # 
    fitStat<- mean(getCV_type2(folds = 5, baseModel))
  }
  
  
  if(fitnessMeasure=="CV"){  
    fitStat <- getCV_CReSS_2D(data, baseModel, dists,invInd[aR],radii,radiusIndices)
  }
  # 
  
  if(fitnessMeasure=="PRESS"){  
    fitStat <- getPRESS_CReSS(data, baseModel)
  }
  
  
  if(fitnessMeasure=="QICb"){       
    fitStat <- QICb(baseModel)
  }
  
  # Hardin and Hilbe AIC statistic.  See Hilbe 2014 modelling count data book
  if(fitnessMeasure=="AICh"){
    fitStat<-AICh(baseModel)
  }
  
  
  #cat("Evaluating new fit: ", fitStat, "\n")
  if(is.na(fitStat)){
    fitStat<- fitStat + 1000
    cat("Change Fit due to NA: ", fitStat, "\n")
  }
  if(getDispersion(baseModel)>initDisp){
    fitStat<- tempMeasure + 1000
    cat("Change Fit due to large dispersion: ",getDispersion(out.lm), ', init: ', initDisp, "\n")
  }
  # output<-fit.thinPlate_2d(fitnessMeasure,dists,invInd[aR],radii,baseModel,radiusIndices,models)
  
  output = fit.thinPlate_2d(fitnessMeasure, dists,invInd[aR],radii, baseModel,radiusIndices,models, fitStat, interactionTerm, data, initDisp)
  
  out.lm<-output$currentModel
  models<-output$models
  print("Initial model fitted...")
  #point <- mapInd[-posKnots]
  point <- mapInd[-invInd[posKnots]]
  
  position<- cbind()
  for (j in 1:length(knotPoint)) {
    position[knotPoint[j]]<- 0
  }    
  for (j in 1:length(point)) {
    position[point[j]]<- j
  }
  
  measures = 0
  
  BIC<-get.measure_2d(fitnessMeasure,measures,out.lm, data,  dists, invInd[aR],radii,radiusIndices, initDisp)$fitStat
  
  
  #print(BIC[length(BIC)])
  
  print("Fitting Initial Radii")
  out<-choose.radii(BIC,1:length(radiusIndices),radiusIndices,radii,out.lm,dists,invInd[aR],baseModel,fitnessMeasure,response,models, interactionTerm, data, initDisp)
  BIC=out$BIC
  radiusIndices=out$radiusIndices
  out.lm=out$out.lm
  models = out$models
  
  print("initialising complete")
  
  
  ####track <- rbind(track,cbind("init",t(aR),BIC[length(BIC)],adjRsq[length(adjRsq)],GCV[length(GCV)]))
  list(point=point,knotPoint=knotPoint,position=position,aR=aR,BIC=BIC,track=track,out.lm=out.lm,invInd=invInd,
       radiusIndices=radiusIndices,models=models)
  
  
  
}