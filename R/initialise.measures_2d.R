initialise.measures_2d<- function(knotDist,maxIterations,gap,radii,dists,explData,startKnots, knotgrid, response, baseModel,radiusIndices, initialise, initialKnots, initialaR, fitnessMeasure, interactionTerm, data, knot.seed, initDisp, cv.opts, basis, hdetest=FALSE){
  
  if (isS4(baseModel)) {
    attributes(baseModel@misc$formula)$.Environment<-environment()
    baseModel<-update(baseModel, data=data)
    splineParams<-baseModel@splineParams  
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
    baseModel<-update(baseModel, data=data)
    splineParams<-baseModel$splineParams
  }
  
  print("******************************************************************************")
  print("Initialising...")
  print("******************************************************************************")
  #greedy pick from the 1st grid point
  fuse=0
  
  mapKnotPoint <- c()
  legPos <- c()
  numRand <- 0
  
  track <- cbind()
 
  
  #if initialise is TRUE:
  if (initialise) {
    require(fields)
    numNeeded = startKnots
    
    print("Space-filling knots....")
    set.seed(knot.seed)

    #   if(ncol(knotgrid)>2){
    #     duppointid<-c()
    #     chunks<-unique(knotgrid[,3])
    #     for(z in 1:length(chunks)){
    #       dupPoints <-paste(knotgrid[,1][which(knotgrid[,3]==chunks[z])], knotgrid[,2][which(knotgrid[,3]==chunks[z])], sep='E')
    #       duppointid<-c(duppointid, which(duplicated(dupPoints)==T))
    #     }
    # }else{
      options('warn'=-1)
      dupPoints <-paste(knotgrid[,1], knotgrid[,2], sep='E')
      duppointid<-which(duplicated(dupPoints)==T)
    # }
    
    
    if(nrow(knotgrid)<1000){
      
      if(length(duppointid)>0){
        spacefillresult<- cover.design((knotgrid[,1:2])[-duppointid,], nd=numNeeded, nruns=1)
      }
      if(length(duppointid)==0){
        spacefillresult<- cover.design(knotgrid[,1:2], nd=numNeeded, nruns=1)
      }
      initialKnots<-spacefillresult$design
      posKnots<-spacefillresult$best.id
      
    }else{
      SampledPoints<- sample(1:dim(knotgrid)[1], min(1500, dim(knotgrid)[1]))
      
      #space-fill data (subsample - see line above) to get knot locations
      if(length(duppointid)>0){
        spacefillresult<- cover.design((knotgrid[,1:2])[SampledPoints,][-duppointid,], nd=numNeeded, nruns=1)
      }
      if(length(duppointid)==0){
        spacefillresult<- cover.design((knotgrid[,1:2])[SampledPoints,], nd=numNeeded, nruns=1)
      }
      initialKnots<-spacefillresult$design
      posKnots<-spacefillresult$best.id
    }
    
    options('warn'=0)
    
    if (dim(initialKnots)[1]<numNeeded) {
      print("WARNING: less knots positioned than desired")
    }

    
      baseModel$splineParams[[1]]$initialKnots <- initialKnots

  }else{
    if(length(initialaR)>0){
      posKnots<-initialaR
    }else{
      numNeeded = nrow(initialKnots)
      knots2<-rbind(initialKnots, knotgrid)
      posKnots<-which(duplicated(knots2)[(nrow(initialKnots)+1):nrow(knots2)]==T)
    }
    
  } # end of ifelse statement related to initialise  = T/F
  
 
  knotPoint<- posKnots
  # print(c('knots: ',knotPoint))
  aR <- knotPoint
  # print(c('actual knots: ',invInd[aR]))
  radiusIndices <-rep((1:length(radii))[ceiling(length(radii)/2)],length(aR))
  
  if (isS4(baseModel)) {
    baseModel@splineParams[[1]]$knotPos<-aR
    baseModel@splineParams[[1]]$radiusIndices<-radiusIndices
    baseModel@splineParams[[1]]$radii<-radii
  } else {
    baseModel$splineParams[[1]]$knotPos<-aR
    baseModel$splineParams[[1]]$radiusIndices<-radiusIndices
    baseModel$splineParams[[1]]$radii<-radii
  }
  # baseModel$splineParams[[1]]$mapInd<-mapInd
  # baseModel$splineParams[[1]]$invInd<-invInd
  
  print("Initialising model...")
  models = vector("list",0)
  
  if(fitnessMeasure=="AIC"){
    fitStat <- AIC(baseModel)}
  
  if(fitnessMeasure=="AICc"){
    fitStat <- AICc(baseModel)}
  
  if(fitnessMeasure=="BIC"){
    fitStat <- BIC(baseModel)}

  if (fitnessMeasure == "newCrit") {
    if (isS4(baseModel)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      fitStat <- mean((residuals(baseModel)/(1-influence(baseModel)$h))**2)
    }
  }
  
  if(fitnessMeasure=="QAIC"){
    if (isS4(baseModel)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      if(baseModel$family[1]=="quasipoisson"){
        PoisMod<-update(baseModel, round(.)~., family=poisson)
        fitStat <- MuMIn::QAIC(PoisMod, chat = initDisp)}
      if(baseModel$family[1]=="quasibinomial"){
        BinMod<-update(baseModel, family=binomial)
        fitStat <- MuMIn::QAIC(BinMod, round(.)~., chat = initDisp)}
    }
  }
  
  if(fitnessMeasure=="QAICc"){
    if (isS4(baseModel)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      if(baseModel$family[1]=="quasipoisson"){
        PoisMod<-update(baseModel, family=poisson)
        fitStat <- MuMIn::QAICc(PoisMod, chat = initDisp)}
      if(baseModel$family[1]=="quasibinomial"){
        BinMod<-update(baseModel, family=binomial)
        fitStat <- MuMIn::QAICc(BinMod, chat = initDisp)}
    }
  }
  
  if(fitnessMeasure=="QBIC"){
    if (isS4(baseModel)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      if(baseModel$family[1]=='quasipoisson'){
        fitStat <- MuMIn::QAIC(update(baseModel,  round(response) ~ ., family=poisson), chat = initDisp, k=log(nrow(baseModel$data)))
      }
      if(baseModel$family[1]=='quasibinomial'){
        fitStat <- MuMIn::QAIC(update(baseModel, family=binomial), chat = initDisp,k=log(nrow(baseModel$data)))
      }
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
    if (isS4(baseModel)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      fitStat<- mean(getCV_type2(folds = 5, baseModel))
    }
  }
  
  
  if(fitnessMeasure=="CV"){
    if (isS4(baseModel)) {
      fitStat <- getCV_CReSS_2D(data, baseModel, dists,invInd[aR],radii,radiusIndices)
    } else {
      fitStat <- getCV_CReSS_2D(data, baseModel, dists,invInd[aR],radii,radiusIndices)
    }
  }
  #
  
  if(fitnessMeasure=="PRESS"){
    if (isS4(baseModel)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      fitStat <- getPRESS_CReSS(data, baseModel)
    }
  }
  
  
  if(fitnessMeasure=="QICb"){
    if (isS4(baseModel)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      fitStat <- QICb(baseModel)
    }
  }
  
  # Hardin and Hilbe AIC statistic.  See Hilbe 2014 modelling count data book
  if(fitnessMeasure=="AICh"){
    if (isS4(baseModel)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      fitStat<-AICh(baseModel)
    }
  }
  
  if(fitnessMeasure=="cv.gamMRSea"){
      set.seed(cv.opts$cv.gamMRSea.seed)
      fitStat<-cv.gamMRSea(data, baseModel, K=cv.opts$K, cost=cv.opts$cost)$delta[2]
      #fitStat<-Inf
    }
  
 
  #cat("Evaluating new fit: ", fitStat, "\n")
  if(is.na(fitStat)){
    # fitStat <- fitStat + 10000000
    fitStat <- 10000000
    cat("Change Fit due to NA: ", fitStat, "\n")
  }
  
  if(baseModel$splineParams[[1]]$modelType!='pointProcess'){
    if(getDispersion(baseModel)>initDisp){
    fitStat<- tempMeasure + 10000000
    cat("Change Fit due to large dispersion: ",getDispersion(out.lm), ', init: ', initDisp, "\n")
	}
	}
  
  output = fit.thinPlate_2d(fitnessMeasure, dists,aR,radii, baseModel,radiusIndices,models, fitStat, interactionTerm, data, initDisp, cv.opts, basis)
  out.lm<-output$currentModel
  models<-output$models
  print("Initial model fitted...")
  #point <- mapInd[-posKnots]
  #point <- mapInd[-invInd[posKnots]]
  point <- (1:nrow(knotgrid))[-posKnots]
  
  position<- cbind()
  
  for (j in 1:length(knotPoint)) {
    position[knotPoint[j]]<- 0
  }
  for (j in 1:length(point)) {
    position[point[j]]<- j
  }
  
  measures = 0
  
  BIC<-output$fitStat
  #BIC<-get.measure_2d(fitnessMeasure,measures,out.lm, data,  dists, aR,radii,radiusIndices, initDisp)$fitStat
  
  
  #print(BIC[length(BIC)])
  
  print("Fitting Initial Radii")
  out<-choose.radii(BIC,1:length(radiusIndices),radiusIndices,radii,out.lm,dists,aR,baseModel,fitnessMeasure,response,models, interactionTerm, data, initDisp, cv.opts, basis)
  BIC=out$BIC
  radiusIndices=out$radiusIndices
  out.lm=out$out.lm
  models = out$models
  
  print("initialising complete")
  
  
  ####track <- rbind(track,cbind("init",t(aR),BIC[length(BIC)],adjRsq[length(adjRsq)],GCV[length(GCV)]))
  list(point=point,knotPoint=knotPoint,position=position,aR=aR,BIC=BIC,track=track,out.lm=out.lm, radiusIndices=radiusIndices,models=models)
  
  
  
}
