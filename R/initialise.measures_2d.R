initialise.measures_2d<- function(knotDist,maxIterations,gap,radii,dists,explData,startKnots, knotgrid, response, baseModel,radiusIndices, initialise, initialKnots, initialaR, fitnessMeasure, interactionTerm, data, knot.seed, initDisp, cv.opts, basis, printout){
  
    attributes(baseModel$formula)$.Environment<-environment()
    baseModel<-update(baseModel, data=data)
    splineParams<-baseModel$splineParams
  
  if(printout){
    print("******************************************************************************")
    print("Initialising...")
    print("******************************************************************************")
  }
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
    
    if(printout){
        print("Space-filling knots....")
    }
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
        warning("WARNING: less knots positioned than desired")
    }
    
    baseModel$splineParams[[1]]$initialKnots <- initialKnots
    
  }else{
    if(length(initialaR)>0){
      posKnots<-initialaR
    }else{
      numNeeded = nrow(initialKnots)
      knots2<-rbind(initialKnots, knotgrid)
      #posKnots<-which(duplicated(knots2)[(nrow(initialKnots)+1):nrow(knots2)]==T)
      # find the nearest knots to the initial locations and remove duplicates
      kd <- as.matrix(dist(knots2))[(nrow(initialKnots) +1) : nrow(knots2), 1:nrow(initialKnots)]
      posKnots <- unique(as.vector(apply(kd,2 ,function(x) which(x == min(x)))))
    }
    
  } # end of ifelse statement related to initialise  = T/F
  
 
  ########### check initialise conforms to gap ##################
  for(i in 1:length(posKnots)){
    toocloseid<-which(knotDist[,posKnots[i]] > 0 & knotDist[,posKnots[i]]<gap)
    if(length(toocloseid)>0){
      posKnots <- posKnots[-toocloseid]
      initialKnots <- initialKnots[-toocloseid]
    }
  }
  
  
  knotPoint<- posKnots
  # print(c('knots: ',knotPoint))
  aR <- knotPoint
  # print(c('actual knots: ',invInd[aR]))
  if(length(radii)>1){
    radiusIndices <-rep((1:length(radii))[(length(radii)/2)],length(aR))  
  }else{
    radiusIndices <-rep(1,length(aR))
  }
  
  
  if (isS4(baseModel)) {
    baseModel@splineParams[[1]]$knotPos<-aR
    baseModel@splineParams[[1]]$radiusIndices<-radiusIndices
    baseModel@splineParams[[1]]$radii<-radii
  } else {
    baseModel$splineParams[[1]]$knotPos<-aR
    baseModel$splineParams[[1]]$radiusIndices<-radiusIndices
    baseModel$splineParams[[1]]$radii<-radii
  }
  
  if(printout){
    print("Initialising model...")
  }
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
        fitStat <- MuMIn::QAIC(PoisMod, chat = initDisp)}
      if(baseModel$family[1]=="quasibinomial"){
        BinMod<-update(baseModel, family=binomial)
        fitStat <- MuMIn::QAIC(BinMod, round(.)~., chat = initDisp)}
  }
  
  if(fitnessMeasure=="QAICc"){
      if(baseModel$family[1]=="quasipoisson"){
        PoisMod<-update(baseModel, family=poisson)
        fitStat <- MuMIn::QAICc(PoisMod, chat = initDisp)}
      if(baseModel$family[1]=="quasibinomial"){
        BinMod<-update(baseModel, family=binomial)
        fitStat <- MuMIn::QAICc(BinMod, chat = initDisp)}
  }
  
  if(fitnessMeasure=="QBIC"){
      if(baseModel$family[1]=='quasipoisson'){
        fitStat <- MuMIn::QAIC(update(baseModel,  round(response) ~ ., family=poisson), chat = initDisp, k=log(nrow(baseModel$data)))
      }
      if(baseModel$family[1]=='quasibinomial'){
        fitStat <- MuMIn::QAIC(update(baseModel, family=binomial), chat = initDisp,k=log(nrow(baseModel$data)))
    }
  }
  
  # fitStat <- get.measure_2d(fitnessMeasure, NULL, baseModel)
  if(fitnessMeasure=="CV.offset"){
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
  
  if(fitnessMeasure=="cv.gamMRSea"){
      fitStat<-cv.gamMRSea(data, baseModel, K=cv.opts$K, cost=cv.opts$cost, s.eed = cv.opts$cv.gamMRSea.seed)$delta[2]
  }
  
  if(fitnessMeasure=="AICtweedie"){
    fitStat<-tweedie::AICtweedie(baseModel)
  }
  
  if(fitnessMeasure=="BICtweedie"){
    fitStat<-tweedie::AICtweedie(baseModel, k=log(nrow(baseModel$data)))
  }
  
  #cat("Evaluating new fit: ", fitStat, "\n")
  if(is.na(fitStat)){
    # fitStat <- fitStat + 10000000
    fitStat <- 10000000
    if(printout){
      cat("Change Fit due to NA: ", fitStat, "\n")  
    }
  }
  
  if(getDispersion(baseModel)>initDisp){
    fitStat<- tempMeasure + 10000000
    if(printout){
      cat("Change Fit due to large dispersion: ",getDispersion(out.lm), ', init: ', initDisp, "\n")
    }
  }

  output = fit.thinPlate_2d(fitnessMeasure, dists,aR,radii, baseModel,radiusIndices,models, fitStat, interactionTerm, data, initDisp, cv.opts, basis)
  out.lm<-output$currentModel
  models<-output$models
  if(printout){
    print("Initial model fitted...")
  }
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
  if(printout){
    print("Fitting Initial Radii")  
  }
  out<-choose.radii(BIC,1:length(radiusIndices),radiusIndices,radii,out.lm,dists,aR,baseModel,fitnessMeasure,response,models, interactionTerm, data, initDisp, cv.opts, basis, printout)
  BIC=out$BIC
  radiusIndices=out$radiusIndices
  out.lm=out$out.lm
  models = out$models
  
  if(printout){
    print("initialising complete")
  }
  
  list(point=point,knotPoint=knotPoint,position=position,aR=aR,BIC=BIC,track=track,out.lm=out.lm, radiusIndices=radiusIndices,models=models)
  
}

