initialise.measures_2d<- function(knotDist,maxIterations,gap,radii,dists,explData,startKnots, knotgrid, response, baseModel,radiusIndices, initialise, initialKnots, initialaR, fitnessMeasure, interactionTerm, data, knot.seed, initDisp, cv.opts, basis, hdetest=FALSE, minKnots=1){
  
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
  # # recall: gridResp has x coordinates of possible knot locations
  # #         or NA if knot location is outside legal region
  # oInd<- 1:length(gridResp)
  # resp<- na.omit(gridResp)
  # 
  # # determine  index of legal knotpoints in gridResp (and hence gridData)
  # if (length(gridResp) > length(resp)) {
  #   mapInd<- oInd[-na.action(resp)]
  # } else {
  #   mapInd<- oInd
  # }
  # # make pointer to determine where in mapInd each legal knot is stored - else = 0
  # # e.g. if 693 grid points but only 500 legal positions, if the 300th knot poistion is
  # # the 200th legal position, then invInd[300] = 200
  # invInd <- rep(0,length(gridResp))
  # for (i in 1:length(mapInd)) {
  #   invInd[mapInd[i]] <- i
  # }
  
  #if initialise is TRUE:
  if (initialise) {
    require(fields)
    numNeeded = startKnots
    
    set.seed(knot.seed)
    
    options('warn'=-1)
    dupPoints <-paste(knotgrid[,1], knotgrid[,2], sep='E')
    duppointid<-which(duplicated(dupPoints)==T)
    
    if(nrow(knotgrid)<1000){
      
      if(length(duppointid)>0){
        if (ncol(knotgrid)>2 & isS4(baseModel)) {
          knotClasses <- knotgrid[,3]
          classes <- sort(unique(knotClasses))
          n_classes <- length(classes)
          num_per_cl <- round(numNeeded / (n_classes))
          if (num_per_cl < minKnots) {
            num_per_cl <- minKnots
          }
          spacefillresultlist <- vector("list",(n_classes))
          knotpozlist <- vector("list",(n_classes))
          for (cl in 1:(n_classes)){
            class_mask = knotClasses == classes[cl]
            knotpoz <- seq(nrow(knotgrid))[class_mask]
            knotpozlist[[cl]] <- knotpoz
            spacefillresult <- cover.design((knotgrid)[class_mask,][-duppointid,1:2], nd=num_per_cl, nruns=1)
            spacefillresultlist[[cl]] <- spacefillresult
          }
        } else {
          spacefillresult<- cover.design((knotgrid)[-duppointid,], nd=numNeeded, nruns=1)
        }
      }
      if(length(duppointid)==0){
        if (ncol(knotgrid)>2 & isS4(baseModel)) {
          knotClasses <- knotgrid[,3]
          classes <- sort(unique(knotClasses))
          n_classes <- length(classes)
          num_per_cl <- round(numNeeded / (n_classes))
          if (num_per_cl < minKnots) {
            num_per_cl <- minKnots
          }
          spacefillresultlist <- vector("list",(n_classes))
          knotpozlist <- vector("list",(n_classes))
          for (cl in 1:(n_classes)){
            class_mask = knotClasses == classes[cl] 
            knotpoz <- seq(nrow(knotgrid))[class_mask]
            knotpozlist[[cl]] <- knotpoz
            spacefillresult<- cover.design(knotgrid[class_mask,1:2], nd=num_per_cl, nruns=1)
            spacefillresultlist[[cl]] <- spacefillresult
          }
        } else {
          spacefillresult<- cover.design(knotgrid, nd=numNeeded, nruns=1)
        }
      }
      if (ncol(knotgrid)>2 & isS4(baseModel)){
        initialKnots <- c()
        posKnots <- c()
        class_kt <- c()
        for (cl in 1:(n_classes)) {
          iK<-spacefillresultlist[[cl]]$design
          pK<-spacefillresultlist[[cl]]$best.id
          pK <- knotpozlist[[cl]][pK]
          initialKnots <- rbind(initialKnots, iK)
          posKnots <- c(posKnots, pK)
          class_kt <- c(class_kt, rep(cl, num_per_cl))
        }
        # posKnots <- list(posKnots=posKnots, classKnots= class_kt)
      } else {
        initialKnots<-spacefillresult$design
        posKnots<-spacefillresult$best.id
      }
    }else{
      SampledPoints<- sample(1:dim(knotgrid)[1], min(1500, dim(knotgrid)[1]))
      
      if(length(duppointid)>0){
        if (ncol(knotgrid)>2 & isS4(baseModel)) {
          knotClasses <- knotgrid[,3]
          n_classes <- sort(unique(knotClasses))
          num_per_cl <- round(numNeeded / (n_classes))
          if (num_per_cl < minKnots) {
            num_per_cl <- minKnots
          }
          spacefillresultlist <- vector("list",(n_classes))
          knotpozlist <- vector("list",(n_classes))
          for (cl in 1:(n_classes)){
            class_mask = knotClasses == classes[cl] 
            knotpoz <- seq(nrow(knotgrid))[class_mask]
            knotpozlist[[cl]] <- knotpoz
            spacefillresult<- cover.design((knotgrid)[SampledPoints,1:2][class_mask,][-duppointid,], nd=num_per_cl, nruns=1)
            spacefillresultlist[[cl]] <- spacefillresult
          }
        } else {
          spacefillresult<- cover.design((knotgrid)[SampledPoints,][-duppointid,], nd=numNeeded, nruns=1)
        }
      }
      if(length(duppointid)==0){
        if (ncol(knotgrid)>2 & isS4(baseModel)) {
          knotClasses <- knotgrid[,3]
          n_classes <- sort(unique(knotClasses))
          num_per_cl <- round(numNeeded / (n_classes))
          if (num_per_cl < minKnots) {
            num_per_cl <- minKnots
          }
          spacefillresultlist <- vector("list",(n_classes))
          knotpozlist <- vector("list",(n_classes))
          for (cl in 1:(n_classes)){
            class_mask = knotClasses == classes[cl] 
            knotpoz <- seq(nrow(knotgrid))[class_mask]
            knotpozlist[[cl]] <- knotpoz
            spacefillresult<- cover.design((knotgrid)[SampledPoints,1:2][class_mask,], nd=num_per_cl, nruns=1)
            spacefillresultlist[[cl]] <- spacefillresult
          }
        } else {
          spacefillresult<- cover.design((knotgrid)[SampledPoints,], nd=numNeeded, nruns=1)
        }
      }
      if (ncol(knotgrid)>2 & isS4(baseModel)){
        initialKnots <- c()
        posKnots <- c()
        class_kt <- c()
        for (cl in 1:(n_classes)) {
          iK<-spacefillresultlist[[cl]]$design
          pK<-spacefillresultlist[[cl]]$best.id
          pK <- knotpozlist[[cl]][pK]
          initialKnots <- rbind(initialKnots, iK)
          posKnots <- c(posKnots, pK)
          class_kt <- c(class_kt, rep(cl, num_per_cl))
        }
        # posKnots <- list(posKnots=posKnots, classKnots= class_kt)
      } else {
        initialKnots<-spacefillresult$design
        posKnots<-spacefillresult$best.id
      }
    }
    
    options('warn'=0)
    
    if (dim(initialKnots)[1]<numNeeded) {
      print("WARNING: less knots positioned than desired")
    }
    
    if (isS4(baseModel)) {
      baseModel@splineParams[[1]]$initialKnots <- initialKnots
    } else {
      baseModel$splineParams[[1]]$initialKnots <- initialKnots
    }
    
  }else{
    numNeeded = nrow(initialKnots)
    knots2<-rbind(initialKnots, knotgrid)
    if(length(initialaR)>0){
      posKnots<-initialaR
    }else{
      posKnots<-which(duplicated(knots2)[(nrow(initialKnots)+1):nrow(knots2)]==T)
    }
  } 
  
  # print(c('knots: ',knotPoint))
  knotPoint <- posKnots
  aR <- posKnots
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
        fitStat <- QAIC(PoisMod, chat = initDisp)}
      if(baseModel$family[1]=="quasibinomial"){
        BinMod<-update(baseModel, family=binomial)
        fitStat <- QAIC(BinMod, round(.)~., chat = initDisp)}
    }
  }
  
  if(fitnessMeasure=="QAICc"){
    if (isS4(baseModel)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      if(baseModel$family[1]=="quasipoisson"){
        PoisMod<-update(baseModel, family=poisson)
        fitStat <- QAICc(PoisMod, chat = initDisp)}
      if(baseModel$family[1]=="quasibinomial"){
        BinMod<-update(baseModel, family=binomial)
        fitStat <- QAICc(BinMod, chat = initDisp)}
    }
  }
  
  if(fitnessMeasure=="QBIC"){
    if (isS4(baseModel)) {
      stop('Fitness measure not supported for multinomial.  Please use AIC, AICc or BIC')
    } else {
      if(baseModel$family[1]=='quasipoisson'){
        fitStat <- QAIC(update(baseModel,  round(response) ~ ., family=poisson), chat = initDisp, k=log(nrow(baseModel$data)))
      }
      if(baseModel$family[1]=='quasibinomial'){
        fitStat <- QAIC(update(baseModel, family=binomial), chat = initDisp,k=log(nrow(baseModel$data)))
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
    if (isS4(baseModel)) {
      set.seed(cv.opts$cv.gamMRSea.seed)
      fitStat<-cv.vglmMRSea(baseModel, K=cv.opts$K, cost=cv.opts$cost)$delta[2]
      #fitStat<-Inf
    } else {
      set.seed(cv.opts$cv.gamMRSea.seed)
      fitStat<-cv.gamMRSea(data, baseModel, K=cv.opts$K, cost=cv.opts$cost)$delta[2]
      #fitStat<-Inf
    }
  }
  
  # calculate accuracy for vglm based multinomial
  if(fitnessMeasure=="mn.accuracy"){ 
    if (isS4(baseModel)) {
      fitStat <- mn.accuracy(baseModel)
    } else {
      stop('Fitness measure only supported for multinomial with vglm')
    }
  }

  #cat("Evaluating new fit: ", fitStat, "\n")
  if(is.na(fitStat)){
    # fitStat <- fitStat + 10000000
    fitStat <- 10000000
    cat("Change Fit due to NA: ", fitStat, "\n")
  }
  if(getDispersion(baseModel)>initDisp){
    fitStat<- tempMeasure + 10000000
    cat("Change Fit due to large dispersion: ",getDispersion(out.lm), ', init: ', initDisp, "\n")
  }

    # output<-fit.thinPlate_2d(fitnessMeasure,dists,invInd[aR],radii,baseModel,radiusIndices,models)
  
  # check for Hauck donner effect 
  if (hdetest) {
    if (isS4(baseModel)){
      hde_check <- hdeff(baseModel)
      if (sum(hde_check) > 0){
        fitStat <- fitStat + 10000000
      }
    }
  }

  output = fit.thinPlate_2d(fitnessMeasure, dists,aR,radii, baseModel,radiusIndices,models, fitStat, interactionTerm, data, initDisp, cv.opts, basis, hdetest)
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
  out<-choose.radii(BIC,1:length(radiusIndices),radiusIndices,radii,out.lm,dists,aR,baseModel,fitnessMeasure,response,models, interactionTerm, data, initDisp, cv.opts, basis,hdetest)
  BIC=out$BIC
  radiusIndices=out$radiusIndices
  out.lm=out$out.lm
  models = out$models
  
  print("initialising complete")
  
  
  ####track <- rbind(track,cbind("init",t(aR),BIC[length(BIC)],adjRsq[length(adjRsq)],GCV[length(GCV)]))
  list(point=point,knotPoint=knotPoint,position=position,aR=aR,BIC=BIC,track=track,out.lm=out.lm, radiusIndices=radiusIndices,models=models)
  
  
  
}


