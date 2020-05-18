"move.knot_2D" <- function(radii,dists,explData,index,fitnessMeasure,BIC,aR,point,
                           response,knotgrid,out.lm,improve,improveEx,track, maxKnots,tol=0,baseModel,radiusIndices,models, interactionTerm, data, initDisp, cv.opts, basis){
  if (isS4(baseModel)) {
    attributes(baseModel@misc$formula)$.Environment<-environment()
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
  }
  print("******************************************************************************")
  print("Moving knot...")
  print("******************************************************************************")
  # cat('Current Fit in: ', BIC, '\n')
  fitStat<-BIC[1]
  newRadii = radiusIndices
  for (i in 1:length(aR)) {
    tempR<-aR
    tempR[i]<-point[index[1]]
    #       improveR=1
    
    #       while (improveR) {
    #         improveR = 0

    tempRadii = radiusIndices
    
    # give new moved knot the mean radius
    tempRadii[i] = ceiling(length(radii)/2)

    output = fit.thinPlate_2d(fitnessMeasure, dists,tempR,radii,baseModel,tempRadii,models, fitStat, interactionTerm, data, initDisp, cv.opts, basis)
    initModel = output$currentModel
    models = output$models
    #initBIC = get.measure_2d(fitnessMeasure,fitStat,initModel, data,  dists, tempR,radii, tempRadii, initDisp)$fitStat
    initBIC = output$fitStat

    out<-choose.radii(initBIC,1:length(radiusIndices),tempRadii,radii,initModel,dists,tempR,baseModel, fitnessMeasure,response,models, interactionTerm, data, initDisp, cv.opts, basis)

    tempRadii=out$radiusIndices
    tempOut.lm=out$out.lm
    models=out$models
    # output<-get.measure_2d(fitnessMeasure,fitStat,tempOut.lm, data,  dists, tempR,radii, tempRadii, initDisp)
    tempMeasure<-out$BIC

    if (tempMeasure+tol < fitStat) {
      out.lm <- tempOut.lm
      fitStat<-tempMeasure
      print("move ***********************************")
      newR <- tempR
      tempKnot <- i
      improve <- 1
      improveEx <- 1
      newRadii <- tempRadii
      tempindex <- index[1]
    }

  }
  if (length(aR)<maxKnots) {
    print("Adding knot...")
    for(i in 1:length(index)){
      tempR<-c(aR,point[index[i]])  
      tempRadii = c(radiusIndices,(1:length(radii))[ceiling(length(radii)/2)])

    output = fit.thinPlate_2d(fitnessMeasure, dists,tempR,radii,baseModel,tempRadii,models, fitStat, interactionTerm, data, initDisp, cv.opts, basis)
    initModel = output$currentModel
    models = output$models
    # initBIC = get.measure_2d(fitnessMeasure,fitStat,initModel, data,  dists, tempR,radii, tempRadii, initDisp)$fitStat
    initBIC = output$fitStat

    out<-choose.radii(initBIC,1:length(tempRadii),tempRadii,radii,initModel,dists,tempR,baseModel, fitnessMeasure,response,models, interactionTerm, data, initDisp, cv.opts, basis)

    tempRadii=out$radiusIndices
    tempOut.lm=out$out.lm
    models=out$models
    # output<-get.measure_2d(fitnessMeasure,fitStat,tempOut.lm, data,  dists, tempR,radii, tempRadii, initDisp)
    tempMeasure<-out$BIC

    if (tempMeasure +tol < fitStat) {
      out.lm <- tempOut.lm
      fitStat<-tempMeasure
      print("knot added ***********************************")
      newR <- tempR
      tempKnot <- length(aR) + 1
      improve <- 1
      improveEx <- 1
      newRadii <- tempRadii
      tempindex <- index[i]
    }
    }}
  
  # cat('Current Fit out: ', fitStat, '\n')
  if (improve){
    ####track<-rbind(track,cbind("move",t(newR),fitStat,tempaRSQ,tempGCV))
    list(fitStat=fitStat,newR=newR,tempKnot=tempKnot,improve=improve,improveEx=improveEx,track=track, out.lm=out.lm,radiusIndices=newRadii,models=models, index=tempindex)
  }else{list(improve=0,improveEx=0,models=models)}
}



### Hierarchical version

"move.knot_2D.hr" <- function(radii,dists,explData,index,fitnessMeasure,BIC,aR,point,
                           response,knotgrid,out.lm,improve,improveEx,track, maxKnots,tol=0,baseModel,radiusIndices,models, interactionTerm, data, initDisp, cv.opts, basis, subselect, submodel, submodels){
  if (isS4(baseModel)) {
    attributes(baseModel@misc$formula)$.Environment<-environment()
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
  }
  print("******************************************************************************")
  print("Moving knot...")
  print("******************************************************************************")
  # cat('Current Fit in: ', BIC, '\n')
  fitStat<-BIC[1]
  newRadii = radiusIndices
  for (i in 1:length(aR)) {
    tempR<-aR
    tempR[i]<-point[index[1]]
    #       improveR=1
    
    #       while (improveR) {
    #         improveR = 0
    
    tempRadii = radiusIndices
    
    # give new moved knot the mean radius
    tempRadii[i] = ceiling(length(radii)/2)
    
    output = fit.thinPlate_2d.hr(fitnessMeasure, dists,tempR,radii,baseModel,tempRadii,models, fitStat, interactionTerm, data, initDisp, cv.opts, basis, subselect, submodel, submodels)
    initModel = output$currentModel
    models = output$models
    model_sub <- output$subMod
    submodels <- output$submodels
    #initBIC = get.measure_2d(fitnessMeasure,fitStat,initModel, data,  dists, tempR,radii, tempRadii, initDisp)$fitStat
    initBIC = output$fitStat + model_sub$subStat
    
    out<-choose.radii.hr(initBIC,1:length(radiusIndices),tempRadii,radii,initModel,dists,tempR,baseModel, fitnessMeasure,response,models, interactionTerm, data, initDisp, cv.opts, basis, subselect, submodel, submodels, model_sub)
    
    tempRadii=out$radiusIndices
    tempOut.lm=out$out.lm
    models=out$models
    model_sub <- out$model_sub
    submodels <- out$submodels
    # output<-get.measure_2d(fitnessMeasure,fitStat,tempOut.lm, data,  dists, tempR,radii, tempRadii, initDisp)
    tempMeasure<-out$BIC
    
    if (tempMeasure+tol < fitStat) {
      out.lm <- tempOut.lm
      fitStat<-tempMeasure
      print("move ***********************************")
      newR <- tempR
      tempKnot <- i
      improve <- 1
      improveEx <- 1
      newRadii <- tempRadii
      tempindex <- index[1]
    }
    
  }
  if (length(aR)<maxKnots) {
    print("Adding knot...")
    for(i in 1:length(index)){
      print(paste("index-",i,"-", length(index)))
      tempR<-c(aR,point[index[i]])  
      tempRadii = c(radiusIndices,(1:length(radii))[ceiling(length(radii)/2)])
      
      output = fit.thinPlate_2d.hr(fitnessMeasure, dists,tempR,radii,baseModel,tempRadii,models, fitStat, interactionTerm, data, initDisp, cv.opts, basis, subselect,submodel, submodels)
      initModel = output$currentModel
      models = output$models
      model_sub <- output$subMod
      submodels <- output$submodels
      # initBIC = get.measure_2d(fitnessMeasure,fitStat,initModel, data,  dists, tempR,radii, tempRadii, initDisp)$fitStat
      initBIC = output$fitStat + model_sub$subStat
      print(paste("ir", initBIC))
      out<-choose.radii.hr(initBIC,1:length(tempRadii),tempRadii,radii,initModel,dists,tempR,baseModel, fitnessMeasure,response,models, interactionTerm, data, initDisp, cv.opts, basis,subselect,submodel,submodels,model_sub)
      
      tempRadii=out$radiusIndices
      tempOut.lm=out$out.lm
      models=out$models
      model_sub <- out$model_sub
      submodels <- out$submodels
      # output<-get.measure_2d(fitnessMeasure,fitStat,tempOut.lm, data,  dists, tempR,radii, tempRadii, initDisp)
      tempMeasure<-out$BIC
      print(paste("tm", tempMeasure))
      if (tempMeasure +tol < fitStat) {
        out.lm <- tempOut.lm
        fitStat<-tempMeasure
        print("knot added ***********************************")
        newR <- tempR
        tempKnot <- length(aR) + 1
        improve <- 1
        improveEx <- 1
        newRadii <- tempRadii
        tempindex <- index[i]
      }
    }}
  
  # cat('Current Fit out: ', fitStat, '\n')
  if (improve){
    ####track<-rbind(track,cbind("move",t(newR),fitStat,tempaRSQ,tempGCV))
    list(fitStat=fitStat,newR=newR,tempKnot=tempKnot,improve=improve,improveEx=improveEx,track=track, out.lm=out.lm,radiusIndices=newRadii,models=models, index=tempindex,submod=model_sub,submodels=submodels)
  }else{list(improve=0,improveEx=0,models=models, submodels=submodels)}
}
