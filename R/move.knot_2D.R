"move.knot_2D" <- function(radii,dists,explData,index,fitnessMeasure,BIC,aR,point,
                           response,knotgrid,out.lm,improve,improveEx,track, maxKnots,tol=0,baseModel,radiusIndices,models, interactionTerm, data, initDisp, fit.opts, basis, printout){
  if (isS4(baseModel)) {
    attributes(baseModel@misc$formula)$.Environment<-environment()
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
  }
  if(printout){
    print("******************************************************************************")
    print("Moving knot...")
    print("******************************************************************************")
  }
    # cat('Current Fit in: ', BIC, '\n')
  fitStat<-BIC[1]
  newRadii = radiusIndices
  for (i in 1:length(aR)) {
    tempR<-aR
    tempR[i]<-point[index[1]]
    tempRadii = radiusIndices
    
    # give new moved knot the mean radius
    tempRadii[i] = ceiling(length(radii)/2)

    output = fit.thinPlate_2d(fitnessMeasure, dists,tempR,radii,baseModel,tempRadii,models, fitStat, interactionTerm, data, initDisp, fit.opts, basis, printout)
    initModel = output$currentModel
    models = output$models
    initBIC = output$fitStat

    out<-choose.radii(initBIC,1:length(radiusIndices),tempRadii,radii,initModel,dists,tempR,baseModel, fitnessMeasure,response,models, interactionTerm, data, initDisp, fit.opts, basis, printout)

    tempRadii=out$radiusIndices
    tempOut.lm=out$out.lm
    models=out$models

    tempMeasure<-out$BIC

    if (tempMeasure+tol < fitStat) {
      out.lm <- tempOut.lm
      fitStat<-tempMeasure
      if(printout){
        print("move ***********************************")
      }
      newR <- tempR
      tempKnot <- i
      improve <- 1
      improveEx <- 1
      newRadii <- tempRadii
      tempindex <- index[1]
    }

  }
  if (length(aR)<maxKnots) {
    if(printout){
      print("Adding knot...")
    }
    for(i in 1:length(index)){
      tempR<-c(aR,point[index[i]])  
      tempRadii = c(radiusIndices,(1:length(radii))[ceiling(length(radii)/2)])

    output = fit.thinPlate_2d(fitnessMeasure, dists,tempR,radii,baseModel,tempRadii,models, fitStat, interactionTerm, data, initDisp, fit.opts, basis, printout)
    initModel = output$currentModel
    models = output$models
    initBIC = output$fitStat

    out<-choose.radii(initBIC,1:length(tempRadii),tempRadii,radii,initModel,dists,tempR,baseModel, fitnessMeasure,response,models, interactionTerm, data, initDisp, fit.opts, basis, printout)

    tempRadii=out$radiusIndices
    tempOut.lm=out$out.lm
    models=out$models

    tempMeasure<-out$BIC

    if (tempMeasure +tol < fitStat) {
      out.lm <- tempOut.lm
      fitStat<-tempMeasure
      if(printout){
        print("knot added ***********************************")
      }
      newR <- tempR
      tempKnot <- length(aR) + 1
      improve <- 1
      improveEx <- 1
      newRadii <- tempRadii
      tempindex <- index[i]
    }
    }}
  if (improve){
    list(fitStat=fitStat,newR=newR,tempKnot=tempKnot,improve=improve,improveEx=improveEx,track=track, out.lm=out.lm,radiusIndices=newRadii,models=models, index=tempindex)
  }else{list(improve=0,improveEx=0,models=models)}
}
