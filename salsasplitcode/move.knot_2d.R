#' Function to move knots around and re-fit model to find best solution
#' 
#' @author Cameron Walker, Department of Enginering Science, University of Auckland.
#' 
#' @export
#' 



"move.knot_2D" <- function(radii,invInd,dists,explData,index,fitnessMeasure,BIC,aR,point,
                           response,explanatory,out.lm,improve,improveEx,track, maxKnots,tol=0,baseModel,radiusIndices,models, interactionTerm, data){
  attributes(baseModel$formula)$.Environment<-environment()
  print("******************************************************************************")
  print("Moving knot...")
  print("******************************************************************************")
  # cat('Current Fit in: ', BIC, '\n')
  fitStat<-BIC[1]
  for (i in 1:length(aR)) {
    tempR<-aR
    tempR[i]<-point[index]
    #       improveR=1
    newRadii = radiusIndices
    #       while (improveR) {
    #         improveR = 0
    
    tempRadii = radiusIndices
    tempRadii[i] = ceiling(length(radii)/2)
    output = fit.thinPlate_2d(fitnessMeasure, dists,invInd[tempR],radii,baseModel,tempRadii,models, fitStat, interactionTerm, data)
    initModel = output$currentModel
    models = output$models
    initBIC = get.measure_2d(fitnessMeasure,fitStat,initModel, data,  dists, invInd[tempR],radii, tempRadii)$fitStat
    #         out<-choose.radii(initBIC,i,tempRadii,radii,initModel,dists,invInd[tempR],baseModel,fitnessMeasure,response)
    out<-choose.radii(initBIC,1:length(radiusIndices),tempRadii,radii,initModel,dists,invInd[tempR],baseModel, fitnessMeasure,response,models, interactionTerm, data)
    tempRadii=out$radiusIndices
    tempOut.lm=out$out.lm
    models=out$models
    output<-get.measure_2d(fitnessMeasure,fitStat,tempOut.lm, data,  dists, invInd[tempR],radii, tempRadii)
    #fitStat<-output$tempMeasure
    tempMeasure<-output$fitStat
    if (tempMeasure+tol < fitStat) {
      out.lm <- tempOut.lm
      fitStat<-tempMeasure
      print("move ***********************************")
      #print(length(as.vector(coefficients(out.lm))))
      #print(tempR)
      #print(fitStat)
      newR <- tempR
      tempKnot <- i
      improve <- 1
      improveEx <- 1
      newRadii = tempRadii
    }
    
  }
  if (length(aR)<maxKnots) {
    tempR<-c(aR,point[index])
    tempRadii = c(radiusIndices,(1:length(radii))[ceiling(length(radii)/2)])
    print("Adding knot...")
    #out<-choose.radii(fitStat,length(radiusIndices)+1,tempRadii,radii,out.lm,dists,invInd[tempR],baseModel,fitnessMeasure,response)
    out<-choose.radii(fitStat,1:length(tempRadii),tempRadii,radii,out.lm,dists,invInd[tempR],baseModel,fitnessMeasure,response,models, interactionTerm, data)
    tempRadii=out$radiusIndices
    tempOut.lm=out$out.lm
    models = out$models
    #tempOut.lm<-fit.thinPlate_2d(dists,invInd[tempR],radii,baseModel,radiusIndices)
    output<-get.measure_2d(fitnessMeasure,fitStat,tempOut.lm, data,  dists, invInd[tempR],radii, tempRadii)
    #fitStat<-output$tempMeasure
    tempMeasure<-output$fitStat
    if (tempMeasure +tol < fitStat) {
      out.lm <- tempOut.lm
      fitStat<-tempMeasure
      print("knot added ***********************************")
      #print(length(as.vector(coefficients(out.lm))))
      #print(tempR)
      #print(fitStat)
      newR <- tempR
      tempKnot <- length(aR) + 1
      improve <- 1
      improveEx <- 1
      newRadii = tempRadii
    }
  }
  # cat('Current Fit out: ', fitStat, '\n')
  if (improve){
    ####track<-rbind(track,cbind("move",t(newR),fitStat,tempaRSQ,tempGCV))
    list(fitStat=fitStat,newR=newR,tempKnot=tempKnot,improve=improve,improveEx=improveEx,track=track, out.lm=out.lm,radiusIndices=newRadii,models=models)
  }else{list(improve=0,improveEx=0,models=models)}
}
