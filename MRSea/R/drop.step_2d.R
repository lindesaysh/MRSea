#' Function that tries dropping knots to find an improvement in fit
#' 
#' @author Cameron Walker, Department of Enginering Science, University of Auckland.
#' 
#' @export
#' 



"drop.step_2d" <- function(radii,invInd,dists,explData,response,explanatory,maxIterations,fitnessMeasure,
                           point,knotPoint,position,aR,BIC,track,out.lm,improveDrop,minKnots,tol=0,baseModel,radiusIndices,models, interactionTerm, data, initDisp) {
  
  attributes(baseModel$formula)$.Environment<-environment()
  
  print("******************************************************************************")
  print("Simplifying model...")
  print("******************************************************************************")
  # cat('Current Fit in: ', BIC, '\n')
  improve<-0
  fitStat<-BIC[length(BIC)]
  newRadii = radiusIndices
  for (i in 1:length(aR)) {
    tempR <- aR
    tempR <- tempR[-i]
    tempRadii = radiusIndices[-i]
    output<-fit.thinPlate_2d(fitnessMeasure, dists,invInd[tempR],radii,baseModel,tempRadii,models, fitStat, interactionTerm, data, initDisp)
    initModel<-output$currentModel
    models<-output$models
    initBIC<-get.measure_2d(fitnessMeasure,BIC,initModel, data,  dists, invInd[tempR],radii, tempRadii, initDisp)$fitStat
    out<-choose.radii(initBIC,1:length(tempRadii),tempRadii,radii,initModel,dists,invInd[tempR],baseModel,fitnessMeasure,response,models, interactionTerm, data, initDisp)
    tempRadii=out$radiusIndices
    tempOut.lm=out$out.lm
    models=out$models
    output<-get.measure_2d(fitnessMeasure,fitStat,tempOut.lm, data,  dists, invInd[tempR],radii,tempRadii, initDisp)
    
    #fitStat<-output$tempMeasure
    tempMeasure<-output$fitStat
    if (tempMeasure +tol < fitStat) {
      out.lm <- tempOut.lm
      fitStat<-tempMeasure
      print("drop **********************************")
      ####print(length(as.vector(coefficients(out.lm))))
      ####print(tempR)
      #print(fitStat)
      newR <- tempR
      newRadii = tempRadii
      tempKnot <- i
      improve <- 1
      improveDrop <- 1
    }      
  }
  
  if (improve) {
    aR <- newR
    point <- c(point,knotPoint[tempKnot])
    position[knotPoint[tempKnot]]<-length(point)
    knotPoint <- knotPoint[-tempKnot]
  }
  ###   }
  list(point=point,knotPoint=knotPoint,position=position,aR=aR,BIC=fitStat,track=track,out.lm=out.lm,improveDrop=improveDrop,
       radiusIndices=newRadii,models=models)
}

