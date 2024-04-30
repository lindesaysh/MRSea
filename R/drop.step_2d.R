#' Function that tries dropping knots to find an improvement in fit
#' 
#' @author Cameron Walker, Department of Engineering Science, University of Auckland.
#' 
#' @export
#' 



"drop.step_2d" <- function(radii,invInd,dists,explData,response,knotgrid,maxIterations,fitnessMeasure,
                           point,knotPoint,position,aR,BIC,track,out.lm,improveDrop,minKnots,tol=0,baseModel,radiusIndices,models, interactionTerm, data, initDisp, cv.opts, basis, printout) {
  
  if (isS4(baseModel)){
    attributes(baseModel@misc$formula)$.Environment<-environment()
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
  }
  if(printout){
    print("******************************************************************************")
    print("Simplifying model...")
    print("******************************************************************************")
  }
  fuse <- 0
  improve <- 1
  fitStat<-BIC[length(BIC)]
  newRadii = radiusIndices
  while ((improve) & (fuse < maxIterations) ) {

    fuse <- fuse + 1
    improve <- 0

    if (length(aR) > minKnots) {
      for (i in 1:length(aR)) {
        tempR <- aR
        tempR <- tempR[-i]
        tempRadii = radiusIndices[-i]
        output<-fit.thinPlate_2d(fitnessMeasure, dists,tempR,radii,baseModel,tempRadii,models, fitStat, interactionTerm, data, initDisp, cv.opts, basis)
        initModel<-output$currentModel
        models<-output$models
        initBIC<-output$fitStat
        out<-choose.radii(initBIC,1:length(tempRadii),tempRadii,radii,initModel,dists,tempR,baseModel,fitnessMeasure,response,models, interactionTerm, data, initDisp, cv.opts, basis, printout)
        tempRadii=out$radiusIndices
        tempOut.lm=out$out.lm
        models=out$models

        tempMeasure<-out$BIC

        if (tempMeasure +tol < fitStat) {
          out.lm <- tempOut.lm
          fitStat<-tempMeasure
          if(printout){
            print("drop **********************************")
          }
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
        radiusIndices <- newRadii
      }
    }
  }
  
  
  ###   }
  list(point=point,knotPoint=knotPoint,position=position,aR=aR,BIC=fitStat,track=track,out.lm=out.lm,improveDrop=improveDrop,
       radiusIndices=radiusIndices,models=models)
}

