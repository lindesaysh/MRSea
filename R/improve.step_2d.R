#' Function to move knots to neighbours to see if there is any improvement in fit
#' 
#' @author Cameron Walker, Department of Enginering Science, University of Auckland.
#' 
#' @export
#' 


####################################################################################################################

"improve.step_2d" <- function(gap,knotDist,radii,dists,explData,num,response,knotgrid,maxIterations,fitnessMeasure, point,knotPoint,position,aR,BIC,track,out.lm,improveNudge,tol=0,baseModel,radiusIndices,models, interactionTerm, data, initDisp, cv.opts, basis, printout){
  if (isS4(baseModel)){
    attributes(baseModel@misc$formula)$.Environment<-environment()
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
  }
  if(printout){
    print("Improving...")
    print("******************************************************************************")
  }
    # cat('Current Fit in: ', BIC, '\n')
  improve <- 1
  fuse <- 0
  newRadii = radiusIndices
  while ( (improve) & (fuse < maxIterations) ) {
    fuse <- fuse + 1
    improve <- 0
    fitStat<-BIC[length(BIC)]
    for (i in 1:num) {
      nhbrs<-c()
      otherKnots=aR[-i]
      # find 6 nearest knot points to current knot point
      nhbrs<-order(knotDist[knotPoint[i],])[2:6]
      nhbrs<-na.omit(nhbrs)
      # remove options if already a knot point      
      if(length(na.omit(match(nhbrs, knotPoint)))>0){
        id<-na.omit(match(knotPoint, nhbrs))
        nhbrs<-nhbrs[-id]
      }

      for (j in nhbrs){
        if ((length(otherKnots)==0) || ( min(knotDist[j,otherKnots])>=gap)) {
          tempR<-aR
          tempR[i]<-j
          
          output<-fit.thinPlate_2d(fitnessMeasure, dists,tempR,radii,baseModel,radiusIndices,models, fitStat, interactionTerm, data, initDisp, cv.opts, basis)
          initModel<-output$currentModel
          models<-output$models
          initBIC<-output$fitStat
          tempRadii=radiusIndices
          tempOut.lm=output$currentModel
          tempMeasure<-initBIC
          
          if (tempMeasure + tol < fitStat) {
            out.lm <- tempOut.lm
            fitStat<-tempMeasure
            if(printout){
              print("improve *************************************")
            }
            
            newR <- tempR
            newRadii = tempRadii
            tempKnot <- i
            adjNode <- j
            improve <- 1
            improveNudge <- 1
          }      
        }
      }    
    }
    if (improve) {
      point[position[adjNode]] <- knotPoint[tempKnot]
      position[knotPoint[tempKnot]] <- position[adjNode]
      position[adjNode] <- 0
      knotPoint[tempKnot] <- adjNode
      aR<-newR
      BIC <- fitStat
    }                
  }
  list(point=point,knotPoint=knotPoint,position=position,aR=aR,BIC=BIC,track=track,out.lm=out.lm,improveNudge=improveNudge,
       radiusIndices=newRadii,models=models)
}