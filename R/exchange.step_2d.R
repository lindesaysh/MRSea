#' Function for exchanging knot locations and re-fitting model to find best one
#'
#'
#' @author Cameron Walker, Department of Enginering Science, University of Auckland.
#'
#' @export
#'

################################################################################################################
"exchange.step_2d" <- function(gap,knotDist,radii,invInd,dists,explData,response,knotgrid,maxIterations,fitnessMeasure, point,knotPoint,position,aR,BIC,track,out.lm,improveEx, maxKnots,tol=0,baseModel,radiusIndices,models, interactionTerm, data, initDisp){

  #attributes(baseModel$formula)$.Environment<-environment()

  # Loop - fuse used to ensure algorithm terminates
  print("******************************************************************************")
  print("Exchanging...")
  print("******************************************************************************")
  # cat('Current Fit in: ', BIC, '\n')
  fuse <- 0
  improve <- 1
  while ( (improve) & (fuse < maxIterations) ) {
    fuse <- fuse + 1
    improve <- 0
    index<-max.col(t(abs(resid(out.lm,type="pearson"))))
    ####fix for grid approach#####
    new<-scale(knotgrid[point,],center=c(explData[index,1],explData[index,2]))
    ####Pick nearest grid point that is also far enough away from another knot
    legPos<-position[which(apply(knotDist[point,aR],1,min)>=gap)]
    if (length(legPos)>0) {
      index<-which.min(abs(new[,1])+abs(new[,2]))
      if (!(any(knotDist[point[index],aR]<gap))) {
        output <- move.knot_2D(radii,invInd,dists,explData,index,fitnessMeasure,BIC,aR,point,
                               response,knotgrid,out.lm,improve,improveEx, track,
                               maxKnots,tol,baseModel,radiusIndices,models, interactionTerm, data, initDisp)
        improve <- output$improve
        improveEx <- output$improveEx
        models <- thinModels(output$models)
        if (1 - improve) break
        out.lm<-output$out.lm
        track <- output$track
        tempKnot <- output$tempKnot
        if (tempKnot <= length(knotPoint)) {
          position[knotPoint[tempKnot]] <- index
          position[point[index]] <- 0
          buff <- knotPoint[tempKnot]
          knotPoint[tempKnot] <- point[index]
          point[index] <- buff
        } else {
          knotPoint<-c(knotPoint,point[index])
          position[point[index]]<-0
          point<-point[-index]
          if (length(point) >= index) {
            for (i in index:length(point)){
              position[point[i]]<-position[point[i]]-1
            }
          }
        }
        aR <- output$newR
        BIC <- output$fitStat
        radiusIndices <- output$radiusIndices
        ####track<- rbind(track,cbind("exchange",t(aR),BIC[length(BIC)],adjRsq[length(adjRsq)],GCV[length(GCV)]))
      }
    }else {
      print("no legal knot positions available")
    }
  }
  # cat('Current Fit out: ', BIC, '\n')
  list(point=point,knotPoint=knotPoint,position=position,aR=aR,BIC=BIC,track=track,out.lm=out.lm,
       improveEx=improveEx,radiusIndices=radiusIndices,models=models)
}

###############################
