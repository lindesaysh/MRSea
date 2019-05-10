#' Function for exchanging knot locations and re-fitting model to find best one
#'
#'
#' @author Cameron Walker, Department of Enginering Science, University of Auckland.
#'
#' @export
#'

################################################################################################################
"exchange.step_2d" <- function(gap,knotDist,radii,dists,explData,response,knotgrid,maxIterations,fitnessMeasure, point,knotPoint,position,aR,BIC,track,out.lm,improveEx, maxKnots,tol=0,baseModel,radiusIndices,models, interactionTerm, data, initDisp, cv.opts, basis){

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
    if (isS4(baseModel)){
      indexdat<-order(rowSums(abs(residuals(baseModel, type='pearson'))), decreasing = TRUE)[1:5]
    } else {
      indexdat<-order(abs(residuals(baseModel, type='pearson')), decreasing = TRUE)[1:5]
    }
    #### Find available knots
    legPos<-position[which(apply(knotDist[point,aR],1,min)>=gap)]
    index<-c()
     for(i in 1:5){
      new<-scale(knotgrid[point,],center=c(explData[indexdat[i],1],explData[indexdat[i],2]))
      index[i]<-which.min(abs(new[,1])+abs(new[,2]))
    }
   index<-unique(index)

    if (length(legPos)>0) {
      #index<-which.min(abs(new[,1])+abs(new[,2]))
      if (!(any(knotDist[point[index[1]],aR]<gap))) {
        output <- move.knot_2D(radii,dists,explData,index,fitnessMeasure,BIC,aR,point,
                               response,knotgrid,out.lm,improve,improveEx, track,
                               maxKnots,tol,baseModel,radiusIndices,models, interactionTerm, data, initDisp, cv.opts, basis)

        improve <- output$improve
        improveEx <- output$improveEx
        models <- thinModels(output$models)
        if (1 - improve) break
        out.lm<-output$out.lm
        track <- output$track
        tempKnot <- output$tempKnot
        if (tempKnot <= length(knotPoint)) {
          position[knotPoint[tempKnot]] <- output$index
          position[point[output$index]] <- 0
          buff <- knotPoint[tempKnot]
          knotPoint[tempKnot] <- point[output$index]
          point[output$index] <- buff
        } else {
          knotPoint<-c(knotPoint,point[output$index])
          position[point[output$index]]<-0
          point<-point[-output$index]
          if (length(point) >= output$index) {
            for (i in output$index:length(point)){
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
