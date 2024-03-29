#' Function for exchanging knot locations and re-fitting model to find best one
#'
#'
#' @author Cameron Walker, Department of Enginering Science, University of Auckland and Lindesay Scott-Hayward (University of St Andrews)
#'
#' @export
#'

################################################################################################################
"exchange.step_2d" <- function(gap,knotDist,radii,dists,explData,response,knotgrid,maxIterations,fitnessMeasure, point,knotPoint,position,aR,BIC,track,out.lm,improveEx, maxKnots,tol=0,baseModel,radiusIndices,models, interactionTerm, data, initDisp, cv.opts, basis, printout){

  #attributes(baseModel$formula)$.Environment<-environment()

  # Loop - fuse used to ensure algorithm terminates
  if(printout){
    print("******************************************************************************")
    print("Exchanging...")
    print("******************************************************************************")
  }
  fuse <- 0
  improve <- 1

  while ( (improve) & (fuse < maxIterations) ) {
    fuse <- fuse + 1
    improve <- 0
    if (isS4(baseModel)){
      indexdat<-order(rowSums(abs(residuals(out.lm, type='pearson'))), decreasing = TRUE)[1:5]
    } else {
      indexdat<-order(abs(residuals(out.lm, type='pearson')), decreasing = TRUE)[1:5]
    }
    #### Find available knots
    legPos<-position[which(apply(knotDist[point,aR],1,min)>=gap)]
    if(ncol(knotgrid)>2){
      nm<-names(knotgrid)[3]
      residchunk<-eval(parse(text=paste('data$', nm, '[indexdat[1]]')))
      knotchunkid<-which(knotgrid[point,nm]==residchunk)
      new<-scale(knotgrid[point[knotchunkid],1:2],center=c(explData[indexdat[1],1],explData[indexdat[1],2]))
      # which aR are in residchunk
      aRresidchunk<-aR[which(knotgrid[aR,nm]==residchunk)]
      legPos1<-position[which(knotgrid[,nm]==residchunk)]
      legPos<-legPos1
      index<-knotchunkid[which.min(abs(new[,1])+abs(new[,2]))]
    }else{
      new<-scale(knotgrid[point,1:2],center=c(explData[indexdat[1],1],explData[indexdat[1],2]))
      index<-which.min(abs(new[,1])+abs(new[,2]))
    }
    # ggplot() + geom_point(data=knotgrid, aes(X1, X2)) + facet_wrap(~yearmonth) +
    #   geom_point(data=knotgrid[point[knotchunkid],], aes(X1, X2), shape=2, size=2)+
    #   geom_point(data=knotgrid[aR,], aes(X1, X2), shape=3, size=2, col='blue') +
    #   geom_point(data=knotgrid[aRresidchunk,], aes(X1, X2), shape=3, size=4, col='blue') +
    #   geom_point(data=knotgrid[point[legPos1],], aes(X1, X2), shape=3, size=2, col='red')+
    #   geom_point(data=knotgrid[point[legPos1[legPos2]],], aes(X1, X2), shape=4, size=2, col='green') +
    #   geom_point(data=knotgrid[point[legPos],], aes(X1, X2), shape=5, size=3, col='maroon') + coord_equal() + geom_point(data=data[indexdat[1],], aes(x.pos, y.pos), size=3, col='thistle') + geom_point(data=knotgrid[point[index],], aes(X1, X2), shape=7, size=3, col='black')
    # 
    
    if (length(legPos)>0) {
      #index<-which.min(abs(new[,1])+abs(new[,2]))
      if (!(any(knotDist[point[index[1]],aR]<gap))) {
        output <- move.knot_2D(radii,dists,explData,index,fitnessMeasure,BIC,aR,point,
                               response,knotgrid,out.lm,improve,improveEx, track,
                               maxKnots,tol,baseModel,radiusIndices,models, interactionTerm, data, initDisp, cv.opts, basis, printout)

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
      if(printout){
        print("no legal knot positions available")
      }
    }
  }
  # cat('Current Fit out: ', BIC, '\n')
  list(point=point,knotPoint=knotPoint,position=position,aR=aR,BIC=BIC,track=track,out.lm=out.lm,
       improveEx=improveEx,radiusIndices=radiusIndices,models=models)
}

###############################
