#' Function to move knots to neighbours to see if there is any improvement in fit
#' 
#' @author Cameron Walker, Department of Enginering Science, University of Auckland.
#' 
#' @export
#' 


####################################################################################################################

"improve.step_2d" <- function(gap,knotDist,radii,invInd,dists,explData,num,response,knotgrid,maxIterations,fitnessMeasure, point,knotPoint,position,aR,BIC,track,out.lm,improveNudge,tol=0,baseModel,radiusIndices,models, interactionTerm, data, initDisp){
  attributes(baseModel$formula)$.Environment<-environment()
  print("Improving...")
  print("******************************************************************************")
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
#       if (grid[knotPoint[i],1] > 1) {
#         if (!(is.na(gridResp[knotPoint[i] - 1]))) {
#           nhbrs<- c(nhbrs, knotPoint[i] - 1)
#         }
#         if ((grid[knotPoint[i],2] > 1) && !(is.na(gridResp[knotPoint[i] - 1 - xVals]))) {
#           nhbrs<-c(nhbrs,knotPoint[i] - 1 - xVals)
#         }
#         if ((grid[knotPoint[i],2] < yVals) && !(is.na(gridResp[knotPoint[i] - 1 + xVals]))) {
#           nhbrs<-c(nhbrs,knotPoint[i] - 1 + xVals)
#         }
#       }
#       if (grid[knotPoint[i],1] < xVals) {
#         if (!(is.na(gridResp[knotPoint[i] + 1]))) {
#           nhbrs<- c(nhbrs, knotPoint[i] + 1)
#         }
#         if ((grid[knotPoint[i],2] > 1) && !(is.na(gridResp[knotPoint[i] + 1 - xVals]))) {
#           nhbrs<-c(nhbrs,knotPoint[i] + 1 - xVals)
#         }
#         if ((grid[knotPoint[i],2] < yVals) && !(is.na(gridResp[knotPoint[i] + 1 + xVals]))) {
#           nhbrs<-c(nhbrs,knotPoint[i] + 1 + xVals)
#         }
#       }
#       if ((grid[knotPoint[i],2] > 1) && !(is.na(gridResp[knotPoint[i] - xVals]))) {
#         nhbrs<- c(nhbrs, knotPoint[i] - xVals)
#       }
#       if ((grid[knotPoint[i],2] < yVals) && !(is.na(gridResp[knotPoint[i] + xVals]))) {
#         nhbrs<- c(nhbrs, knotPoint[i] + xVals)
#       }

      # find 6 nearest knot points to current knot point
      nhbrs<-point[order(knotDist[i,])[2:6]]
      nhbrs<-na.omit(nhbrs)
      # remove options if already a knot point      
      if(length(na.omit(match(nhbrs, knotPoint)))>0){
        id<-na.omit(match(nhbrs, knotPoint))
        nhbrs<-nhbrs[-id]
      }

      for (j in nhbrs){
#         if (is.na(position[j])) {
#           print(j)
#           print(knotPoint[i])
#           print(xVals)
#           print(nhbrs)
#         }
        #cat(invInd[otherKnots], '\n', invInd[j], '\n', aR, '\n')
        #cat('parta: ', length(otherKnots), '\n')
        #cat('partb: ', min(knotDist[invInd[j],invInd[otherKnots]]), '\n')
#         cat('other ', otherKnots, '\n')
#         cat('invind other ', invInd[otherKnots], '\n')
#         cat('j ', j, '\n')
#         cat('invindj ', invInd[j], '\n')
        
        if ((length(otherKnots)==0) || ( min(knotDist[j,otherKnots])>=gap)) {
          tempR<-aR
          tempR[i]<-j
          
          output<-fit.thinPlate_2d(fitnessMeasure, dists,tempR,radii,baseModel,radiusIndices,models, fitStat, interactionTerm, data, initDisp)
          initModel<-output$currentModel
          models<-output$models
          initBIC<-get.measure_2d(fitnessMeasure,fitStat,initModel, data,  dists, tempR,radii,radiusIndices, initDisp)$fitStat
          
          
          out<-choose.radii(initBIC,i,radiusIndices,radii,initModel,dists,tempR,baseModel,fitnessMeasure,response,models, interactionTerm, data, initDisp)
          tempRadii=out$radiusIndices
          tempOut.lm=out$out.lm
          models=out$models
          output<-get.measure_2d(fitnessMeasure,fitStat,tempOut.lm, data,  dists, tempR,radii,tempRadii, initDisp)
          
          #fitStat<-output$tempMeasure
          tempMeasure<-output$fitStat
          #### print(length(as.vector(coefficients(tempOut.lm))))
          if (tempMeasure + tol < fitStat) {
            out.lm <- tempOut.lm
            fitStat<-tempMeasure
            print("improve *************************************")
            #print(fitStat)
            newR <- tempR
            newRadii = tempRadii
            tempKnot <- i
            adjNode <- j
            improve <- 1
            improveNudge <- 1
          }      
        }
        ##}
        # cat(improve, '\n')
      }    
    }
    if (improve) {
      #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      point[position[adjNode]] <- knotPoint[tempKnot]
      position[knotPoint[tempKnot]] <- position[adjNode]
      position[adjNode] <- 0
      knotPoint[tempKnot] <- adjNode
      #tempR<-aR
      #tempR[tempKnot] <- adjNode
      #aR <- tempR
      aR<-newR
      BIC <- fitStat
      ####track<-rbind(track,cbind("improve",t(tempR),BIC[length(BIC)],adjRsq[length(adjRsq)],GCV[length(GCV)]))
    }                
  }
  list(point=point,knotPoint=knotPoint,position=position,aR=aR,BIC=BIC,track=track,out.lm=out.lm,improveNudge=improveNudge,
       radiusIndices=newRadii,models=models)
}