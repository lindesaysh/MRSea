"drop.step_2d_badfit" <- function(radii,invInd,dists,explData,response,knotgrid,maxIterations,fitnessMeasure,
                           point,knotPoint,position,aR,BIC,track,out.lm,improveDrop,minKnots,tol=0,baseModel,radiusIndices,models, 
						   interactionTerm, data, initDisp, cv.opts, basis) {

  
  if (isS4(baseModel)){
    attributes(baseModel@misc$formula)$.Environment<-environment()
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
  }
  
  print("******************************************************************************")
  print("Finding suitable initialise...")
  print("******************************************************************************")
  # cat('Current Fit in: ', BIC, '\n')
  badfit <- 1
  fuse<-1
  improvebadDrop<-0
  fitStat<-BIC[length(BIC)]
  newRadii = radiusIndices
  
  while (badfit) {
    print(paste("remove bad knot", fuse))
    fuse <- fuse + 1
    badfit <- 0
    print(fitStat)
  
      if (length(aR) > minKnots) {
        
        twoDcoeffid <- grep("LRF.", names(coefficients(out.lm)))
        length(twoDcoeffid)/length(aR)
        
        badknots<-data.frame(knots=rep(aR, by=2), abscoeffs = abs(coef(out.lm))[twoDcoeffid], ses=sqrt(diag(summary(out.lm)$cov.robust))[twoDcoeffid])
        
        badknots$dif<-badknots$abscoeffs - badknots$ses
        i <- which(badknots$dif==min(badknots$dif))
        badknotid<-which(aR==badknots[i,1])
        tempR <- aR
        tempR <- tempR[-badknotid]
        tempRadii = radiusIndices[-badknotid]
        output<-fit.thinPlate_2d(fitnessMeasure, dists,tempR,radii,baseModel,tempRadii,models, fitStat, interactionTerm, data, initDisp, cv.opts, basis)
        initModel<-output$currentModel
        models<-output$models
        initBIC<-output$fitStat
        #get.measure_2d(fitnessMeasure,BIC,initModel, data,  dists, tempR,radii, tempRadii, initDisp)$fitStat
        out<-choose.radii(initBIC,1:length(tempRadii),tempRadii,radii,initModel,dists,tempR,baseModel,fitnessMeasure,response,models, interactionTerm, data, initDisp, cv.opts, basis)
        tempRadii=out$radiusIndices
        tempOut.lm=out$out.lm
        models=out$models
        #output<-get.measure_2d(fitnessMeasure,fitStat,tempOut.lm, data,  dists, tempR,radii,tempRadii, initDisp)
        #fitStat<-output$tempMeasure
        tempMeasure<-out$BIC
        print(paste(tempMeasure, fitStat, length(aR), badfit, improvebadDrop))
        
          out.lm <- tempOut.lm
          fitStat<-tempMeasure
          print("bad knot removed **********************************")
          ####print(length(as.vector(coefficients(out.lm))))
          ####print(tempR)
          #print(fitStat)
          newR <- tempR
          newRadii = tempRadii
          tempKnot <- badknotid
        if (summary(tempOut.lm)$dispersion>initDisp) {  
          badfit <- 1
          improvebadDrop <- 0
        }else{
          badfit <- 0
          improvebadDrop <- 1
        }     

        aR <- newR
        point <- c(point,knotPoint[tempKnot])
        position[knotPoint[tempKnot]]<-length(point)
        knotPoint <- knotPoint[-tempKnot]
        radiusIndices <- newRadii
    }
  }
  
  
  ###   }
  list(point=point,knotPoint=knotPoint,position=position,aR=aR,BIC=fitStat,track=track,out.lm=out.lm,improveDrop=improvebadDrop,
       radiusIndices=radiusIndices,models=models)
}
