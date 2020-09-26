"drop.step_2d_badfit" <- function(radii,invInd,dists,explData,response,knotgrid,maxIterations,fitnessMeasure,
                           point,knotPoint,position,aR,BIC,track,out.lm,improveDrop,minKnots,tol=0,baseModel,radiusIndices,models, interactionTerm, data, initDisp, cv.opts, basis,hdetest) {
  
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
        if (isS4(out.lm)) {
          n_predictors <- dim(out.lm@predictors)[2]
          hdeff_vals <- matrix(hdeff(out.lm), ncol=2, byrow=TRUE)
          if (nrow(hdeff_vals) > length(aR)) {
            hdeff_vals <- hdeff_vals[2:nrow(hdeff_vals),]
          }
          badknots <-  seq(length(aR))[apply(hdeff_vals, 1, any)]
          i <- sample(length(badknots), size=1)
        } else {
          badknots<-data.frame(knots=aR, abscoeffs = abs(coef(out.lm))[-1], ses=sqrt(diag(summary(out.lm)$cov.robust))[-1])
          badknots$dif<-badknots$abscoeffs - badknots$ses
          i <- which(badknots$dif==min(badknots$dif))
        }
        tempR <- aR
        tempR <- tempR[-i]
        tempRadii = radiusIndices[-i]
        output<-fit.thinPlate_2d(fitnessMeasure, dists,tempR,radii,baseModel,tempRadii,models, fitStat, interactionTerm, data, initDisp, cv.opts, basis, hdetest)
        initModel<-output$currentModel
        models<-output$models
        initBIC<-output$fitStat
        #get.measure_2d(fitnessMeasure,BIC,initModel, data,  dists, tempR,radii, tempRadii, initDisp)$fitStat
        out<-choose.radii(initBIC,1:length(tempRadii),tempRadii,radii,initModel,dists,tempR,baseModel,fitnessMeasure,response,models, interactionTerm, data, initDisp, cv.opts, basis, hdetest)
        tempRadii=out$radiusIndices
        tempOut.lm=out$out.lm
        models=out$models
        #output<-get.measure_2d(fitnessMeasure,fitStat,tempOut.lm, data,  dists, tempR,radii,tempRadii, initDisp)
        #fitStat<-output$tempMeasure
        tempMeasure<-out$BIC
        
          out.lm <- tempOut.lm
          fitStat<-tempMeasure
          print("bad knot removed **********************************")
          ####print(length(as.vector(coefficients(out.lm))))
          ####print(tempR)
          #print(fitStat)
          newR <- tempR
          newRadii = tempRadii
          tempKnot <- i
        if (isS4(tempOut.lm)){
          n_predictors <- dim(tempOut.lm@predictors)[2]
          hdeff_vals <- matrix(hdeff(tempOut.lm), ncol=2, byrow=TRUE)
          if (nrow(hdeff_vals) > length(tempR)) {
            hdeff_vals <- hdeff_vals[2:nrow(hdeff_vals),]
          }
          badfit_test <- ifelse(hdetest, sum(apply(hdeff_vals, 1, any)), FALSE)
        } else {
          badfit_test <- summary(tempOut.lm)$dispersion > initDisp
        }
        if (badfit_test) {  
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
