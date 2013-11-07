#' Code for adaptively spacing knots for a spatial smooth. The smoothing process uses a CReSS basis.
#' 
#' @param splineParams List object where the first element [[1]] contains a list of objects for the 2D SALSA fitting process: \code{knotDist}, \code{radii}, \code{dist}, \code{gridresp}, \code{grid}, \code{datacoords}, \code{response}, \code{knotgrid}, \code{minKnots}, \code{maxKnots}, \code{gap}
#' @param startKnots number of space-filled knots to start with (between minKnots and maxKnots)
#' @param fitnessMeasure (default=BIC). Measure used to evaluate the fit. Other options are AIC, AICc, BIC, QICb (Quasi-Likelihood Information Criterion with log(n) penalty)
#' @param maxIterations exchange/improve heuristic will terminate after \code{maxIterations} if still running
#' @param baseModel starting model for SALSA to use.  Must not already contain a spatial smooth.
#' @param radiusIndices vector of length \code{startKnots} identifying which radii (\code{splineParams[[1]]$radii}) will be used to initialise the model
#' @param initialise (default = TRUE). Logical stating whether or not to start with space-filled knots (TRUE) or user specified locations (FALSE)
#' @param initialKnots If \code{initialise=FALSE} then the start locations for the knots are specified in \code{initialKnots}. Must be coordinates.
#' @param interactionTerm (default=NULL). Specifies which term in \code{baseModel} the spatial smooth will interact with.  If NULL no interaction term is fitted
#' @param winHalfWidth Half-width of window used to calculate region with biggest average residual magnitude
#' @param tol Tolerance for difference between fit measures.  E.g. \code{tol=2} means that the calculated fitness measures must be 2 units apart to be considered different
#' 
#' 
#' @details The following are the details of the splineParams[[1]] objects.  Note.  If salsa1D has been run then details for those covariates will sit in splineParams[[2]] and onward.
#' 
#' \code{knotDist}: matrix of knot to knot distances (k x k).  May be Euclidean or geodesic distances. Must be square and the same dimensions as \code{nrows(na.omit(knotgrid))}
#' 
#' \code{radii} Sequence of range parameters for the CReSS basis from local (small) to global (large).  Determines the range of the influence of each knot.
#' 
#' \code{dist}: matrix of distances between data locations and knot locations (n x k). May be Euclidean or geodesic distances.
#' 
#' \code{gridresp} The first column of knotgrid
#' 
#' \code{grid} Index of knotgrid locations.  Should be same length as \code{knotgrid} but with x=integer values from 1 to number of unique x-locations and y= integer values from 1 to number of unique y-locations.
#' 
#' \code{datacoords}: Coordinates of the data locations
#' 
#' \code{response}: vector of response data for the modelling process
#' 
#' \code{knotgrid}: grid of legal knot locations.  Must be a regular grid with c(NA, NA) for rows with an illegal knot
#' 
#' \code{minKnots}: minimum number of knots to fit
#' 
#' \code{maxKnots}: maximum number of knots to fit
#' 
#' \code{gap}: Minimum gap between knots (in unit of measurement of \code{datacoords})
#' 

"return.reg.spline.fit.2d" <- function(splineParams, startKnots, winHalfWidth,fitnessMeasure="BIC", maxIterations=100, tol=0, baseModel=NULL, radiusIndices=NULL, initialise=TRUE, initialKnots=NULL, interactionTerm=NULL, knot.seed=10){

#Where am I?
# requires splines library and mgcv ibrary to be loaded!!

# PARAMETERS
# minKnots:        minimum number of knots to fit
# maxKnots:        maximum number of knots to fit
# gap:             minimum gap between knots (i.e. number of data points)
# winHalfWidth:    half-width of window used to calculate region with biggest average residual magnitude
# fitnessMeasure:  measure used to evaluate the fit (value = 1--4)
# fitnessMeasure=="AIC" uses AIC       
# fitnessMeasure=="AICc" uses AICc       
# fitnessMeasure=="BIC" uses BIC       
# fitnessMeasure=="QAIC" uses QAIC       
# fitnessMeasure=="QAICc" uses QAICc       
# fitnessMeasure=="CV.offset" uses CV where an offset is allowed
# fitnessMeasure== "CV" uses CV where blocking of the data in the folds is allowed.  offset also allowed
  
# maxIterations:   exchange/improve heuristic will terminate after maxIterations if still running

# output:          aR            - list of knot points for each iteration
#                  arSq          - list of adjusted r-squareds for each it.
#                  finaR         - list of optimal knots for each number of knots
#                  finarSq       - list of adjusted r-squareds for each number of knots
#                  finBIC        - list of BICs for each number of knots

# pointers:        knotPoint     - the index of the knot points (i, where explanatory[i] is a knot)
#                  point         - the index of the other points
#                  position      - the index in point of ith data point
#                                  0 otherwise (position[i] = j,where point[j] = i)

# x and y should be indices here - x = 1, 2, ...xVals, 1, 2, ...
#                                  y = 1, 1,...1, 2, 2,...,yVals, yVals,...

# data:           grid           - row and column index of possible knot locations
#                 gridResp       - x coordinate of possible knot locations
#                 gridData       - x and y coordinate of possible knot locations
#                 explData       - x and y coordinates of observations
#                 response       - values of observations
#                 explanatory    - same as gridData?
#                 radii          - possible values of radius
#                 radiusIndices  - index in radii for each knot radius (e.g. [2,3,1,2,1]
#                 interactionTerm - allows interaction between space and another term "ns(Year, knots =   
                  # splineParams[[6]][[2]])", for example.  If this is NULL, no interaction is done
 
  # split out spline parameter object into its pieces
  knotDist <- splineParams[[1]]$knotDist
  radii <- splineParams[[1]]$radii
  dists <- splineParams[[1]]$dist
  gridResp <- splineParams[[1]]$gridResp
  grid <- splineParams[[1]]$grid
  explData <- splineParams[[1]]$datacoords
  response <- splineParams[[1]]$response
  explanatory <- splineParams[[1]]$knotgrid
  minKnots <- splineParams[[1]]$minKnots
  maxKnots <- splineParams[[1]]$maxKnots
  gap <- splineParams[[1]]$gap
  
  
  
  attributes(baseModel$formula)$.Environment<-environment()
  data<- baseModel$data
  baseModel<-update(baseModel, data=data)
  
  ###########################triangulation of points############################
  x<-as.vector(grid[,1])
      y<-as.vector(grid[,2])
      #Get dimensions of grid
      xvals <- max(x)
      yvals <- max(y)
      ###########################initialisation######################################
      output <- initialise.measures_2d(knotDist,maxIterations,gap,radii,dists,gridResp,explData,startKnots,xvals, yvals, explanatory, response, baseModel, radiusIndices, initialise, initialKnots,fitnessMeasure, interactionTerm, data, knot.seed)

      
      
      point <- output$point
      knotPoint <- output$knotPoint
      position <- output$position
      aR <- output$aR
      BIC <- output$BIC
      track <- output$track
      out.lm <- output$out.lm
      invInd <- output$invInd
      models <- (output$models)
      radiusIndices <- output$radiusIndices

    ####################################algorithm loop#############################
    improveEx <- 1
    improveNudge <- 1
    improveDrop <- 1
    overallImprove = 0
    while (improveEx | improveNudge | improveDrop) {
      improveEx <- 0
      improveNudge <- 0
      improveDrop <- 0
    ####################################exchange step#############################
     ####track <- rbind(track,cbind("exchanging",t(aR),BIC[length(BIC)],adjRsq[length(adjRsq)],GCV[length(GCV)]))
      output <- exchange.step_2d(gap,knotDist,radii,invInd,dists,explData,response,explanatory,maxIterations,fitnessMeasure, point,knotPoint,position,aR,BIC,track,out.lm,improveEx,maxKnots,tol,baseModel,radiusIndices,models, interactionTerm, data)
    #  ####print("here e")
      point <- output$point
      knotPoint <- output$knotPoint
      position <- output$position
      aR <- output$aR
      BIC <- output$BIC
      track <- output$track
      models <- output$models
      out.lm <- output$out.lm
      radiusIndices <- output$radiusIndices
      improveEx <- output$improveEx
    ######################################improve step############################
      ####track <- rbind(track,cbind("improving",t(aR),BIC[length(BIC)],adjRsq[length(adjRsq)],GCV[length(GCV)]))
      ####print("here im")
      output <- improve.step_2d(gap,knotDist,radii,invInd,dists,gridResp,grid,explData,xvals, yvals, length(aR),response,explanatory,maxIterations,fitnessMeasure, point,knotPoint,position,aR,BIC,track,out.lm,improveNudge,tol,baseModel,radiusIndices,models, interactionTerm, data)
      ####print("here im")
      point <- output$point
      knotPoint <- output$knotPoint
      position <- output$position
      aR <- output$aR
      BIC <- output$BIC
      track <- output$track
      models <- thinModels(output$models)
      out.lm <- output$out.lm
      radiusIndices <- output$radiusIndices
      improveNudge <- output$improveNudge
    ###################################drop step#################################
      if (length(aR) > minKnots) {
         output <- drop.step_2d(radii,invInd,dists,explData,response,explanatory,maxIterations,fitnessMeasure,point,knotPoint,position,aR,BIC,track,out.lm,improveDrop,minKnots,tol,baseModel,radiusIndices,models, interactionTerm, data)
         ####print("here e")
         point <- output$point
         knotPoint <- output$knotPoint
         position <- output$position
         aR <- output$aR
         BIC <- output$BIC
         ####track <- output$track
         models <- thinModels(output$models)
         out.lm <- output$out.lm
         radiusIndices <- output$radiusIndices
         improveDrop <- output$improveDrop      
         }
       if ((improveEx) | (improveNudge) | (improveDrop)) overallImprove = 1
      }
   ####################################write to file###############################
    ####track <- rbind(track,cbind("writing",t(aR),BIC[length(BIC)],adjRsq[length(adjRsq)],GCV[length(GCV)]))
  ####print("here fin")
  print("And we're done...")
  return(list(outputFS=c(length(aR),BIC[length(BIC)],aR),aR=aR,track=track, radiusIndices = radiusIndices,
            out.lm=out.lm,invInd=invInd,models=models,actualKnotIndices=invInd[aR],improve=overallImprove))

}


######################################################################################################################
"initialise.measures_2d" <- function(knotDist,maxIterations,gap,radii,dists,gridResp,explData,startKnots,xvals, yvals, explanatory, response, baseModel,radiusIndices, initialise, initialKnots,fitnessMeasure, interactionTerm, data, knot.seed){
  
  attributes(baseModel$formula)$.Environment<-environment()
  baseModel<-update(baseModel, data=data)
  
  print("******************************************************************************")
  print("Initialising...")
  print("******************************************************************************")
  #greedy pick from the 1st grid point
  fuse=0

  mapKnotPoint <- c()
  legPos <- c()
  numRand <- 0
    
  track <- cbind()
  # recall: gridResp has x coordinates of possible knot locations
  #         or NA if knot location is outside legal region
  oInd<- 1:length(gridResp)
  resp<- na.omit(gridResp)
  
  # determine  index of legal knotpoints in gridResp (and hence gridData)
  if (length(gridResp) > length(resp)) {
     mapInd<- oInd[-na.action(resp)]       
  } else {
     mapInd<- oInd
  }
  # make pointer to determine where in mapInd each legal knot is stored - else = 0
  # e.g. if 693 grid points but only 500 legal positions, if the 300th knot poistion is
  # the 200th legal position, then invInd[300] = 200
  invInd <- rep(0,length(gridResp))
  for (i in 1:length(mapInd)) {
    invInd[mapInd[i]] <- i
  }

  #if initialise is TRUE:
  if (initialise) {
    require(fields)
    numNeeded = startKnots
    numGot = 0
   # while ((numGot < numNeeded) && (fuse < maxIterations)) {
   #  fuse = fuse + 1
   #  legPos = mapInd
   #   posKnots = cbind()
    #  for (i in 1:numNeeded) {
    #     newKnot=legPos[floor(runif(1,1,length(legPos)+1))]
    #     posKnots = c(posKnots,newKnot)
    #     numGot=length(posKnots)
    #     if (length(legPos>1)) {
    #       entries=which(apply(as.matrix(knotDist[invInd[legPos],invInd[posKnots]]),1,min)>=gap)
    #       badEntries=which(apply(as.matrix(knotDist[invInd[legPos],invInd[newKnot]]),1,min)<gap)
    #     } else {
    #       if (any(knotDist[invInd[legPos],invInd[posKnots]]<gap)) {
    #         entries=c();badEntries=invInd[legPos]
    #       } else {
    #         badEntries=c();entries=invInd[legPos]
    #       }
    #     }     
    #     fixFault<-na.omit(match(mapInd[invInd[legPos][badEntries]],legPos))
    #     legPos<-legPos[-fixFault]
    #     if (length(legPos)==0) break
    #  }
   #}
   #if (numGot < numNeeded) print("WARNING: less knots fitted than desired")
   #knotPoint<- posKnots
   #print(knotPoint)
   #aR <- knotPoint
   #radiusIndices <-rep((1:length(radii))[ceiling(length(radii)/2)],length(aR))
    print("Space-filling knots....")
    set.seed(knot.seed)
  if(nrow(explData)<1000){initialKnots<- cover.design(explData, nd=numNeeded, nruns=1)$design
                          }else{SampledPoints<- sample(1:dim(explData)[1], min(1000, dim(explData)[1]))
      #space-fill data (subsample - see line above) to get knot locations
      # remove any duplicated points as doesnt work in cover design
       dupPoints <-paste(explData[SampledPoints,1], explData[SampledPoints,2], sep='E')
                                
    if(length(which(duplicated(dupPoints)==T))>0){                                
  initialKnots<- cover.design((explData)[SampledPoints,][-which(duplicated(dupPoints)==T),], nd=numNeeded, nruns=1)$design}
                                
  if(length(which(duplicated(dupPoints)==T))==0){                                
initialKnots<- cover.design((explData)[SampledPoints,], nd=numNeeded, nruns=1)$design}
        }
    
  posKnots = cbind()
  legPos=mapInd
  for (i in 1:(dim(initialKnots)[1])) {
     new<-scale(explanatory[legPos,],center=c(initialKnots[i,1],initialKnots[i,2]))
     ####Pick nearest grid point that is also far enough away from another knot
     ind<-which.min(abs(new[,1])+abs(new[,2]))
     newKnot = legPos[ind]
     posKnots = c(posKnots,newKnot)
     if (length(legPos>1)) {
       entries=which(apply(as.matrix(knotDist[invInd[legPos],invInd[posKnots]]),1,min)>=gap)
       badEntries=which(apply(as.matrix(knotDist[invInd[legPos],invInd[newKnot]]),1,min)<gap)
     } else {
       if (any(knotDist[invInd[legPos],invInd[posKnots]]<gap)) {
         entries=c();badEntries=invInd[legPos]
       } else {
         badEntries=c();entries=invInd[legPos]
       }
     }     
     fixFault<-na.omit(match(mapInd[invInd[legPos][badEntries]],legPos))
     if(length(fixFault)>0){
      legPos<-legPos[-fixFault]
     }
     if (length(legPos)==0) break
   }
    
    
    
   numGot=length(posKnots)
    
    
   if (numGot < dim(initialKnots)[1]) {
     print("WARNING: less knots positioned than desired")
   }
   knotPoint<- posKnots
   # print(c('knots: ',knotPoint))
   aR <- knotPoint
  # print(c('actual knots: ',invInd[aR]))
  radiusIndices <-rep((1:length(radii))[ceiling(length(radii)/2)],length(aR))
    
} else {
  posKnots = cbind()
  legPos=mapInd
  for (i in 1:(dim(initialKnots)[1])) {
     new<-scale(explanatory[legPos,],center=c(initialKnots[i,1],initialKnots[i,2]))
     ####Pick nearest grid point that is also far enough away from another knot
     ind<-which.min(abs(new[,1])+abs(new[,2]))
     newKnot = legPos[ind]
     posKnots = c(posKnots,newKnot)
     if (length(legPos>1)) {
       entries=which(apply(as.matrix(knotDist[invInd[legPos],invInd[posKnots]]),1,min)>=gap)
       badEntries=which(apply(as.matrix(knotDist[invInd[legPos],invInd[newKnot]]),1,min)<gap)
     } else {
       if (any(knotDist[invInd[legPos],invInd[posKnots]]<gap)) {
         entries=c();badEntries=invInd[legPos]
       } else {
         badEntries=c();entries=invInd[legPos]
       }
     }     
     fixFault<-na.omit(match(mapInd[invInd[legPos][badEntries]],legPos))
     if(length(fixFault)>0){
       legPos<-legPos[-fixFault]
     }
     if (length(legPos)==0) break
   }
   numGot=length(posKnots)
   if (numGot < length(initialKnots)) {
     print("WARNING: less knots positioned than desired")
     radiusIndices=radiusIndices[1:numGot]
   }
   knotPoint<- posKnots
   #print(c('knots: ',knotPoint))
   aR <- knotPoint
    # print(c('actual knots: ',invInd[aR]))
  } # end of ifelse statement related to initialise  = T/F
  
  
  print("Initialising model...")
  models = vector("list",0)
  
   if(fitnessMeasure=="AIC"){       
           fitStat <- AIC(baseModel)}
  
   if(fitnessMeasure=="AICc"){       
           fitStat <- AICc(baseModel)}
  
   if(fitnessMeasure=="BIC"){       
           fitStat <- BIC(baseModel)}
  
   if(fitnessMeasure=="QAIC"){       
           if(baseModel$family[1]=="quasipoisson"){
     PoisMod<-update(baseModel, round(.)~., family=poisson)
           fitStat <- QAIC(PoisMod, chat = summary(baseModel)$dispersion)}
         if(baseModel$family[1]=="quasibinomial"){
     BinMod<-update(baseModel, family=binomial)
           fitStat <- QAIC(BinMod, round(.)~., chat = summary(baseModel)$dispersion)}
   }
  
   if(fitnessMeasure=="QAICc"){       
           if(baseModel$family[1]=="quasipoisson"){
     PoisMod<-update(baseModel, family=poisson)
           fitStat <- QAICc(PoisMod, chat = summary(baseModel)$dispersion)}
         if(baseModel$family[1]=="quasibinomial"){
     BinMod<-update(baseModel, family=binomial)
           fitStat <- QAICc(BinMod, chat = summary(baseModel)$dispersion)}
   }
  
 # fitStat <- get.measure_2d(fitnessMeasure, NULL, baseModel)
 if(fitnessMeasure=="CV.offset"){    
#    if(dim(model.matrix(baseModel))[2]==1){
#       data2<- data.frame(response=response)
#       textForEval<- "tempCVFit<-glm(response~1, data=data2)" 
#             
#       }
#     
#     if(dim(model.matrix(baseModel))[2]>1){
#       data2<- data.frame(response=response, model.matrix(baseModel)[,2:length(coefficients(baseModel))])
#       
#   names(data2)<- c("response", paste("V", 1:(length(coefficients(baseModel))-1), sep=""))
#   textForEval<- paste("tempCVFit<-glm( response ~ ", paste("V", 1:(length(coefficients(baseModel))-1), sep="", collapse="+"), ", data=data2)")
     #}
        
#   eval(parse(text=textForEval))  
#   require(boot)
#   fitStat<-cv.glm(data2,tempCVFit, K=5)$delta[2]
# 
  fitStat<- mean(getCV_type2(folds = 5, baseModel))
  }
  

  if(fitnessMeasure=="CV"){  
    
    nfolds<-length(unique(data$foldid))
    store<- matrix(0, nrow=nfolds, ncol=1)
    d2k<-dists
    for(f in 1:nfolds){
      dists<-d2k[data$foldid!=f,]
      foldedFit<- update(baseModel, .~., data=data[data$foldid!=f,])
      if(length(coef(foldedFit))==1){
        dists<-d2k[data$foldid==f,]
        predscv<- exp(as.matrix(model.matrix(baseModel)[data$foldid==f])%*%coef(foldedFit)) * exp(baseModel$offset)[data$foldid==f]
      }else{  
        dists<-d2k[data$foldid==f,]
        predscv<- exp(model.matrix(baseModel)[data$foldid==f,]%*%coef(foldedFit)) * exp(baseModel$offset)[data$foldid==f]
      }        
      store[f]<- mean((data$response[data$foldid==f]-predscv)**2)
    }
    dists<-d2k
    fitStat<- mean(store)
  }
# 
  
  if(fitnessMeasure=="QICb"){       
          fitStat <- QICb(baseModel$data$response, fitted(baseModel), length(baseModel$coeff), length(baseModel$data$response))
    }
  
  
  
  #cat("Evaluating new fit: ", fitStat, "\n")
  if(is.na(fitStat)){
    fitStat<- tempMeasure + 1000
    cat("Change Fit due to NA: ", fitStat, "\n")
  }
 # output<-fit.thinPlate_2d(fitnessMeasure,dists,invInd[aR],radii,baseModel,radiusIndices,models)
 
  output = fit.thinPlate_2d(fitnessMeasure, dists,invInd[aR],radii, baseModel,radiusIndices,models, fitStat, interactionTerm, data)
  
  out.lm<-output$currentModel
  models<-output$models
  print("Initial model fitted...")
  point <- mapInd[-posKnots]
  position<- cbind()
  for (j in 1:length(knotPoint)) {
    position[knotPoint[j]]<- 0
  }    
  for (j in 1:length(point)) {
    position[point[j]]<- j
  }
  
  measures = 0
  BIC<-get.measure_2d(fitnessMeasure,measures,out.lm, data,  dists, invInd[aR],radii,radiusIndices)$fitStat
  
  #print(BIC[length(BIC)])

print("Fitting Initial Radii")
out<-choose.radii(BIC,1:length(radiusIndices),radiusIndices,radii,out.lm,dists,invInd[aR],baseModel,fitnessMeasure,response,models, interactionTerm, data)
BIC=out$BIC
radiusIndices=out$radiusIndices
out.lm=out$out.lm
models = out$models

  print("initialising complete")
  
  
  ####track <- rbind(track,cbind("init",t(aR),BIC[length(BIC)],adjRsq[length(adjRsq)],GCV[length(GCV)]))
  list(point=point,knotPoint=knotPoint,position=position,aR=aR,BIC=BIC,track=track,out.lm=out.lm,invInd=invInd,
        radiusIndices=radiusIndices,models=models)
  
  
  
}


################################################################################################################
"exchange.step_2d" <- function(gap,knotDist,radii,invInd,dists,explData,response,explanatory,maxIterations,fitnessMeasure, point,knotPoint,position,aR,BIC,track,out.lm,improveEx, maxKnots,tol=0,baseModel,radiusIndices,models, interactionTerm, data){

	attributes(baseModel$formula)$.Environment<-environment()
   
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
      new<-scale(explanatory[point,],center=c(explData[index,1],explData[index,2]))
      ####Pick nearest grid point that is also far enough away from another knot
      legPos<-position[which(apply(knotDist[invInd[point],invInd[aR]],1,min)>=gap)]
      if (length(legPos)>0) {
         index<-which.min(abs(new[,1])+abs(new[,2]))
         if (!(any(knotDist[invInd[point[index]],invInd[aR]]<gap))) {
           output <- move.knot_2D(radii,invInd,dists,explData,index,fitnessMeasure,BIC,aR,point,
                                      response,explanatory,out.lm,improve,improveEx, track,
                                      maxKnots,tol,baseModel,radiusIndices,models, interactionTerm, data)
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
           BIC <- rbind(BIC,output$fitStat)
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

################################################################################################################

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
              print("move ***********************************")
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


####################################################################################################################

"improve.step_2d" <- function(gap,knotDist,radii,invInd,dists,gridResp,grid,explData,xVals,yVals,num,response,explanatory,maxIterations,fitnessMeasure, point,knotPoint,position,aR,BIC,track,out.lm,improveNudge,tol=0,baseModel,radiusIndices,models, interactionTerm, data){
	attributes(baseModel$formula)$.Environment<-environment()
   print("******************************************************************************")
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
          if (grid[knotPoint[i],1] > 1) {
             if (!(is.na(gridResp[knotPoint[i] - 1]))) {
                nhbrs<- c(nhbrs, knotPoint[i] - 1)
             }
             if ((grid[knotPoint[i],2] > 1) && !(is.na(gridResp[knotPoint[i] - 1 - xVals]))) {
                nhbrs<-c(nhbrs,knotPoint[i] - 1 - xVals)
             }
             if ((grid[knotPoint[i],2] < yVals) && !(is.na(gridResp[knotPoint[i] - 1 + xVals]))) {
	          nhbrs<-c(nhbrs,knotPoint[i] - 1 + xVals)
             }
          }
          if (grid[knotPoint[i],1] < xVals) {
             if (!(is.na(gridResp[knotPoint[i] + 1]))) {
                nhbrs<- c(nhbrs, knotPoint[i] + 1)
             }
             if ((grid[knotPoint[i],2] > 1) && !(is.na(gridResp[knotPoint[i] + 1 - xVals]))) {
	          nhbrs<-c(nhbrs,knotPoint[i] + 1 - xVals)
	       }
	       if ((grid[knotPoint[i],2] < yVals) && !(is.na(gridResp[knotPoint[i] + 1 + xVals]))) {
	          nhbrs<-c(nhbrs,knotPoint[i] + 1 + xVals)
             }
          }
          if ((grid[knotPoint[i],2] > 1) && !(is.na(gridResp[knotPoint[i] - xVals]))) {
             nhbrs<- c(nhbrs, knotPoint[i] - xVals)
          }
          if ((grid[knotPoint[i],2] < yVals) && !(is.na(gridResp[knotPoint[i] + xVals]))) {
             nhbrs<- c(nhbrs, knotPoint[i] + xVals)
          }
          for (j in nhbrs){
             if (is.na(position[j])) {
                #print(j)
                #print(knotPoint[i])
                #print(xVals)
             }
#cat(invInd[otherKnots], '\n', invInd[j], '\n', aR, '\n')
           
                if ((length(otherKnots)==0) || ( min(knotDist[invInd[j],invInd[otherKnots]])>=gap)) {
                   tempR<-aR
                   tempR[i]<-j

                   output<-fit.thinPlate_2d(fitnessMeasure, dists,invInd[tempR],radii,baseModel,radiusIndices,models, fitStat, interactionTerm, data)
                   initModel<-output$currentModel
                   models<-output$models
                   initBIC<-get.measure_2d(fitnessMeasure,fitStat,initModel, data,  dists, invInd[tempR],radii,radiusIndices)$fitStat


         out<-choose.radii(initBIC,i,radiusIndices,radii,initModel,dists,invInd[tempR],baseModel,fitnessMeasure,response,models, interactionTerm, data)
         tempRadii=out$radiusIndices
         tempOut.lm=out$out.lm
         models=out$models
         output<-get.measure_2d(fitnessMeasure,fitStat,tempOut.lm, data,  dists, invInd[tempR],radii,tempRadii)

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
 
 ####################################################################################################################################
 
 "drop.step_2d" <- function(radii,invInd,dists,explData,response,explanatory,maxIterations,fitnessMeasure,
                              point,knotPoint,position,aR,BIC,track,out.lm,improveDrop,minKnots,tol=0,baseModel,radiusIndices,models, interactionTerm, data) {
       
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
          output<-fit.thinPlate_2d(fitnessMeasure, dists,invInd[tempR],radii,baseModel,tempRadii,models, fitStat, interactionTerm, data)
          initModel<-output$currentModel
          models<-output$models
          initBIC<-get.measure_2d(fitnessMeasure,BIC,initModel, data,  dists, invInd[tempR],radii, tempRadii)$fitStat
         out<-choose.radii(initBIC,1:length(tempRadii),tempRadii,radii,initModel,dists,invInd[tempR],baseModel,fitnessMeasure,response,models, interactionTerm, data)
         tempRadii=out$radiusIndices
         tempOut.lm=out$out.lm
         models=out$models
         output<-get.measure_2d(fitnessMeasure,fitStat,tempOut.lm, data,  dists, invInd[tempR],radii,tempRadii)

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
 
 
 ###################################################################################
 
  "get.measure_2d" <- function(fitnessMeasure,measures,out.lm, data, dists,aR,radii,radiusIndices){
  #print("ooooooooooooooooooooooooooooooooooooooo")
  #print("Getting measure...")
  #print("ooooooooooooooooooooooooooooooooooooooo")
  
  attributes(out.lm$formula)$.Environment<-environment()

  tempMeasure <- measures[1]
  if(fitnessMeasure=="AIC"){       
          fitStat <- AIC(out.lm)}
  
  if(fitnessMeasure=="AICc"){       
          fitStat <- AICc(out.lm)}
  
  if(fitnessMeasure=="BIC"){       
          fitStat <- BIC(out.lm)}
  
  if(fitnessMeasure=="QAIC"){       
    
    if(baseModel$family[1]=="quasipoisson"){
    PoisMod<-update(out.lm, round(.)~., family=poisson)
          fitStat <- QAIC(PoisMod, chat = summary(out.lm)$dispersion)}
    
    if(baseModel$family[1]=="quasibinomial"){
    BinMod<-update(out.lm, family=binomial)
          fitStat <- QAIC(BinMod, chat = summary(out.lm)$dispersion)}
  }
    
  if(fitnessMeasure=="QAICc"){       
     if(baseModel$family[1]=="quasipoisson"){
    PoisMod<-update(out.lm, family=poisson)
          fitStat <- QAICc(PoisMod, chat = summary(out.lm)$dispersion)}
    
  if(baseModel$family[1]=="quasibinomial"){
    BinMod<-update(out.lm, family=binomial)
          fitStat <- QAICc(BinMod, chat = summary(out.lm)$dispersion)}
     
  }
  
  if(fitnessMeasure=="cv.offset"){
#     if(dim(model.matrix(out.lm))[2]==1){
#       data2<- data.frame(response=response)
#       textForEval<- "tempCVFit<-glm(response~1, data=data2, family=family(out.lm))" 
#     }
#     if(dim(model.matrix(out.lm))[2]>1){
#       data2<- data.frame(response=response, model.matrix(out.lm)[,2:length(coefficients(out.lm))], offset = exp(baseModel$offset))
#       names(data2)<- c("response", paste("V", 1:(length(coefficients(out.lm))-1), sep=""), "offset")
#       textForEval<- paste("tempCVFit<-glm(round(response) ~ ", paste("V", 1:(length(coefficients(out.lm))-1), sep="", collapse="+"), ", family=family(out.lm), data=data2, offset = log(offset))")
#     }
#     eval(parse(text=textForEval))  
#     require(boot)
    #fitStat<-cv.glm(data2,tempCVFit, K=5)$delta[2]
    fitStat<- mean(getCV_type2(folds = 5, out.lm))
  }
    
  
  if(fitnessMeasure=="CV"){  
    
    nfolds<-length(unique(data$foldid))
    store<- matrix(0, nrow=nfolds, ncol=1)
    d2k<-dists
    for(f in 1:nfolds){
      dists<-d2k[data$foldid!=f,]
      foldedFit<- update(out.lm, .~., data=data[data$foldid!=f,])
      if(length(coef(foldedFit))==1){
        dists<-d2k[data$foldid==f,]
        predscv<- exp(as.matrix(model.matrix(out.lm)[data$foldid==f])%*%coef(foldedFit)) * exp(baseModel$offset)[data$foldid==f]
      }else{  
        dists<-d2k[data$foldid==f,]
        predscv<- exp(model.matrix(out.lm)[data$foldid==f,]%*%coef(foldedFit)) * exp(out.lm$offset)[data$foldid==f]
      }        
      store[f]<- mean((data$response[data$foldid==f]-predscv)**2)
    }
    
    fitStat<- mean(store)
  }
  
  # calculate a QIC with bayesian penalty
  if(fitnessMeasure=="QICb"){       
    fitStat <- QICb(out.lm$data$response, fitted(out.lm), length(out.lm$coeff), length(out.lm$data$response))
  }
  
 # cat("Evaluating new fit: ", fitStat, "\n")
  if(is.na(fitStat)){
    fitStat<- tempMeasure + 1000
    cat("Change Fit due to fitStat=NA: ", fitStat, "\n")
  }
  #if(length(which(is.na(out.lm$coefficients)))>0){
  #  fitStat<- tempMeasure + 1000
  #  cat("Change Fit due to NA coefficients: ", fitStat, "\n")
  #}
  
  #print("ooooooooooooooooooooooooooooooooooooooo")
  #print("Got measure...")
  #print("ooooooooooooooooooooooooooooooooooooooo")
  list(tempMeasure=tempMeasure,fitStat=fitStat)
 }
 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# getCVids <- function(data, folds){
#   N <- 1:nrow(data)
#   n_cv <- round(length(N)/folds)
#   set.seed(1234)
#   id_cv <- sample(rep(1:folds, n_cv), n_cv*folds)
#   id_cv <- id_cv[1:length(N)]
#   return(id_cv)
# }

getCV_type2<- function(folds, baseModel){
  data<- baseModel$data
  data$response<- baseModel$y
  if(baseModel$family[[1]]=="poisson"){data$response<- round(baseModel$y)}
  data$density <- data$response/exp(baseModel$offset)
  cvscore <- vector(length=folds)
  id_cv<- getCVids(data, folds)
  
  for(i in 1:folds){
    #cat('fold: ', i, '\n')
    tempid<- which(id_cv!=i)
    data2<- data.frame(response=baseModel$y[tempid], model.matrix(baseModel)[tempid,2:length(coefficients(baseModel))], offset=baseModel$offset[tempid])
    names(data2)<- c("response", paste("V", 1:(length(coefficients(baseModel))-1), sep=""), "offset")
    textForEval<- paste("tempCVFit<-glm(response ~ ", paste("V", 1:(length(coefficients(baseModel))-1), sep="", collapse="+"), ", family=family(baseModel), data=data2, offset = offset)")
    eval(parse(text=textForEval))
    # make predictions to new model  
    idpred<- which(id_cv==i)
    newdata<- data.frame(model.matrix(baseModel)[idpred,2:length(coefficients(baseModel))], offset = rep(1, length(idpred)))
    names(newdata)<- c(paste("V", 1:(length(coefficients(baseModel))-1), sep=""), "offset")
    # make predictions to values ==i
    preds<- predict(tempCVFit, newdata, type='response')
    # find cv score between the data - turned into a density using offset and the predictions (calc as a density)
    cvscore[i] <- sum((data$density[idpred] - preds)^2)/length(preds)
  }
  return(cvscore)
}
#QICinternal<- function(model,response){
#  y<-response[,1]
#  p<-predict(model, type="response")
#  n<- response[,2]+response[,1]
#  if (!(length(y)==length(p))) {
#    ohmy<<-model
#    print("*****")
#  }
#  ql<- sum(y*log(p)+(n-y)*log(1-p))
#  npar<- length(coef(model))
#  -2*ql+2*npar
#}
#
###########################################################################

  "fit.thinPlate_2d" <- function(fitnessMeasure, dists,aR,radii,baseModel,radiusIndices,models, currentFit, interactionTerm, data) {

   attributes(baseModel$formula)$.Environment<-environment()
  
  #print("ooooooooooooooooooooooooooooooooooooooo")
  #print("Fitting Model...")
  #print("ooooooooooooooooooooooooooooooooooooooo")
    #aR<- aR
    #radiusIndices<- radiusIndices

  #cat("R Indices: ", radiusIndices, "\n")                         

  
 # print(aR)
  bspl<- "LocalRadialFunction(radiusIndices, dists, radii, aR)"
  
 
  if(is.null(interactionTerm)){  
    test<-paste("update(baseModel, .  ~ . + ",bspl, ")",sep="")
    currentModel<-eval(parse(text=test))
  }else{
    test<-paste("update(baseModel, .  ~  . + ",bspl, "*",interactionTerm, ")", sep="")
    currentModel<-eval(parse(text=test)) 
  }
  

  tempFit <- get.measure_2d(fitnessMeasure, currentFit, currentModel,data, dists,aR,radii,radiusIndices)$fitStat
  if(tempFit <= (currentFit+10)){
     models[[length(models)+1]] = list(aR,radiusIndices, radii, tempFit)
  }
 modelinprogress<<-currentModel
  #print("ooooooooooooooooooooooooooooooooooooooo")
  #print("Model fitted...")
  #print("ooooooooooooooooooooooooooooooooooooooo")
 return(list(currentModel=currentModel,models=models))
 
}
###########################################################################
# local radial function in another file

##############################################################################
choose.radii <- function(currentFit,indices,radiusIndices,radii,out.lm,dists,
                         aR,baseModel,fitnessMeasure,response,models, interactionTerm, data) {
  #print("+++++++++++++++++++++++++++")
  #print("Fitting Radii")
  #print("+++++++++++++++++++++++++++")
  iterations = 0
  bestRadii=radiusIndices
  currentModel=out.lm
  improving = 1
  last = rep(0,length(radiusIndices))
  if(length(radii)> 1){
  #print("Fitting Radii")
  #cat('Start: ', radiusIndices, '\n')
  while (improving) {
    iterations = iterations + 1
    #print(c("iteration: ",iterations))
    improving = 0
    for (i in indices) {
      tempRadii = radiusIndices
      if ((tempRadii[i] > 1)&(last[i] <= 0)) {
        thisImprove = 1
        while (thisImprove) {
          thisImprove = 0
          tempRadii[i] = tempRadii[i] - 1
          output<- fit.thinPlate_2d(fitnessMeasure, dists,aR,radii,baseModel,tempRadii,models, currentFit, interactionTerm, data)
          tempModel<-output$currentModel
          models<-output$models
          tempFit <- get.measure_2d(fitnessMeasure,currentFit,tempModel, data, dists,aR,radii, tempRadii)$fitStat
          #print(c(i,"Down",tempFit))
          if (tempFit < currentFit) {
            #print("UPDATING")
            bestRadii = tempRadii
            currentFit = tempFit
            currentModel = tempModel
            improving = 1
            if (tempRadii[i] > 1) thisImprove = 1
            last = rep(0,length(radiusIndices))
            last[i] = -1
          }
        }
      }
      tempRadii = radiusIndices
      if ((tempRadii[i]<length(radii))&(last[i] >= 0)) {
        thisImprove = 1
        while (thisImprove) {
          thisImprove = 0
          tempRadii[i] = tempRadii[i] + 1
          output<- fit.thinPlate_2d(fitnessMeasure, dists,aR,radii,baseModel,tempRadii,models, currentFit, interactionTerm, data)
          tempModel<-output$currentModel
          models<-output$models
          tempFit <- get.measure_2d(fitnessMeasure,currentFit,tempModel, data,  dists,aR, radii, tempRadii)$fitStat
          #print(c(i,"Up",tempFit))
          if (tempFit < currentFit) {
           # print("UPDATING")
            bestRadii = tempRadii
            currentFit = tempFit
            currentModel = tempModel
            improving = 1
            if (tempRadii[i]<length(radii)) thisImprove = 1
            last = rep(0,length(radiusIndices))
            last[i] = 1
          }
        }
      }
    }
    #print(c("Current fit:",currentFit))
    radiusIndices = bestRadii
  }
  }# if length radii > 1 loop
 # print("+++++++++++++++++++++++++++")
  #print("Fitted Radii")
  #print("+++++++++++++++++++++++++++")
  list(BIC=currentFit,out.lm=currentModel,radiusIndices=radiusIndices,models=models)
}
##############################################
# this function takes a list object containing models fitted and reduces the 
# number to only those within a delta FitStat (eg aic, bic...) of 10

thinModels<- function(models){

      bic = vector(length=length(models))

      for(i in 1:length(models)){
        bic[i] = models[[i]][[4]]
		  }
	
      id = which(bic<=(min(bic)+10))
      tempModels = models[id]
      
      return(models = tempModels)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~ QICb ~~~~~~~~~~~~~~~~~~~~~~
# function to calculate a quasi version of a bic
# data inputs:
# data: vector of data
# fits: vector of fitted values from a model
# k: number of parameters estimated (length(coefficients))
# n: number of data points

QICb<-function(data, fits, k, n){
  ql<-sum(data*log(fits) - fits)
  -2*ql + k*log(n)
}
