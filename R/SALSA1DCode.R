#' Code for adaptively spacing knots for a given covariate.
#' 
#' @param response vector of response data for the modelling process
#' @param explanatory vector of covariate to find knots for
#' @param degree degree of the spline to be used
#' @param minKnots minimum number of knots to fit
#' @param maxKnots maximum number of knots to fit
#' @param startKnots number of equally spaced knots to start with (between minKnots and maxKnots)
#' @param gap minimum gap between knots (in unit of measurement of \code{explanatory})
#' @param fitnessMeasure (default=BIC). Measure used to evaluate the fit. Other options are AIC, AICc, BIC, QAIC, QAICc, QICb (Quasi-Likelihood Information Criterion with log(n) penalty)
#' @param maxIterations exchange/improve heuristic will terminate after \code{maxIterations} if still running
#' @param initialise (default = TRUE). Logical stating whether or not to start with equally spaced knots (TRUE) or user specified locations (FALSE)
#' @param initialKnots If \code{initialise=FALSE} then the start locations for the knots are specified in \code{initialKnots}
#' @param baseModel starting model for SALSA to use.  Must not contain the covariate in \code{explanatory}
#' @param bd the x-coordinate of the boundary knots of \code{explanatory}
#' @param spl "bs" uses b-spline, "cc" uses cyclic cubic, "ns" uses natural cubic spline for fitting smooth to \code{explanatory}
#' @param winHalfWidth Half-width of window used to calculate region with biggest average residual magnitude
#' @param interactionTerm character stating the variable to use as an interaction
#' @param suppress.printout \code{default=FALSE}. If TRUE, progress is printed into the workspace. If FALSE, a .log file is created in the working directory.  
#' @param cv.opts A list object containing options for \code{cv.gamMRSea}.
#' 
#' @author Cameron Walker, Department of Engineering Science, University of Auckland, University of Auckland), Lindesay Scott-Hayward (University of St Andrews)
#' 
#' 
#'@export
#'

"return.reg.spline.fit" <- function(response,explanatory,degree,minKnots,maxKnots,startKnots,gap,winHalfWidth,fitnessMeasure="BIC", maxIterations=100, initialise = TRUE, initialKnots = NULL, baseModel=NULL, bd, spl,interactionTerm=interactionTerm, suppress.printout=FALSE, cv.opts, splineParams){

  # if(suppress.printout){
  #   sink(file = 'salsa1d.log')
  # }
  # 
  varWinHW=5
  computeWt=0
  wts=rep(1,length(explanatory))
  maxSites=200
# requires splines library and mgcv library to be loaded!!

# PARAMETERS
# minKnots:        minimum number of knots to fit
# maxKnots:        maximum number of knots to fit
# startKnots:      number of equally spaced knots to sstart with (between minKnots and maxKnots)
# gap:             minimum gap between knots (i.e. number of data points)
# winHalfWidth:    half-width of window used to calculate region with biggest average residual magnitude
# fitnessMeasure:  measure used to evaluate the fit (value = 1--4)
# fitnessMeasure=="AIC" uses AIC       
# fitnessMeasure=="AICc" uses AICc       
# fitnessMeasure=="BIC" uses BIC       
# fitnessMeasure=="QAIC" uses QAIC       
# fitnessMeasure=="QAICc" uses QAICc
# fitnessMeasure== "CV.offset" uses CV with an offset - no blocking
# fitnessMeasure== "CV.glm" uses CV with no offset included
# fitnessMeasure== "CV" uses CV with an offset with blocking structure
#
# maxIterations:   exchange/improve heuristic will terminate after maxIterations if still running
# varWinHW:        used for determining heteroscedastic weights
# computeWt:       whether to use weights for heteroscedasticity - default is No
# wts:             for heteroscedasticity, and bootstrapping
# bd:              the x-coordinate of the boundary knots
# spl:             "bs" uses b-spline, "cc" uses cyclic cubic, "ns" uses natural cubic
 #
# output:          aR        - list of knot points for each iteration
#                  arSq      - list of adjusted r-squareds for each it.
#                  finaR     - list of optimal knots for each number of knots
#                  finarSq   - list of adjusted r-squareds for each number of knots
#                  finBIC    - list of BICs for each number of knots

# pointers:        knotPoint - the index of the knot points (i, where explanatory[i] is a knot)
#                  point     - the index of the other points
#                  position  - the index in point of ith data point
#                              0 otherwise (position[i] = j,where point[j] = i)

# sort out baseModel formulae environment issue
  if (isS4(baseModel)) {
    attributes(baseModel@misc$formula)$.Environment<-environment()
    data <- baseModel@data
    data$explanatory<-explanatory
    baseModel@data <- data
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
    data<- baseModel$data
    data$explanatory<-explanatory
    baseModel<-update(baseModel, data=data)
  }

# LSH 12/3/15 added dispersion parameter calc
initDisp<-getDispersion(baseModel)
if(length(unique(response))!=2){
  print(paste('initialDispersion ', initDisp, sep=''))
}
    
####deal with multiple unordered x-values
knotSites <- cbind(sort(explanatory), rep(1, length(explanatory)))
knotSites <- knotSites[which(duplicated(knotSites)==F),]
if (nrow(knotSites) > maxSites) {
 if(nrow(knotSites)>800){
   knotSites <- knotSites[sample(1:nrow(knotSites), 800),]
   knotSites = sort(cover.design(knotSites, nd=maxSites)$design[,1])
 }else{
   knotSites = sort(cover.design(knotSites, nd=maxSites)$design[,1]) 
 }
  # LSH updated 19/2/15 so that the candidate knot locations may only be at data locations.
  #quantile(knotSites,probs=seq(0,1,length=maxSites),na.rm=TRUE,names=FALSE)
}
#print(knotSites)

# remove locations for knots if they are also boundary knots????
if(length(which(knotSites==bd[1]))>0){
   knotSites <- knotSites[-(which(knotSites<=(bd[1] + gap)))]
}
if(length(which(knotSites==bd[2]))>0){
   knotSites <- knotSites[-(which(knotSites>=(bd[2] - gap)))]
}

#print(knotSites)

    ###########################initialisation######################################
    output <- initialise.measures(startKnots, explanatory, response, degree, wts, initialise, initialKnots,baseModel,knotSites, bd, spl, fitnessMeasure, interactionTerm, initDisp, cv.opts, splineParams=splineParams)
    point <- output$point
    knotPoint <<- output$knotPoint
    position <- output$position
    aR <- output$aR
    measures <- output$measures # comes from update measures
    out.lm <- output$out.lm
    models <- output$models

        print("^^^^^^^^^^^Initial^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        print(get.measure(fitnessMeasure,measures,out.lm, initDisp, cv.opts)$fitStat)

    ###################################algorithm loop#############################
    improveEx <- 1
    improveNudge <- 1
    improveDrop <- 1
    overallImprove = 0
    while (improveEx | improveNudge | improveDrop) {
      improveEx <- 0
      improveNudge <- 0
      improveDrop <- 0
    ###################################exchange step#############################
      output <- exchange.step(degree, gap, response,explanatory,maxIterations,fitnessMeasure,point,knotPoint,position,aR,measures,
                                 out.lm,improveEx,maxKnots,winHalfWidth,wts, baseModel,knotSites,models, bd, spl, interactionTerm , initDisp, cv.opts, splineParams=splineParams)
      point <- output$point
      knotPoint <- output$knotPoint
      position <- output$position
      aR <- output$aR
      measures <- output$measures
      out.lm <- output$out.lm
      models <- output$models
      improveEx <- output$improveEx
    
      if (improveEx) {
        print("^^^^^^^^^^^Exchange^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
       # print(get.measure(fitnessMeasure,measures,out.lm)$fitStat)
      }
    #####################################improve step############################
      output <- improve.step(degree, gap, length(aR), response,explanatory,maxIterations,fitnessMeasure,point,knotPoint,position,
                                aR,measures,out.lm,improveNudge, wts, baseModel,knotSites,models, bd, spl, interactionTerm, initDisp, cv.opts, splineParams=splineParams)
      point <- output$point
      knotPoint <- output$knotPoint
      position <- output$position
      aR <- output$aR
      measures <- output$measures
      out.lm <- output$out.lm
      models <- output$models
      improveNudge <- output$improveNudge
      if (improveNudge) {
        print("^^^^^^^^^^^Improve^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
       # print(get.measure(fitnessMeasure,measures,out.lm)$fitStat)
      }
    ###################################drop step#################################
      if (length(aR) > minKnots) {
         output <- drop.step(degree, response,explanatory,maxIterations,fitnessMeasure,point,knotPoint,position,aR,measures,
                                  out.lm,improveDrop,minKnots, wts, baseModel,models, bd, spl, interactionTerm, initDisp, cv.opts, splineParams=splineParams)
         point <- output$point
         knotPoint <- output$knotPoint
         position <- output$position
         aR <- output$aR
         measures <- output$measures
         out.lm <- output$out.lm
         models <- output$models
         improveDrop <- output$improveDrop    
      if (improveDrop) {
        print("^^^^^^^^^^^Drop^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
       # print(get.measure(fitnessMeasure,measures,out.lm)$fitStat)
      }
         }
       if ((improveEx) | (improveNudge) | (improveDrop)) overallImprove = 1
      }
  print("And we're done...")
  list(output=c(length(aR),measures,aR),aR=aR,out.lm=out.lm,improve=overallImprove,knotSites=knotSites,models=models)
}

########################################################################################################################

"initialise.measures" <- function(num,explanatory, response, degree, wts,  initialise, initialKnots,baseModel,knotSites, bd, spl, fitnessMeasure, interactionTerm, initDisp, cv.opts, splineParams=splineParams){
   
  if (isS4(baseModel)){
    attributes(baseModel@misc$formula)$.Environment<-environment()
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
  }
  
    print("Initialising...")
    models = vector("list",0)
    if (initialise) {
      step<-(knotSites[length(knotSites )]-knotSites[1])/(num + 1)
      aR <- c()
      knotPoint<-c() 
         for (j in 1:num) {
            knotPoint[j] <- floor(j*length(knotSites )/(num+1))
            aR <- c(aR,knotSites [knotPoint[j]])
         }

      } else {
        aR = initialKnots
        unsortedKP=c()
        for (i in 1:length(aR)) {
          newKnot = which(abs(knotSites-aR[i])==min(abs(knotSites-aR[i])))
          unsortedKP = c(unsortedKP,newKnot)
        }
        knotPoint = sort(unique(unsortedKP))
        aR = knotSites[knotPoint]
        num=length(aR)
      }
      point<<-1:length(knotSites)
      point<-point[-knotPoint]
      position<-1:(knotPoint[1]-1)
      if (num > 1) {
         for (j in 2:num) {
           if((knotPoint[j]-knotPoint[j-1])>1){
             position<-c(position,0,(knotPoint[j-1]-(j-2)):(knotPoint[j]-j))
           } else {
             position<-c(position,0)
           }
         }
      }
    position<-c(position,0,(knotPoint[num]-(num-1)):(length(knotSites )-num))
    output <- fit.model(explanatory,degree,aR,baseModel,models, bd, spl, fitnessMeasure, interactionTerm, initDisp, cv.opts, splineParams)
    out.lm<-output$currentModel
    models<-output$models
    #model.out<<-out.lm
    measures<- get.measure(fitnessMeasure,measures=NA,out.lm, initDisp, cv.opts)$fitStat
    #measures <- update.measures(out.lm)
cat("Initial fit = ", measures, aR,"\n")
print("initialisation complete...")
    list(point=point,knotPoint=knotPoint,position=position,aR=aR,measures=measures,out.lm=out.lm)
}

######################################################################################################################

"exchange.step" <- function(degree, gap, response,explanatory,maxIterations,fitnessMeasure,point,knotPoint,position,aR,
                               measures,out.lm,improveEx,maxKnots,winHalfWidth,wts, baseModel,knotSites,models, bd, spl, interactionTerm, initDisp, cv.opts,
                            splineParams){
  
  if (isS4(baseModel)){
    attributes(baseModel@misc$formula)$.Environment<-environment()
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
  }
  
  # Loop - fuse used to ensure algorithm terminates
  print("Exchanging...")
  fuse <- 0
  improve <- 1
  while ( (improve) & (fuse < maxIterations) ) {
    fuse <- fuse + 1
    improve <- 0
    output <- locate.max.res(point,position,gap,response,explanatory, bd, winHalfWidth,out.lm,knotPoint,aR,wts, knotSites, spl)
    index <- output$index
    #if (length(index)>1) browser()
    if (length(index)>0) {
      if (index > 0) {
        output <- move.knot(degree, index,fitnessMeasure,measures,aR,point,response,explanatory,out.lm,improve,improveEx,
                            maxKnots, wts,  baseModel,knotSites,models, bd, spl, interactionTerm, initDisp, cv.opts, splineParams=splineParams)
        improve <- output$improve
        improveEx <- output$improveEx
        models <-output$models
      }
    }
    if (1 - improve) break
    out.lm<-output$out.lm
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
    measures <- output$tempMeasures
  }
  print("Exchanging done...")
  list(point=point,knotPoint=knotPoint,position=position,aR=aR,
       measures=measures,out.lm=out.lm,improveEx=improveEx,models=models)
}



###################################################################################################################

"locate.max.res" <- function(point,position,gap,response,explanatory, bd,winHalfWidth,out.lm,knotPoint,aR,wts,knotSites, spl){
   
  print("Locating maximum residual......")
  tempRes <<-residuals(out.lm,type="pearson")
  index <- NULL
  for (i in 1:length(knotPoint)) {
    if (isS4(out.lm)){
      tempRes[(explanatory <= knotSites[knotPoint[i]]+gap)&(explanatory >= knotSites[knotPoint[i]]-gap),] = 0
      tempRes[(explanatory <= bd[1] + gap),] = 0
      tempRes[(explanatory >= bd[2] - gap),] = 0
    } else {
      tempRes[(explanatory <= knotSites[knotPoint[i]]+gap)&(explanatory >= knotSites[knotPoint[i]]-gap)] = 0
	    tempRes[(explanatory <= bd[1] + gap)] = 0
	    tempRes[(explanatory >= bd[2] - gap)] = 0
    }
  }
  # 8/8/12 LSH
  # add if statement to prevent random allocation of maxInd if tempRes is all zeros
  #
  if (max(abs(tempRes))==0) {
		print("All Residuals Zero - no move")
	} else {
	  if (isS4(out.lm)){
	    maxInd = which.max(tempRes)
	  } else {
	    maxInd = max.col(t(abs(tempRes)))
	  }
    siteIndex = which(abs(knotSites-explanatory[maxInd])==min(abs(knotSites[position>0]-explanatory[maxInd])))[1]
    # added the [1] to only take the first nearest point.  only an issue if selected data location to place knot is exactly equidistant to knot locations available eitherside
    #if(length(siteIndex>1)) browser()
    index = which(point==siteIndex)
    point<<-point
    position<<-position
    knotSites<<-knotSites
    print("Maximum residual found...")
	}
  list(index=index)
}

################################################################################################################

"move.knot" <- function(degree, index,fitnessMeasure,measures,aR,point,response,explanatory,out.lm,improve,improveEx,maxKnots,
                            wts, baseModel,knotSites,models, bd, spl, interactionTerm, initDisp, cv.opts, splineParams){
  
  if (isS4(out.lm)) {
    attributes(baseModel@misc$formula)$.Environment<-environment()
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
  }
  
  print("Moving knot...")
  chck<-c()
  tempMeasures <- measures
  
  for (i in 1:length(aR)) {
    tempR<-aR
    tempR[i]<-knotSites[point[index]]
    chck<-rbind(chck,tempR)
    output <- fit.model(explanatory,degree,tempR,baseModel,models, bd, spl, fitnessMeasure, interactionTerm, initDisp, cv.opts, splineParams)
       
    if (isS4(out.lm)) {
      converge <- output$currentModel@iter < output$currentMode@control$maxit
    } else {
      converge <- output$currentModel$converge
    }
       
    if (converge) {
      tempOut.lm<-output$currentModel
      models<-output$models
      output<-get.measure(fitnessMeasure,measures,tempOut.lm, initDisp, cv.opts)
      tempMeasure<-output$tempMeasure
      fitStat<-output$fitStat
      chck<-rbind(chck,fitStat)
      if (tempMeasure > fitStat ) { ## Lindesay: removed convergence check - already tested 7 lines above)
      #print("t measure:")
      #print(tempMeasure)
      #print(tempR)
      #print("f stat:")
      #print(fitStat)
      out.lm <- tempOut.lm
      tempMeasures<- fitStat
      newR <- tempR
      tempKnot <- i
      improve <- 1
      improveEx <- 1
      }
    }      
  }
    
  if (length(aR)<maxKnots) {
    tempR<-c(aR,knotSites[point[index]])
    output <- fit.model(explanatory,degree,tempR,baseModel,models, bd, spl, fitnessMeasure, interactionTerm, initDisp, cv.opts, splineParams)
      
    if (isS4(out.lm)) {
      converge <- output$currentModel@iter < output$currentMode@control$maxit
    } else {
      converge <- output$currentModel$converge
    }	 
      
    if (converge) {
      tempOut.lm<-output$currentModel
      models<-output$models
      output<-get.measure(fitnessMeasure,measures,tempOut.lm, initDisp, cv.opts)
      tempMeasure<-output$tempMeasure
      fitStat<-output$fitStat
      if (tempMeasure > fitStat) {
        #print("t measure:")
        #print(tempMeasure)
        #print(tempR)
        #print("f stat:")
        #print(fitStat)
        out.lm <- tempOut.lm
        tempMeasures<-fitStat
        newR <- tempR
        tempKnot <- length(aR) + 1
        improve <- 1
        improveEx <- 1
      }
    }
  }

  print("Knot moved...")
  if (improve) {
    list(tempMeasures=tempMeasures,newR=newR,tempKnot=tempKnot,improve=improve,
         improveEx=improveEx, out.lm=out.lm,models=models)
  } else {
    list(improve=improve,improveEx=improveEx,models=models)
  }
}

####################################################################################################################

"improve.step" <- function(degree, gap, num,response,explanatory,maxIterations,fitnessMeasure,point,knotPoint,position,aR,measures,out.lm,improveNudge,wts, baseModel,knotSites,models, bd, spl, interactionTerm, initDisp, cv.opts, splineParams){
  if (isS4(baseModel)){
    attributes(baseModel@misc$formula)$.Environment<-environment()
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
  }

  print("Improving...")
   improve <- 1
   fuse <- 0
   while ( (improve) & (fuse < maxIterations) ) {
     fuse <- fuse + 1
     improve <- 0
     for (i in 1:length(aR)) {
       #browser()
       output <- local.shift.up(degree, out.lm,i,knotPoint,gap,position,fitnessMeasure,measures,aR,point,response,explanatory, improve,improveNudge,wts, baseModel,knotSites,models, bd, spl, interactionTerm, initDisp, cv.opts, splineParams)
       point <- output$point
       
       #print(length(point))
       ##print(length(unique(point)))
       
       position <- output$position
       knotPoint <- output$knotPoint
       aR <- output$aR
       measures <- output$measures
       out.lm<- output$out.lm
       models<-output$models
       shouldBreak <- output$shouldBreak
       improve <- output$improve
       improveNudge <- output$improveNudge
       if (shouldBreak) {break}
       output <- local.shift.down(degree, out.lm,i,knotPoint,gap,position,fitnessMeasure,measures,aR,point,response,explanatory, improve,improveNudge,wts,baseModel,knotSites,models, bd, spl, interactionTerm, initDisp, cv.opts, splineParams)
       point <- output$point
       
       #print(length(unique(point)))
       #print(length(point))
       
      position <- output$position
       knotPoint <- output$knotPoint
       aR <- output$aR
       measures <- output$measures
       out.lm<- output$out.lm
       models<-output$models
       shouldBreak <- output$shouldBreak
       improve <- output$improve
       improveNudge <- output$improveNudge
       if (shouldBreak) {break}
       }
   }
print("Improving complete...")
   list(point=point,knotPoint=knotPoint,position=position,aR=aR,measures=measures,
         out.lm=out.lm,improveNudge=improveNudge,models=models)
}

#########################################################################################################################

"drop.step" <- function(degree, response,explanatory,maxIterations,fitnessMeasure,point,knotPoint,position,aR,measures,out.lm,
                          improveDrop,minKnots,wts, baseModel,models, bd, spl, interactionTerm, initDisp, cv.opts, splineParams) {
  
  if (isS4(baseModel)) {
    attributes(baseModel@misc$formula)$.Environment<-environment()
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
  }
  
  print("Dropping...")
  improve<-0
  tempMeasures <- measures
  for (i in 1:length(aR)) {
    tempR <- aR
    tempR <- tempR[-i]
    output <- fit.model(explanatory,degree,tempR,baseModel,models, bd, spl, fitnessMeasure, interactionTerm, initDisp, cv.opts, splineParams)
    
    if (isS4(output$currentModel)) {
      converge <- output$currentModel@iter < output$currentMode@control$maxit
    } else {
      converge <- output$currentModel$converge
    }
    
    if (converge) {
      tempOut.lm<-output$currentModel
      models<-output$models
      output<-get.measure(fitnessMeasure,measures,tempOut.lm, initDisp, cv.opts)
      tempMeasure<-output$tempMeasure
      fitStat<-output$fitStat
      if (tempMeasure > fitStat) {
        #print("Drop t measure:")
        #print(tempMeasure)
        #print(tempR)
        #print("Drop f stat:")
        #print(fitStat)
        out.lm <- tempOut.lm
        tempMeasures<-fitStat
        newR <- tempR
        tempKnot <- i
        improve <- 1
        improveDrop <- 1
      }      
    }
  }
  if (improve) {
    aR <- newR
    point <- c(point,knotPoint[tempKnot])
    position[knotPoint[tempKnot]]<-length(point)
    knotPoint <- knotPoint[-tempKnot]
  }
  print("Dropped...")
  list(point=point,knotPoint=knotPoint,position=position,aR=aR,measures=tempMeasures,
       out.lm=out.lm,improveDrop=improveDrop,models=models)
}

         
###################################################################################

"local.shift.up" <- function(degree, out.lm, i,knotPoint,gap,position,fitnessMeasure,measures,aR,point,response,explanatory,
                                  improve,improveNudge,wts, baseModel,knotSites,models,bd , spl, interactionTerm, initDisp, cv.opts, splineParams){
  if (isS4(baseModel)){
    attributes(baseModel@misc$formula)$.Environment<-environment()
  }  else {
    attributes(baseModel$formula)$.Environment<-environment()
  }
  print("Shifting up...")
  shouldBreak <- 0
  move <- 1
  tempR <- aR
  while(move) {
    move <- 0
    #if (knotPoint[i] + gap + 1 <= length(knotSites)) {#have already removed illegal end knot sites at start of solve
    if (knotPoint[i] + 1 <= length(knotSites)) {#don't go over end!
      #check <- 1
      #for (j in 1:gap){
      #   check<-check*(position[knotPoint[i]+j]!=0)
      #}
      #if(knotPoint[i]==83) browser()
      otherPoints=knotPoint[-i]
      check = !any(abs(knotSites[otherPoints]-knotSites[knotPoint[i]+1])<=gap)
      # LSH added <= (19/2/15) rather than < so no two knots in same plae
      if (check) {
        tempR[i] <- knotSites[knotPoint[i]+1]
        output <- fit.model(explanatory,degree,tempR,baseModel,models, bd, spl, fitnessMeasure, interactionTerm, initDisp, cv.opts, splineParams)
        
        if (isS4(baseModel)){
          converge <- output$currentModel@iter < output$currentMode@control$maxit
        } else {
          converge <- output$currentModel$converge
        }
        
        if (converge) {
          tempOut.lm<-output$currentModel
          models<-output$models
          output<-get.measure(fitnessMeasure,measures,tempOut.lm, initDisp, cv.opts)
          tempMeasure<-output$tempMeasure
          fitStat<-output$fitStat
          if (fitStat < tempMeasure) {
            #browser()
            #print("Up t measure:")
            #print(tempMeasure)
            #print(tempR)
            #print("Up f stat:")
            #print(fitStat)
            shouldBreak <- 1
            improve <- 1
            improveNudge <- 1
            point[position[knotPoint[i]+1]] <- knotPoint[i]
            position[knotPoint[i]] <- position[knotPoint[i]+1]
            position[knotPoint[i]+1] <- 0
            knotPoint[i] <- knotPoint[i]+1
            aR <- tempR
            out.lm<- tempOut.lm
            measures<- fitStat
            move<-1
          }
        }
      }
    }
  }
  print("Up done...")
  list(out.lm=out.lm,point=point,position=position,knotPoint=knotPoint,aR=aR,measures=measures,
       shouldBreak=shouldBreak,improve=improve,improveNudge=improveNudge,models=models)
}
 
 ######################################################################################################################################
 
 "local.shift.down" <- function(degree, out.lm, i,knotPoint,gap,position,fitnessMeasure,measures,aR,point,response,explanatory,
                                     improve,improveNudge,wts, baseModel,knotSites,models, bd, spl, interactionTerm, initDisp, cv.opts, splineParams){
  
  if (isS4(baseModel)){
    attributes(baseModel@misc$formula)$.Environment<-environment() 
  } else {
    attributes(baseModel$formula)$.Environment<-environment() 
  }
  print("Shifting down...") 
  shouldBreak <- 0
  move <- 1
  tempR<-aR
  while (move) {
    move <- 0
    #if (knotPoint[i] - gap >= 2) {#have already removed illegal end knot sites at start of solve
    if (knotPoint[i]  >= 2) {#don't go over end!
      otherPoints=knotPoint[-i]
      check = !any(abs(knotSites[otherPoints]-knotSites[knotPoint[i]-1])<=gap)
      # LSH changed to <= from < so that two of the same knots do not get addded.
      #check <- 1
      #for (j in 1:gap){
      #   check<-check*(position[knotPoint[i]-j]!=0)
      #}    
      if (check) {
        tempR[i]<-knotSites[knotPoint[i]-1]
        output <- fit.model(explanatory,degree,tempR,baseModel,models, bd, spl, fitnessMeasure, interactionTerm, initDisp, cv.opts, splineParams)
        
        if (isS4(baseModel)) {
          converge <- output$currentModel@iter < output$currentMode@control$maxit
        } else {
          converge <- output$currentModel$converge
        }
        
        if (converge) {
          tempOut.lm<-output$currentModel
          models<-output$models
          output<-get.measure(fitnessMeasure,measures,tempOut.lm, initDisp, cv.opts)
	        tempMeasure<-output$tempMeasure
          fitStat<-output$fitStat             
          if (fitStat < tempMeasure) {
            #print("Down t measure:")
            #print(tempMeasure)
            #print(tempR)
            #print("Down f stat:")
            #print(fitStat)
            shouldBreak <- 1
            improve <- 1
            improveNudge <- 1
            point[position[knotPoint[i]-1]] <- knotPoint[i]
                 
            #print(length(unique(point)))
            # print(length(point))
            #if(length(unique(point))!=length(point)) browser()
                 
            position[knotPoint[i]] <- position[knotPoint[i]-1]
            position[knotPoint[i]-1] <- 0
            knotPoint[i] <- knotPoint[i]-1
            aR <- tempR
            out.lm<-tempOut.lm
            measures<- fitStat
            move<-1
          }
        }
      }
    }
  }
  print("Down done...")
  list(out.lm=out.lm,point=point,position=position,knotPoint=knotPoint,aR=aR,measures=measures,
       shouldBreak=shouldBreak,improve=improve,improveNudge=improveNudge,models=models)
}
 
 ####################################################################################################################################
 
 # ~~~~~~~~~~~~~~~~~


####
QICinternal<- function(model,response){
  y<-response[,1]
  p<-predict(model, type="response")
  n<- response[,2]+response[,1]
  if (!(length(y)==length(p))) {
    ohmy<<-model
    print("*****")
  }
  mu<- n*p
  ql<- sum(y*log(p)+(n-y)*log(1-p))
  PredictionsBasedOnModel<- predict(model, type="response")
  npar<- length(coef(model)) -2*ql+2*npar
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# functions needed for lsh cv calculation in fitnessmeasure
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
  if(baseModel$family[[1]]=="poisson"){data$response<- round(data$response)}
  data$density <- data$response/exp(baseModel$offset)
  cvscore <- vector(length=folds)
  id_cv<- getCVids(data, folds)
  
  for(i in 1:folds){
    cat('fold: ', i, '\n')
    tempid<- which(id_cv!=i)
    data2<- data.frame(response=data$response[tempid], model.matrix(baseModel)[tempid,2:length(coefficients(baseModel))], offset=baseModel$offset[tempid])
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
 ###################################################################################################################################
 
 "fit.model" <- function(explanatory,degree,aR,baseModel,models, bd, spl, fitnessMeasure, interactionTerm, initDisp, cv.opts, splineParams) {
 
  if (isS4(baseModel)){
    attributes(baseModel@misc$formula)$.Environment<-environment()
    data <- baseModel@data
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
    data<-baseModel$data
  }

   # ~~~~~~~~~~~~~~~~~~~~~
   # b-spline
  if(spl == 'bs'){
    #print("fitting model...")
    bspl<-paste("bs(explanatory, degree=", degree, ",Boundary.knots=c(",bd[1], ",", bd[2],"), knots= c(", sep="")
    if (length(aR)>1) {
        for (i in 1:(length(aR)-1)) {
             bspl<- paste(bspl, aR[i], ",", sep="")
        }
    }
    bspl<-paste(bspl, aR[length(aR)], ")",")",sep="")

    if(is.null(interactionTerm)){  
      test<-paste("update(baseModel, .  ~ . + ",bspl,")", sep="")
      out.lm<-eval(parse(text=test))
      if (isS4(out.lm)){
        out.lm@interactionterm <- interactionTerm
      }
    }else{
      test<-paste("update(baseModel, .  ~  . + ",bspl, "*",interactionTerm, ")", sep="")
      out.lm<-eval(parse(text=test))
      if (isS4(out.lm)){
        out.lm@interactionterm <- interactionTerm
      }
    }
    
    if (isS4(out.lm)) {
      converge <- out.lm@iter < out.lm@control$maxit
    } else {
      converge <- out.lm$converge
    }
    
    if(is.na(converge)) {converge = TRUE}
    #print(aR)
    #print("model fitted...")
    if (converge) {
      tempFit <- get.measure(fitnessMeasure, NA, out.lm, initDisp, cv.opts)$fitStat
      models[[length(models)+1]] = list(aR, tempFit)
    }
    return(list(currentModel=out.lm,models=models))
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~
  # cubic splines
  if(spl == 'cc'){
    #print("fitting model...")
    
    knots<- paste0(sort(c(bd[1],aR, bd[2])), collapse=",")
    ccspl <- paste0("cSplineDes(x=explanatory, knots=c(", knots, "), ord=", degree+1, ")")
    
    out.lm<- eval(parse(text=paste0("update(baseModel,. ~. +", ccspl,  ")")))
    
    if (isS4(out.lm)) {
      converge <- out.lm@iter < out.lm@control$maxit
    } else {
      converge <- out.lm$converge
    }
    
    if(is.na(converge)) {converge = TRUE}
    
    # out.lm<-update(baseModel, .  ~ . + splineString)
    
    #print(aR)
    #print("model fitted...")
    if (converge) {
      tempFit <- get.measure(fitnessMeasure, NA, out.lm, initDisp, cv.opts)$fitStat
      models[[length(models)+1]] = list(aR, tempFit)
    }
    return(list(currentModel=out.lm,models=models))
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~
  # natural cubics
  if(spl=="ns"){
    #print("fitting model...")
    #print(aR)
    bspl<-paste("ns(explanatory,Boundary.knots=c(",bd[1], ",", bd[2], "), knots= c(", sep="")
    if (length(aR)>1) {
      for (i in 1:(length(aR)-1)) {
        bspl<- paste(bspl, aR[i], ",", sep="")
      }
    }
    bspl<-paste(bspl, aR[length(aR)], ")",")",sep="")
    
    test<-paste("update(baseModel, .  ~ . + ",bspl,")", sep="")
    out.lm<-eval(parse(text=test))
    
    if (isS4(out.lm)) {
      converge <- out.lm@iter < out.lm@control$maxit
    } else {
      converge <- out.lm$converge
    }
    
    if(is.na(converge)) {converge = TRUE}

    #print("model fitted...")
    if (converge) {
      tempFit <- get.measure(fitnessMeasure, NA, out.lm, initDisp, cv.opts)$fitStat
      models[[length(models)+1]] = list(aR, tempFit)
    }
    
    # if(suppress.printout){
    #   sink()
    # }
    
    return(list(currentModel=out.lm,models=models))
  }
} # end of function

