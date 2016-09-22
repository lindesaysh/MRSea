#' Wrapper function for running SALSA2D
#'
#'
#'
#' @author Cameron Walker, Department of Enginering Science, University of Auckland.
#'
#' @export
#'


"return.reg.spline.fit.2d" <- function(splineParams, startKnots, winHalfWidth,fitnessMeasure="BIC", maxIterations=10, tol=0, baseModel=NULL, radiusIndices=NULL, initialise=TRUE, initialKnots=NULL, interactionTerm=NULL, knot.seed=10, suppress.printout=FALSE){
  #
  # if(suppress.printout){
  #   sink(file='salsa2d.log')
  # }


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

  # LSH 12/3/15 added dispersion parameter calc
  initDisp<-getDispersion(baseModel)
  print(paste('initialDispersion ', initDisp, sep=''))

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
  output <- initialise.measures_2d(knotDist,maxIterations,gap,radii,dists,gridResp,explData,startKnots,xvals, yvals, explanatory, response, baseModel, radiusIndices, initialise, initialKnots,fitnessMeasure, interactionTerm, data, knot.seed, initDisp)

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

  out.lm$splineParams[[1]]$knotPos<-aR
  out.lm$splineParams[[1]]$radiusIndices<-radiusIndices
  baseModel$splineParams<-out.lm$splineParams

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
    output <- exchange.step_2d(gap,knotDist,radii,invInd,dists,explData,response,explanatory,maxIterations,fitnessMeasure, point,knotPoint,position,aR,BIC,track,out.lm,improveEx,maxKnots,tol,baseModel,radiusIndices,models, interactionTerm, data, initDisp)
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

    out.lm$splineParams[[1]]$knotPos<-aR
    out.lm$splineParams[[1]]$radiusIndices<-radiusIndices
    baseModel$splineParams<-out.lm$splineParams

    ######################################improve step############################
    ####track <- rbind(track,cbind("improving",t(aR),BIC[length(BIC)],adjRsq[length(adjRsq)],GCV[length(GCV)]))
    ####print("here im")
    output <- improve.step_2d(gap,knotDist,radii,invInd,dists,gridResp,grid,explData,xvals, yvals, length(aR),response,explanatory,maxIterations,fitnessMeasure, point,knotPoint,position,aR,BIC,track,out.lm,improveNudge,tol,baseModel,radiusIndices,models, interactionTerm, data, initDisp)
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

    out.lm$splineParams[[1]]$knotPos<-aR
    out.lm$splineParams[[1]]$radiusIndices<-radiusIndices
    baseModel$splineParams<-out.lm$splineParams

    ###################################drop step#################################
    if (length(aR) > minKnots) {
      output <- drop.step_2d(radii,invInd,dists,explData,response,explanatory,maxIterations,fitnessMeasure,point,knotPoint,position,aR,BIC,track,out.lm,improveDrop,minKnots,tol,baseModel,radiusIndices,models, interactionTerm, data, initDisp)
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

      out.lm$splineParams[[1]]$knotPos<-aR
      out.lm$splineParams[[1]]$radiusIndices<-radiusIndices
      baseModel$splineParams<-out.lm$splineParams

    }
    if ((improveEx) | (improveNudge) | (improveDrop)) overallImprove = 1
  }
  ####################################write to file###############################
  ####track <- rbind(track,cbind("writing",t(aR),BIC[length(BIC)],adjRsq[length(adjRsq)],GCV[length(GCV)]))
  ####print("here fin")
  print("And we're done...")

  # if(suppress.printout){
  #   sink()
  # }
  #
    return(list(outputFS=c(length(aR),BIC[length(BIC)],aR),aR=aR,track=track, radiusIndices = radiusIndices, out.lm=out.lm,invInd=invInd,models=models,actualKnotIndices=invInd[aR],improve=overallImprove))

}
