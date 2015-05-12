#' Running SALSA for a spatial smooth with a CReSS basis
#' 
#' This function fits a spatially adaptive two dimensional smooth of spatial coordinates with knot number and location selected by SALSA. 
#' 
#' @param model A model with no spatial smooth
#' @param salsa2dlist Vector of objects required for \code{runSALSA2D}: \code{fitnessMeasure}, \code{knotgrid}, \code{startKnots}, \code{minKnots}, code{maxKnots}, \code{r_seq}, \code{gap}, \code{interactionTerm}.
#' @param d2k (n x k) Matrix of distances between all data points in \code{model} and all valid knot locations specified in \code{knotgrid}
#' #' @param k2k (k x k) Matrix of distances between all valid knot locations specified in \code{knotgrid}
#' @param splineParams (default \code{=NULL}) List object containng output from runSALSA (e.g. knot locations for continuous covariates)
#' 
#' @references Scott-Hayward, L.; M. Mackenzie, C.Donovan, C.Walker and E.Ashe.  Complex Region Spatial Smoother (CReSS). Journal of computational and Graphical Statistics. 2013. doi: 10.1080/10618600.2012.762920
#' 
#' @references Scott-Hayward, L.. Novel Methods for species distribution mapping including spatial models in complex regions: Chapter 5 for SALSA2D methods. PhD Thesis, University of St Andrews. 2013
#' 
#' @details
#' 
#'There must be a column called \code{response} in the data, which is the response variable used in the initial model to be fitted.
#'
#'The object \code{salsa2dlist} contains parameters for the \code{runSALSA2D} function.  
#'
#'    \code{fitnessMeasure}. The criterion for selecting the `best' model.  Available options: AIC, AIC_c, AIC_H, BIC, QIC_b.
#'    
#'    \code{knotgrid}. A grid of legal knot locations.  Must be a regular grid with \code{c(NA, NA)} for rows with an illegal knot.  An illegal knot position may be outside the study region or on land for a marine species for example.
#'    
#'    \code{knotdim}. The dimensions of the knot grid as a vector. (x, y)
#'    
#'    \code{startknots}. Starting number of knots (initialised as spaced filled locations).
#'    
#'    \code{minKnots}. Minimum number of knots to be tried.
#'    
#'    \code{maxKnots}. Maximum number of knots to be tried.
#'    
#'    \code{r_seq}. Sequence of range parameters for the CReSS basis from local (small) to global (large).  Determines the range of the influence of each knot. Sequence made using \code{\link{getRadiiChoices}}.
#'    
#'    \code{gap}. The minimum gap between knots (in unit of measurement of coordinates).
#'    \code{interactionTerm}. Specifies which term in \code{baseModel} the spatial smooth will interact with.  If \code{NULL} no interaction term is fitted.
#'    
#' 
#' @return
#' The spline paramater object that is return now contains a list in the first element (previously reserved for the spatial component).  This list contains the objects required for the SALSA2D fitting process:
#' 
#' \item{knotDist}{Matrix of knot to knot distances (k x k).  May be Euclidean or geodesic distances. Must be square and the same dimensions as \code{nrows(na.omit(knotgrid))}.  Created using \code{\link{makeDists}}.}
#' \item{radii}{Sequence of range parameters for the CReSS basis from local (small) to global (large).  Determines the range of the influence of each knot.}
#' \item{dist}{ Matrix of distances between data locations and knot locations (n x k). May be Euclidean or geodesic distances. Euclidean distances created using \code{\link{makeDists}}.}
#' \item{gridresp}{The first column of knotgrid.}
#' \item{grid}{Index of knotgrid locations.  Should be same length as \code{knotgrid} but with x=integer values from 1 to number of unique x-locations and y= integer values from 1 to number of unique y-locations.}
#' \item{datacoords}{Coordinates of the data locations}
#' \item{response}{Vector of response data for the modelling process}
#' \item{knotgrid}{Grid of legal knot locations.  Must be a regular grid with c(NA, NA) for rows with an illegal knot.}
#' \item{minKnots}{Minimum number of knots to be tried.}
#' \item{maxKnots}{Maximum number of knots to be tried.}
#' \item{gap}{Minimum gap between knots (in unit of measurement of \code{datacoords})}
#' \item{radiusIndices}{Vector of length startKnots identifying which radii (\code{splineParams[[1]]$radii}) will be used for each knot location (\code{splineParams[[1]]$knotPos})}
#' \item{knotPos}{Index of knot locations. The index identifies which knots (i.e. which rows) from \code{knotgrid} were selected by SALSA}
#' \item{invInd}{This is a vector of length the number of rows of \code{knotgrid}.  It is used to translate between \code{knotgrid} (used in SALSA) and \code{na.omit(knotgrid)} (used in \code{dist} and \code{LocalRadialFunction}).}
#' 
#' 
#' @examples
#' # load data
#' data(ns.data.re)
#' # load prediction data
#' data(ns.predict.data.re)
#' # load knot grid data
#' data(knotgrid.ns)
#' 
#' splineParams<-makesplineParams(data=ns.data.re, varlist=c('observationhour'))
#' 
#' #set some input info for SALSA
#' ns.data.re$response<- ns.data.re$birds
#' 
#' # make distance matrices for datatoknots and knottoknots
#' distMats<-makeDists(cbind(ns.data.re$x.pos, ns.data.re$y.pos), na.omit(knotgrid.ns))
#' 
#' # choose sequence of radii
#' r_seq<-getRadiiChoices(8, distMats$dataDist)
#' 
#' # set initial model without the spatial term
#' # (so all other non-spline terms)
#' initialModel<- glm(response ~ as.factor(floodebb) + as.factor(impact) + offset(log(area)), 
#'                    family='quasipoisson', data=ns.data.re)
#' 
#' # make parameter set for running salsa2d
#' salsa2dlist<-list(fitnessMeasure = 'QICb', knotgrid = knotgrid.ns, knotdim = c(7, 9), 
#'                   startKnots=6, minKnots=4, maxKnots=20, r_seq=r_seq, gap=1, 
#'                   interactionTerm="as.factor(impact)")
#' 
#' salsa2dOutput_k6<-runSALSA2D(initialModel, salsa2dlist, d2k=distMats$dataDist, 
#'                             k2k=distMats$knotDist, splineParams=splineParams) 
#'
#'@author Cameron Walker
#'
#' @export
#' 
runSALSA2D<-function(model, salsa2dlist, d2k, k2k, splineParams=NULL, tol=0){
  
  if(class(model)[1]=='glm'){
    data<-model$data  
  }
  if(class(model)[1]=='lm'){
    data<-model$model
  }
  
  attributes(model$formula)$.Environment<-environment()
  
  # check for response variable
  if(is.null(data$response)) stop('data does not contain response column')
  
  #set input 2D input data
  if(is.null(data$x.pos)){ stop('no x.pos in data frame; rename coordinates')}
  explData<- cbind(data$x.pos, data$y.pos)
  #specify the full grid (including NAs where knots cannot go)
  explanatory<- salsa2dlist$knotgrid
  #this sets the index of grid positions
  grid<-expand.grid(1:salsa2dlist$knotdim[1], 1:salsa2dlist$knotdim[2])
  gridResp<-salsa2dlist$knotgrid[,1]
  
  if(length(salsa2dlist$r_seq)>1){
    radii<- salsa2dlist$r_seq[round(length(salsa2dlist$r_seq)/2)]  
  }else{
    radii<- salsa2dlist$r_seq[1]
  }
  winHalfWidth = 0
  
  interactionTerm<-salsa2dlist$interactionTerm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~ SET UP ~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #number of covariates with spline terms  
  if(is.null(splineParams)){
    splineParams<-list(2)    
  }
  
  splineParams[[1]]=list(13) #2D Metadata
  splineParams[[1]][[1]]= NULL                        #starting knot positions
  splineParams[[1]][['knotDist']]= k2k
  splineParams[[1]][['dist']]= d2k
  splineParams[[1]][[4]]= NULL                        #starting radius indices
  splineParams[[1]][['gridResp']]= gridResp
  splineParams[[1]][['grid']]= grid
  splineParams[[1]][['response']]= data$response
  splineParams[[1]][['knotgrid']]= salsa2dlist$knotgrid
  splineParams[[1]][['datacoords']]= explData
  splineParams[[1]][['radii']]= radii
  splineParams[[1]][['minKnots']]= salsa2dlist$minKnots
  splineParams[[1]][['maxKnots']]= salsa2dlist$maxKnots
  splineParams[[1]][['gap']]= salsa2dlist$gap
  splineParams[[1]][[14]]= NULL                       #this will be invInd
  
  if(dim(k2k)[2]==salsa2dlist$minKnots & dim(k2k)[2]==salsa2dlist$maxKnots & dim(k2k)[2]==salsa2dlist$startKnots){stop('Min, Max and Start knots all identical and equal to the total number of valid knots\n Please add more valid knot locations (knotgrid) or reduce min/max/start')}
  
  if(dim(k2k)[2]<salsa2dlist$minKnots | dim(k2k)[2]<salsa2dlist$startKnots){stop('Starting number of knots or min knots more than number of valid knot locations in knotgrid')}
  
  if(salsa2dlist$fitnessMeasure=='BIC' & initialModel$family$family=='quasipoisson' ){stop('Please use fitness Measure appropriate for quasi family, e.g. QICb')}
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~ 2D SALSA RUN ~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  baseModel1D<- model
  baseModel<- baseModel1D
  
  output<-return.reg.spline.fit.2d(splineParams, startKnots=salsa2dlist$startKnots, winHalfWidth,fitnessMeasure=salsa2dlist$fitnessMeasure, maxIterations=10, tol=tol, baseModel=baseModel, radiusIndices=NULL, initialise=TRUE,  initialKnots=NULL, interactionTerm=interactionTerm, knot.seed=10)
  
  baseModel<- output$out.lm
    
  
  if(length(output$models)>0){
    modRes<- c()
    for(m in 1:length(output$models)){
      modelNo<- m
      knotPosition<-output$models[[m]][[1]]
      rIs<- output$models[[m]][[2]]
      r<- output$models[[m]][[3]]
      fitScore<- output$models[[m]][[4]]
      modRes<- rbind(modRes, data.frame(modelNo, knotPosition, rIs, fitScore))
    }
    
    modRes<-modRes[order(modRes$fitScore),]
    bestModNo<- unique(modRes$modelNo[which(modRes$fitScore==min(modRes$fitScore))])[1]
    
    a<-sum(as.vector(output$invInd[output$aR]) - as.vector(unlist(output$models[[bestModNo]][1])))
    print(paste('a = ', a, sep=''))
    if(a!=0) break 
    
    #output$aR<- output$models[[bestModNo]][[1]]
  }
  
#   if(length(output$models)==0){
#     output$aR<- output$invInd[output$aR]  
#   }
  
  # use initialise step to change radii
  radii = salsa2dlist$r_seq
  x <- as.vector(splineParams[[1]]$grid[,1])
  y <- as.vector(splineParams[[1]]$grid[,2])
  xvals <- max(x)
  yvals <- max(y)
  
  if(length(salsa2dlist$r_seq)>1){
    radiusIndices<- rep(round(length(salsa2dlist$r_seq)/2), (length(output$aR)))
  }else{
      radiusIndices<- rep(1, length(output$aR))
  }
  output_radii<- initialise.measures_2d(k2k, maxIterations=10, salsa2dlist$gap, radii, d2k, gridResp, explData, splineParams[[1]]$startKnots, xvals, yvals, explanatory, splineParams[[1]]$response, baseModel, radiusIndices=radiusIndices, initialise=F, initialKnots=salsa2dlist$knotgrid[output$aR,], fitnessMeasure=salsa2dlist$fitnessMeasure, interactionTerm=interactionTerm, data=data, knot.seed=10)
  
  
  splineParams[[1]]$radii= radii
  splineParams[[1]][['knotPos']]= output_radii$aR
  TwoDModelsfirstStage <- output_radii$models
  splineParams[[1]][['radiusIndices']]= output_radii$radiusIndices
  splineParams[[1]][['invInd']] = output_radii$invInd
  modelFit = output_radii$BIC                                         
  
  baseModel<-output_radii$out.lm
  attributes(baseModel$formula)$.Environment<-globalenv()
  
  #save.image(paste("salsa2D_k", splineParams[[1]]$startKnots, ".RData", sep=''))
  return(list(bestModel=baseModel, splineParams=splineParams, fitStat=output_radii$BIC, aR=list(aR1=output$aR, aR2=output_radii$aR)))
}