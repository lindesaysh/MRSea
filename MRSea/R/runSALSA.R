#' Running SALSA for continuous one-dimensional covariates.
#' 
#' This function finds spatially adaptive knot locations for one or more continuous one-dimensional covariates.
#' 
#' @param initialModel The best fitting \code{CReSS} model with no continuous covariates specified
#' @param salsa1dlist Vector of objects required for \code{runSALSA1D}: \code{fitnessMeasure}, \code{minKnots_1d}, \code{maxKnots_1d}, \code{startKnots_1d} \code{degree}, \code{maxIterations} \code{gap}. 
#' @param varlist Vector of variable names for the covariates required for knot selection
#' @param factorlist vector of factor variables specified in \code{initialModel}.  Specified so that a check can be made that there are non-zero counts in all levels of each factor. Uses the function \code{checkfactorlevelcounts}. Default setting is NULL.
#' @param predictionData The data to be predicted for. column names correspond to the data in \code{initialModel}
#' @param varlist_cyclicSplines Vector of variable names for covariates to be modelled with cyclic cubic splines.  This must be a subset of \code{varlist}.The default is \code{NULL}
#' @param splineParams List object containing information for fitting splines to the covariates in \code{varlist}. If not specified (\code{NULL}) this object is created and returned. See \code{\link{makesplineParams}} for details.
#' 
#' @details
#' There must be a column called \code{response} in the data, which is the response variable used in the initial model to be fitted.
#' 
#' The object \code{salsa1dlist} contains parameters for the \code{runSALSA1D} function.
#'   
#' \code{fitnessMeasure}. The criterion for selecting the `best' model.  Available options: AIC, AIC_c, BIC, QIC_b.
#' 
#'    \code{minKnots_1d}. Minimum number of knots to be tried.
#'    
#'    \code{maxKnots_1d}. Maximum number of knots to be tried.
#'    
#'    \code{startKnots_1d}. Starting number of knots (spaced at quantiles of the data).
#'    
#'    \code{degree}. The degree of the B-spline. Does not need to be specified if \code{splineParams} is a parameter in \code{runSALSA1D}.
#'    
#'    \code{maxIterations}.The exchange/improve steps will terminate after maxIterations if still running.
#'    
#'    \code{gaps}. The minimum gap between knots (in unit of measurement of explanatory).
#'    
#'    
#' \code{minKnots_1d}, \code{maxKnots_1d}, \code{startKnots_1d} and \code{gaps} are vectors the same length as \code{varlist}.  This enables differing values of these parameters for each covariate.
#'
#'  The initial model contains all the factor level covariates and any covariates of interest that are not specified in the \code{varlist} argument of \code{runSALSA1D} 
#' 
#' \emph{Note:} The algorithm will not remove variables in \code{varlist}.  If there is no better model than with a knot at the mean, the output will include that covariate with a knot at the mean.  The user must decide if the covariate is required in the model as a linear term instead.
#' 
#' @return
#' A list object is returned containing 4 elements:
#' 
#' \item{bestModel}{A glm model object from the best model fitted}
#' \item{modelFits1D}{A list object with an element for each new term fitted to the model.  The first element is a model fitted with a knot at the mean for each of the covariates in \code{varlist}.  Within the first element, the model term of interest, the current fit, knots and formula.  The second element is the result of SALSA on the first term in \code{varlist}.  Within this element, the knots chosen and the improvement in model fit are presented \code{$modelfits}.  This continues till all covariates in \code{varlist} have been through SALSA.}
#' \item{splineParams}{The updated spline parameter object, with the new (if chosen) knot locations for each covariate in \code{varlist}}
#' \item{fitstat}{The final fit statistic of \code{bestModel}.  The type of statistic was specified in \code{salsa1dlist}.}
#' 
#' 
#' @references Walker, C.; M. Mackenzie, C. Donovan and M. O'Sullivan. SALSA - a Spatially Adaptive Local Smoothing Algorithm. Journal of Statistical Computation and Simulation, 81(2):179-191, 2010
#' 
#' @examples
#' # load data
#' data(ns.data.re)
#' # load prediction data
#' data(ns.predict.data.re)
#' 
#' splineParams<-makesplineParams(data=ns.data.re, varlist=c('observationhour', 'DayOfMonth'))
#' #set some input info for SALSA
#' ns.data.re$response<- ns.data.re$birds
#' 
#' #' # set initial model without the spline terms in there 
#' # (so all other non-spline terms)
#' initialModel<- glm(response ~ as.factor(floodebb) + as.factor(impact) + offset(log(area)), 
#'                     family='quasipoisson',data=ns.data.re)
#'
#' salsa1dlist<-list(fitnessMeasure = 'QICb', minKnots_1d=c(2,2), maxKnots_1d = c(20, 20), 
#'                   startKnots_1d = c(2,2), degree=c(2,2), maxIterations = 10, gaps=c(1,1))
#' # run SALSA
#' salsa1dOutput<-runSALSA1D(initialModel, salsa1dlist, varlist=c('observationhour', 'DayOfMonth'), 
#'                  factorlist=c('floodebb', 'impact'), ns.predict.data.re, splineParams=splineParams)
#' 
#' @export
#' 
runSALSA1D<-function(initialModel, salsa1dlist, varlist, factorlist=NULL, predictionData, varlist_cyclicSplines=NULL, splineParams=NULL){
  
  if(is.null(factorlist)==F){
    # check factor level counts:
  checkfactorlevelcounts(factorlist, initialModel$data, initialModel$y)
  }
  
  # ~~~~~~~~~~~~ SET UP ~~~~~~~~~~~~~~~~
  # set parameters for SALSA
  winHalfWidth = 0
  family=initialModel$family$family
  data<-initialModel$data
  
  attributes(initialModel$formula)$.Environment<-environment()
  
  # check for response variable
  if(is.null(data$response)) stop('data does not contain response column')
  
  # check parameters in salsa1dlist are same length as varlist
  if(length(varlist)!=length(salsa1dlist$minKnots_1d)) stop('salsa1dlist$minKnots_1d not same length as varlist')
  if(length(varlist)!=length(salsa1dlist$maxKnots_1d)) stop('salsa1dlist$maxKnots_1d not same length as varlist')
  if(length(varlist)!=length(salsa1dlist$startKnots_1d)) stop('salsa1dlist$startKnots_1d not same length as varlist')
  if(length(varlist)!=length(salsa1dlist$degree)) stop('salsa1dlist$degree not same length as varlist')
  if(length(varlist)!=length(salsa1dlist$gaps)) stop('salsa1dlist$gaps not same length as varlist')
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(splineParams)){
    splineParams<-makesplineParams(data, varlist, predictionData, salsa1dlist$degree)
    varID<-(1:length(varlist))+1
  }else{
    # check whats in varlist and splineParams to see if they match
    if((length(splineParams)-1)!=length(varlist)){
      varID<-vector(length=length(varlist))
      for(i in 1:length(varlist)){
        varID[i]<-grep(varlist[i], splineParams)
      }
    }else{varID<-(1:length(varlist))+1}
  }
  
  terms1D <- list(length(varlist))
  
  for(i in 2:(length(varlist)+1))
    if(varlist[varID[(i-1)]]%in%varlist_cyclicSplines){
      require(mgcv)
      terms1D[[(i-1)]]<- paste("as.matrix(data.frame(gam(response ~ s(", varlist[varID[(i-1)]], ", bs='cc', k=(length(splineParams[[", varID[(i-1)], "]]$knots) +2)), knots = list(c(splineParams[[",varID[(i-1)], "]]$bd[1], splineParams[[", varID[(i-1)], "]]$knots, splineParams[[",varID[(i-1)], "]]$bd[2])),fit=F)$X[,-1]))", sep="")  
    }else{
      terms1D[[(i-1)]]<- paste("bs(", varlist[(i-1)], ", knots = splineParams[[", varID[(i-1)], "]]$knots, degree=splineParams[[", varID[(i-1)], "]]$degree, Boundary.knots=splineParams[[",varID[(i-1)], "]]$bd)", sep='')
    }
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # initial model is year one knot at mean and NO 2D smooth
  require(mgcv)
  baseModel <- eval(parse(text=paste("glm(response ~ ", paste(formula(initialModel)[3],sep=""), "+", paste(terms1D, collapse="+"),", family =", family, ", data = data)", sep='')))
  
  modelFit = BIC(baseModel)
  fitStat<-get.measure(salsa1dlist$fitnessMeasure,'NA',baseModel)$fitStat
  
  modelFits1D <- list((length(varlist)+1))
  modelFits1D[[1]] <- list(term = terms1D[[1]], tempfits = c(fitoutput = fitStat, get.measure(salsa1dlist$fitnessMeasure,'NA',baseModel)$fitStat), knots = splineParams[[2]]$knots, formula = baseModel$call)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~ loop through 1D covar ~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  starttime = proc.time()
  timings<- vector(length=length(varlist))
  knots=NULL
  #baseModel <- eval(parse(text=paste("update(baseModel, ~. -", term, ")", sep="")))
  for (i in 2:(length(varlist)+1)){
    explanatory <- splineParams[[varID[(i-1)]]]$explanatory
    response<-data$response
    bd <- as.numeric(splineParams[[varID[(i-1)]]]$bd)   # i is the location of covar in varid +1 (2d has 1st entry in spline params)
    gap <- (salsa1dlist$gaps[(i-1)])
    term<- terms1D[[(i-1)]]
    baseModel <- eval(parse(text=paste("update(baseModel, ~. -", term, ")", sep="")))
    if(length(grep(varlist[(i-1)], baseModel$formula))>0){stop(paste('Multiple instances of covariate in model. Remove ',splineParams[[varID[(i-1)]]]$covar , ' before proceding', sep=''))}
    
    if(varlist[(i-1)]%in%varlist_cyclicSplines){spl<- "cc"}else{spl="bs"}
    
    if(spl == "cc"){minKnots_1d <- 3; startKnots_1d<-3}
    sttime<- proc.time()[3]
    output <- return.reg.spline.fit(response,explanatory,splineParams[[varID[(i-1)]]]$degree,salsa1dlist$minKnots_1d[(i-1)],salsa1dlist$maxKnots_1d[(i-1)],salsa1dlist$startKnots_1d[(i-1)], gap, winHalfWidth, salsa1dlist$fitnessMeasure, maxIterations=100, baseModel=baseModel, bd=bd, spl=spl)
    
    timings[(i-1)]<- proc.time()[3] - sttime
    
    thisFit = output$output[2] 
    if (thisFit < fitStat) {
      splineParams[[varID[(i-1)]]]$knots= sort(output$aR)
      fitStat = thisFit
    }
    # update best model to have new knot locations and covariate back in model
    baseModel<- eval(parse(text=paste("update(baseModel, ~. +", term, ")", sep="")))
    # calculate a BIC score here too
    thisBIC = BIC(baseModel)
    
    modelFits1D[[i]] <- list(term = term, tempfits = c(fitoutput = thisFit, get.measure(salsa1dlist$fitnessMeasure,'NA',output$out.lm)$fitStat), knots = sort(c(output$aR)), formula = baseModel$call, modelfits = c(fit = thisFit))
  }
  
  attributes(baseModel$formula)$.Environment<-.GlobalEnv
  #save.image("Test.RData")
  return(list(bestModel=baseModel, modelFits1D=modelFits1D, splineParams=splineParams, fitStat=fitStat))
  
}


