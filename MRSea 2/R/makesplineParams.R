#' Constructing an object of spline parameters
#' 
#' This function makes a list object containing all of the information to fit splines to continuous data.  
#' 
#' @param data Data frame containing columns of covariates listed in \code{varlist}.  Column names must match with names in \code{varlist}
#' @param varlist Vector of variable names for the covariates of interest
#' @param predictionData Data frame containing columns of covariates listed in \code{varlist}.  Column names must match with those in \code{varlist}.  This parameter is used to find the maximum range of covariates between the data and prediction data. If \code{predictionData} is \code{NULL} then the range of the data is used.
#' @param degree Vector specifying the degree of the spline. If unspecified, degree 2 is stored.
#' 
#' @details 
#' The information is stored in list slots \code{[[2]]} and onward (slot \code{[[1]]} is reserved for a spatial term). Specifically: 
#' 
#' \code{covar}. Name of covariate.
#' 
#' \code{explanatory}. Vector of covariate data.
#' 
#' \code{knots}. Knot(s) for spline fitting.  This function initialises with a knot at the mean covariate value.
#' 
#' \code{bd}. This specifies the boundary knots.  If \code{predictionData} is \code{NULL} then this is the range of the covariate data.  Otherwise, the boundary knots are the maximum combined range of the data and prediction data.
#' 
#' \code{degree}. The degree of a B-spline. This function retuns 2 by default.
#' 
#' 
#' See \code{\link{runSALSA2D}} for details on the spatial slot (\code{[[1]]})
#' 
#' @examples
#' # load data
#' data(ns.data.re)
#' # load prediction data
#' data(ns.predict.data.re)
#' 
#' splineParams<- makesplineParams(ns.data.re, varlist=c('observationhour', 'DayOfMonth'),
#'                 predictionData=ns.predict.data.re)
#' 
#' @export
#' 
makesplineParams<-function(data, varlist, predictionData=NULL, degree=NULL){
  
  if(is.null(degree)){degree<-c(rep(2, length(varlist)))}
  if(is.null(predictionData)){predictionData<-data}
  # number of covariates with spline terms 
  # make object called splineParams which contains the knots, boundary knots and data needed for each one dimensional covariate
  # at the same time create all the terms needed.  
  splineParams=list(2)
  splineParams[[1]]=list() #2D Metadata
  splineParams[[2]]=list() #1D Metadata
  splineParams[[2]]$covar <- varlist[1]
  splineParams[[2]]$explanatory<-eval(parse(text=paste('data$', varlist[1], sep=''))) #explanatory
  splineParams[[2]]$knots<-mean(splineParams[[2]]$explanatory)
  splineParams[[2]]$bd<- eval(parse(text=paste("c(min(data$",varlist[1], ",predictionData$", varlist[1], "),max(data$", varlist[1], ", predictionData$", varlist[1], "))", sep="")))
  splineParams[[2]]$degree <- degree[1]
  
  for(i in 2:(length(varlist)+1)){
    splineParams[[i]]=list() #1D Metadata
    splineParams[[i]]$covar<-varlist[(i-1)]
    splineParams[[i]]$explanatory<-eval(parse(text=paste('data$', varlist[(i-1)], sep=''))) #explanatory  
    splineParams[[i]]$knots<-mean(splineParams[[i]]$explanatory)
    splineParams[[i]]$bd<- eval(parse(text=paste("c(min(data$",varlist[(i-1)], ",predictionData$", varlist[(i-1)], "),max(data$", varlist[(i-1)], ", predictionData$", varlist[(i-1)], "))", sep="")))
    splineParams[[i]]$degree <- degree[(i-1)]
  }
  
  return(splineParams)
}

