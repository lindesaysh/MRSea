#-----------------------------------------------------------------------------
#' Calculate cumulative residuals and plot.  
#' 
#' The output is plots of cumulative residuals.  
#' 
#' @param model Fitted model object (glm or gam)
#' @param varlist Vector of covariate names (continous covariates only)
#' @param d2k (default=NULL).  Distance matrix of data to knot points. Used only if there is a \code{\link{LocalRadialFunction}} smooth in the model formula
#' @param splineParams (default \code{=NULL}) List object containing output from runSALSA/runSALSA2D required for updating \code{model}.  Used only if there is a \code{LocalRadialFunction} smooth in the model formula. See \code{\link{makesplineParams}} for details of this object.
#' @param label Label printed at the end of the plot name to identify it if \code{save=TRUE}.
#' @param save (\code{default=FALSE}). Logical stating whether plot should be saved into working directory.
#'
#' @return
#' Cumulative residual plots are returned for residuals ordered by each covariate in \code{varlist}, predicted value and index of observations (temporally).
#' The blue dots are the residuals
#' The black line is the line of cumulative residual.
#' On the covariate plots (those in \code{varlist}) the grey line indicates what we would expect from a well fitted covariate. i.e. one that is fitted with excessive knots.
#' 
#' Note: if the covariate is discrete in nature (like the example below), there will be a lot of overplotting of residuals.
#'
#' @examples 
#' # load data
#' data(ns.data.re)
#' 
#' model<-gamMRSea(birds ~ observationhour + as.factor(floodebb) + as.factor(impact), 
#'            family='quasipoisson', data=ns.data.re)
#' 
#' plotCumRes(model, varlist=c('observationhour'))
#' 
#' @export
#' 
plotCumRes<- function(model, varlist, d2k=NULL, splineParams=NULL, label='', save=FALSE){
  
  require(splines)
  
  print("Calculating cumulative residuals")
  
  namesOfx<- c(varlist, c("Predicted", "Index"))
  
  if(class(model)[1]=='geeglm' | class(model)[1]=='glm'){
    dat<- model$data
  }
  
  if(class(model)[1]=='gam'){
    dat<-model$model 
  }
  
  if(class(model)[1]=='gamMRSea'){
    dat<-model$data 
  }
  
  # find which data columns refer to the used variables
  coefpos<-c()
  for(c in 1:(length(namesOfx)-2)){
    coefpos<- c(coefpos, grep(namesOfx[c], names(dat)))
  }
  
  # make matrix (n x number covar + fitted + orderdata)
  xmatrix<- cbind(dat[,coefpos], fitted(model), ord=1:length(fitted(model)))
  
  attributes(model$formula)$.Environment<-environment()
  
  for(c in 1:length(namesOfx)){
    
    if(c < (length(namesOfx)-1)){
      covardat<- dat[model$y>0, coefpos[c]]
      newknots<-as.vector(quantile(covardat, seq(0.05, 0.95, length=7)))
      newknots<-unique(newknots)
      dists<-d2k
      # removal term
      #term<-attributes(terms(model))$term.labels[grep(namesOfx[i], attributes(terms(model))$term.labels)]
      term<-labels(terms(model))[grep(namesOfx[c], labels(terms(model)))]
      newterm<- paste('bs(', namesOfx[c], ', knots=newknots)', sep='')
      eval(parse(text=paste('covarModelUpdate<-update(model, .~. -',term, ' +', newterm,')', sep='')))
      xmatrix_test<- cbind(dat[,coefpos], fitted(covarModelUpdate), ord=1:length(fitted(covarModelUpdate)))
    }
    
    type='response'
    
    maxy<- max(residuals(model),cumsum(residuals(model, type=type)[order(xmatrix[,c])]))
    miny<- min(residuals(model),cumsum(residuals(model, type=type)[order(xmatrix[,c])]))
    
    # make plot
    if(save==T){png(paste("CumRes_", namesOfx[c],label, ".png", sep=''), height=600, width=700)}
    else{devAskNewPage(ask=TRUE)}
    plot(xmatrix[,c][order(xmatrix[,c])],residuals(model, type=type)[order(xmatrix[,c])], ylim=c(miny,maxy), xlab=namesOfx[c], ylab="Response   residuals", pch='.',cex=0.1, main="Cumulative Residuals", cex.lab=1.3, cex.axis=1.3)
    
    if(c < (length(namesOfx)-1)){
      lines(xmatrix_test[,c][order(xmatrix_test[,c])],cumsum(residuals(covarModelUpdate, type=type)[order(xmatrix_test[,c])]), lwd=2, col='grey')  
    }
    points(xmatrix[,c][order(xmatrix[,c])],residuals(model, type=type)[order(xmatrix[,c])], pch=20, col='turquoise4', cex=0.5)
    lines(xmatrix[,c][order(xmatrix[,c])],cumsum(residuals(model, type=type)[order(xmatrix[,c])]), lwd=2)
    abline(h=0)
    if(save==T){dev.off()}else{devAskNewPage(ask=FALSE)}
  }
}