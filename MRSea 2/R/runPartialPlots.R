
#-----------------------------------------------------------------------------
#' Plot partial plots for each of the variables listed in \code{factorlist} or \code{varlist}.
#' 
#' @param model Fitted model object (glm or gam)
#' @param data Data frame of data information used to fit \code{model}
#' @param factorlist (\code{default=NULL}). Vector or names of factor variables
#' @param varlist (\code{default=NULL}). Vector of names of continuous variables
#' @param showKnots (\code{default=FALSE}). Logical stating whether knot locations should be plotted.
#' @param save (\code{default=FALSE}). Logical stating whether plot should be saved into working directory.
#' @param savedata (\code{default=FALSE}). Logical stating whether the data used to make the plots should be saved into the working directory.  The object is called PartialData_'variablename'.RData
#' 
#' @return
#' Partial plots, one for each covariate in \code{factorlist} and \code{varlist}
#' 
#' @examples
#' 
#' #' # load data
#' data(ns.data.re)
#' 
#' model<-glm(birds ~ observationhour + as.factor(floodebb) + as.factor(impact), 
#'            family='quasipoisson', data=ns.data.re)
#' 
#' runPartialPlots(model, ns.data.re, factorlist=c('floodebb', 'impact'), 
#'                 varlist=c('observationhour'))
#'
#' @export
#' 

runPartialPlots<-function(model, data, factorlist=NULL, varlist=NULL, showKnots=FALSE, save=FALSE, savedata=F){
  
  print("Making partial plots")
  
  require(mvtnorm)
  
  if(save==T){png("PartialFitsLink%i.png", height=600, width=650)}else{devAskNewPage(ask=TRUE)}
  par(mfrow=c(1,1))
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(is.null(factorlist)==F){
    for(i in 1:length(factorlist)){
      coeffac<- c(grep(factorlist[i], colnames(model.matrix(model))))
      coefradial<-c(grep('LocalRadialFunction', colnames(model.matrix(model))))
      coefpos <- coeffac[which(is.na(match(coeffac, coefradial)))]
      xvals<-data[,which(names(data)==factorlist[i])] 
      newX<- sort(unique(xvals))
      newX<- newX[2:length(newX)]
      #eval(parse(text=paste(varlist[i], "<- newX", sep="")))
      #response<- rep(1, 500)
      #newBasis<- eval(parse(text=labels(terms(model))[grep(varlist[i], labels(terms(model)))]))
      partialfit<- coef(model)[coefpos]
      rcoefs<- NULL
      try(rcoefs<- rmvnorm(1000,coef(model), summary(model)$cov.scaled), silent=T)
      if(is.null(rcoefs) || length(which(is.na(rcoefs)==T))>0){
        rcoefs<- rmvnorm(1000,coef(model), as.matrix(nearPD(summary(model)$cov.scaled)$mat))}
      
      rpreds<- rcoefs[,coefpos]
      quant.func<- function(x){quantile(x, probs=c(0.025, 0.975))}  
      # factor variables have rpreds vector
      if(is.null(dim(rpreds))){
        cis<-quant.func(rpreds)
        plot(newX, partialfit, pch=20,
             xlab=factorlist[i], ylab="Partial Fit", lwd=2, xaxt="n", ylim=c(range(0,cis)), cex.lab=1.3, cex.axis=1.3)
        axis(1, at=newX, labels=colnames(model.matrix(model))[coefpos])
        segments(newX,cis[1], newX , cis[2], lwd=2)
        abline(h=0, lwd=2)
      }else{
        # cts variables have matrix of rpreds
        cis<- t(apply(rpreds, 2,quant.func))
        plot(1:length(newX), partialfit, pch=20,
             xlab=factorlist[i], ylab="Partial Fit", lwd=2, xaxt="n", ylim=c(range(0,cis)), cex.lab=1.3, cex.axis=1.3)
        axis(1, at=(1:length(newX)), labels=newX)
        segments(1:length(newX),cis[,1], 1:length(newX) , cis[,2], lwd=2)
        abline(h=0, lwd=2)
      }
      
      if(savedata==T){
        partialdata<-data.frame(newX, partialfit, cis[1], cis[2])
      save(partialdata, file=paste('PartialData_', factorlist[i], '.RData', sep=''), compress='bzip2' )  
      }
      
    }
  }
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(is.null(varlist)==F){
    n<-length(varlist)
    for(i in 1:n){
      coefpos<- c(1,grep(varlist[i], colnames(model.matrix(model))))
      xvals<-data[,which(names(data)==varlist[i])] 
      newX<- seq(min(xvals),max(xvals), length=500)
      eval(parse(text=paste(varlist[i], "<- newX", sep="")))
      response<- rep(1, 500)
      
      newBasis<- eval(parse(text=labels(terms(model))[grep(varlist[i], labels(terms(model)))]))
      partialfit<- cbind(rep(1,500),newBasis)%*%coef(model)[coefpos]
      rcoefs<- NULL
      try(rcoefs<- rmvnorm(1000,coef(model), summary(model)$cov.scaled), silent=T)
      if(is.null(rcoefs) || length(which(is.na(rcoefs)==T))>0){
        rcoefs<- rmvnorm(1000,coef(model), as.matrix(nearPD(summary(model)$cov.scaled)$mat))}
      
      rpreds<- cbind(rep(1,500),newBasis)%*%t(rcoefs[,coefpos])
      quant.func<- function(x){quantile(x, probs=c(0.025, 0.975))}  
      cis<- t(apply(rpreds, 1,quant.func))
      
      plot(newX, partialfit, type="l",
           xlab=varlist[i], ylab="Partial Fit", lwd=2, ylim=range(cis), cex.lab=1.3, cex.axis=1.3)
      rug(xvals)
      lines(newX,cis[,1], col="darkred", lty=4)
      lines(newX,cis[,2], col="darkred", lty=4)
      if(showKnots=="TRUE"){
        abline(v=splineParams[[(i+1)]]$knots, lty=4, col="grey")}
      eval(parse(text=paste("rm(", varlist[i], ")", sep="")))
      rm(response)
      
      if(savedata==T){
        partialdata<-data.frame(newX, partialfit,  cis[,1], cis[,2])
        save(partialdata, file=paste('PartialData_', varlist[i], '.RData', sep=''), compress='bzip2' )  
      }
    }
  }
  if(save==T){dev.off()}else{devAskNewPage(ask=FALSE)}
  
  
  
}
  
