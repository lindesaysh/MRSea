#' Plot partial plots for each of the variables listed in \code{factorlist.in} or \code{varlist.in}.
#'
#' @param model Fitted model object (glm or gam)
#' @param data Data frame of data information used to fit \code{model}
#' @param factorlist.in (\code{default=NULL}). Vector or names of factor variables
#' @param varlist.in (\code{default=NULL}). Vector of names of continuous variables
#' @param showKnots (\code{default=FALSE}). Logical stating whether knot locations should be plotted.
#' @param type (\code{default='responss'}).  Character stating whether to return partial plots on the scale of the link function or the response.
#' @param partial.resid (\code{default=FALSE}).  Logical stating whether to include partial residuals on the plot.
#' @param save (\code{default=FALSE}). Logical stating whether plot should be saved into working directory.
#' @param savedata (\code{default=FALSE}). Logical stating whether the data used to make the plots should be saved into the working directory.  The object is called PartialData_'variablename'.RData
#' @param label (\code{default=NULL}).  This enables the user to specify a character label for the plots saved to the working directory. This may also be used to specify an alternative directory.
#' @param includeB0 (\code{default=TRUE}).  Logical stating whether to include the intercept in the partial plots. 
#'
#' @return
#' Partial plots, one for each covariate in \code{factorlist.in} and \code{varlist.in}
#'
#' @examples
#'
#' #' # load data
#' data(ns.data.re)
#'
#' model<-gamMRSea(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),
#'            family='quasipoisson', data=ns.data.re)
#'
#' runPartialPlots(model, ns.data.re, factorlist.in=c('floodebb', 'impact'),
#'                 varlist.in=c('observationhour'))
#' runPartialPlots(model, ns.data.re, factorlist.in=c('floodebb', 'impact'),
#'                 varlist.in=c('observationhour'), type='link')
#' runPartialPlots(model, ns.data.re, factorlist.in=c('floodebb', 'impact'),
#'                 varlist.in=c('observationhour'), partial.resid=TRUE)

#' @author LAS Scott-Hayward, University of St Andrews
#' 
#' @export
#'

runPartialPlots<-function(model, data, factorlist.in=NULL, varlist.in=NULL, showKnots=FALSE, type='response', partial.resid=FALSE, save=FALSE, savedata=F, label=NULL, includeB0=FALSE){

  print("Making partial plots")

  require(mvtnorm)
  require(splines)
  require(Matrix)
  
 splineParams <- model$splineParams
 
 setClass("vglmMRSea", contains=c("vglm"), slots=c(varshortnames="character", panels="ANY", splineParams="list", data="data.frame", cvfolds="numeric", interactionterm="ANY")) -> vglmMRSea
  
  if(save==T){
    if(type=='response'){
      png(paste(label,"PartialFitsResponse%i.png",sep=''), height=600, width=650)}else{
      png(paste(label,"PartialFitsLink%i.png",sep=''), height=600, width=650)
  }}else{devAskNewPage(ask=TRUE)}

  par(mfrow=c(1,1))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(factorlist.in)==F){
      for(i in 1:length(factorlist.in)){
        coeffac<- c(grep(factorlist.in[i], colnames(model.matrix(model))))
        coefradial<-c(grep('LRF', colnames(model.matrix(model))))
        coefpos <- coeffac[which(is.na(match(coeffac, coefradial)))]
        xvals<-data[,which(names(data)==factorlist.in[i])]
        newX<- sort(unique(xvals))
        newX<- newX[2:length(newX)]
        if(class(newX)=='character'){
        newX<-as.factor(newX)  
        }

        partialfit<- coef(model)[c(coefpos)]
        rcoefs<- NULL
        try(rcoefs<- rmvnorm(1000,coef(model), summary(model)$cov.robust), silent=T)
        if(is.null(rcoefs) || length(which(is.na(rcoefs)==T))>0){
          rcoefs<- rmvnorm(1000,coef(model), as.matrix(nearPD(summary(model)$cov.scaled)$mat))}

        rpreds<- rcoefs[,c(coefpos)]
        quant.func<- function(x){quantile(x, probs=c(0.025, 0.975))}
        # factor variables have rpreds vector
        if(is.null(dim(rpreds))){
          cis<-quant.func(rpreds)
          if(type=='response'){
            partialfit <- model$family$linkinv(partialfit)
            cis <- model$family$linkinv(cis)
          }
          if(partial.resid==T){
            par.res <- residuals(model, type='partial')
            par.res <- par.res[, c(grep(factorlist.in[i], colnames(par.res)))]
            if(type=='response'){
              par.res<-model$family$linkinv(par.res)
              y.lab = "Partial Fit (response)"
            }else{y.lab="Partial Fit (link)"}
            plot(xvals, par.res, xlab=factorlist.in[i], ylab=y.lab, ylim=range(c(0, cis, par.res)), xaxt="n", cex.lab=1.3, cex.axis=1.3, col='grey', cex=0.5)
            points(newX, partialfit, pch=20)
          }else{
            if(type=='response'){
              y.lab = "Partial Fit (response)"
            }else{y.lab="Partial Fit (link)"}
          plot(newX, partialfit,xlab=factorlist.in[i], ylab=y.lab,  xaxt="n", ylim=c(range(0,cis)), cex.lab=1.3, cex.axis=1.3)
          points(1:length(newX), partialfit, pch=20)
          }
          axis(1, at=newX, labels=colnames(model.matrix(model))[c(coefpos)])
          segments(1:length(newX),cis[1], 1:length(newX) , cis[2], lwd=2)
          abline(h=0, lwd=2)
        }else{
          # cts variables have matrix of rpreds
          cis<- t(apply(rpreds, 2,quant.func))
          if(type=='response'){
            partialfit <- model$family$linkinv(partialfit)
            cis <- model$family$linkinv(cis)
          }
          if(partial.resid==T){
            par.res <- residuals(model, type='partial')
            par.res<- par.res[, c(grep(factorlist.in[i], colnames(par.res)))]
            if(type=='response'){
              par.res<-model$family$linkinv(par.res)
              y.lab = "Partial Fit (response)"
            }else{y.lab="Partial Fit (link)"}
            plot(xvals, par.res, xlab=factorlist.in[i], ylab=y.lab, ylim=range(c(0, cis, par.res)), xaxt="n", cex.lab=1.3, cex.axis=1.3, col='grey', cex=0.5)
            points(1:length(newX), partialfit, pch=20)
          }else{
            if(type=='response'){
              y.lab = "Partial Fit (response)"
            }else{y.lab="Partial Fit (link)"}
          plot(1:length(newX), partialfit, pch=20,
               xlab=factorlist.in[i], ylab=y.lab, lwd=2, xaxt="n", ylim=range(c(0,cis)), cex.lab=1.3, cex.axis=1.3)
          }
          axis(1, at=(1:length(newX)), labels=newX)
          segments(1:length(newX),cis[,1], 1:length(newX) , cis[,2], lwd=2)
          abline(h=0, lwd=2)
        }

        if(savedata==T){
          partialdata<-data.frame(newX, partialfit, cis[,1], cis[,2])
          save(partialdata, file=paste('PartialData_', factorlist.in[i], '.RData', sep=''), compress='bzip2' )
        }

      }
    } # factor
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Continuous covariates
    if(is.null(varlist.in)==F){
      n<-length(varlist.in)
      for(i in 1:n){
        if(includeB0){
          coefpos<- c(1,grep(varlist.in[i], colnames(model.matrix(model))))  
        }else{
          coefpos<- c(grep(varlist.in[i], colnames(model.matrix(model))))
        }
        xvals<-data[,which(names(data)==varlist.in[i])]
        newX<- seq(min(xvals),max(xvals), length=500)
        eval(parse(text=paste(varlist.in[i], "<- newX", sep="")))

        #if(dim(initialModel$model[[1]])[2]==2){
 #         response<-
#        }else{
          response<- rep(1, 500)
  #      }

        newBasis<- eval(parse(text=labels(terms(model))[grep(varlist.in[i], labels(terms(model)))]))
        
        if(is.null(dim(newBasis))){
          newBasis<-as.matrix(newBasis)
        }
        
        if(includeB0){
          partialfit<- cbind(rep(1,500),newBasis)%*%coef(model)[coefpos]
        }else{
          partialfit<- newBasis%*%coef(model)[coefpos]
        }
        rcoefs<- NULL
        try(rcoefs<- rmvnorm(1000,coef(model), summary(model)$cov.robust), silent=T)
        if(is.null(rcoefs) || length(which(is.na(rcoefs)==T))>0){
          rcoefs<- rmvnorm(1000,coef(model), as.matrix(nearPD(summary(model)$cov.scaled)$mat))}
        if(includeB0){
          rpreds<- cbind(rep(1,500),newBasis)%*%t(rcoefs[,coefpos])
        }
        else{
          rpreds<- newBasis%*%t(rcoefs[,coefpos])
        }
        
        
        quant.func<- function(x){quantile(x, probs=c(0.025, 0.975))}
        cis<- t(apply(rpreds, 1,quant.func))
        if(type=='response'){
          partialfit <- model$family$linkinv(partialfit)
          cis <- model$family$linkinv(cis)
        }
        if(partial.resid==T){
          par.res <- residuals(model, type='partial')
          par.res <- par.res[, c(grep(varlist.in[i], colnames(par.res)))]
          if(type=='response'){
            par.res<-model$family$linkinv(par.res)
            y.lab = "Partial Fit (response)"
          }else{y.lab="Partial Fit (link)"}
          plot(xvals, par.res, xlab=varlist.in[i], ylab=y.lab, ylim=range(c(cis, par.res)), cex.lab=1.3, cex.axis=1.3, col='grey', cex=0.5)
          lines(newX, partialfit, lwd=2)
        }else{
          if(type=='response'){
            y.lab = "Partial Fit (response)"
          }else{y.lab="Partial Fit (link)"}
          plot(newX, partialfit, type="l",
               xlab=varlist.in[i], ylab=y.lab, lwd=2, ylim=range(cis), cex.lab=1.3, cex.axis=1.3)
        }
        if(model$family$family=='binomial'){
          if(length(unique(data$response)) ==2 | length(unique(data$successes))==2){
            rug(xvals[data$response==0], side=1)
            rug(xvals[data$response==1], side=3)
          }else{
            rug(xvals)
          }
        }else{
          rug(xvals)
        }
        lines(newX,cis[,1], col="darkred", lty=4)
        lines(newX,cis[,2], col="darkred", lty=4)

        if(showKnots=="TRUE"){
          spid<-grep(varlist.in[i], sapply(splineParams, '[[','covar'))
          if(length(spid)==0){
            XX<-FALSE
          }else{
            if(length(splineParams[[spid]]$knots)==1 & is.character(splineParams[[spid]]$knots)){
              XX<-FALSE
            }else{XX<-TRUE}  
          }
          

          if(XX==TRUE){
           
            abline(v=splineParams[[spid]]$knots, lty=4, col="grey")}
        }
        eval(parse(text=paste("rm(", varlist.in[i], ")", sep="")))
        rm(response)

        if(savedata==T){
          partialdata<-data.frame(newX, partialfit,  cis[,1], cis[,2])
          save(partialdata, file=paste('PartialData_', varlist.in[i], '.RData', sep=''), compress='bzip2' )
        }
      }
    } # varlist.in

  if(save==T){dev.off()}else{devAskNewPage(ask=FALSE)}


}



runPartialPlots_mn<-function(model, data, factorlist.in=NULL, varlist.in=NULL, showKnots=FALSE, type='link', partial.resid=FALSE, save=FALSE, savedata=F, label=NULL, includeB0=FALSE){
  
  print("Making partial plots")
  
  require(mvtnorm)
  require(splines)
  require(Matrix)
  require(abind)
  
  if (isS4(model)) {
    splineParams <- model@splineParams
  } else {
    splineParams <- model$splineParams
  }
  
  if(save==T){
    if(type=='response'){
      png(paste(label,"PartialFitsResponse%i.png",sep=''), height=600, width=650)}else{
        png(paste(label,"PartialFitsLink%i.png",sep=''), height=600, width=650)
      }}else{devAskNewPage(ask=TRUE)}
  
  par(mfrow=c(1,1))
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(factorlist.in)==F){
    for(i in 1:length(factorlist.in)){
      coeffac<- c(grep(factorlist.in[i], colnames(model.matrix(model))))
      coefradial<-c(grep('LRF', colnames(model.matrix(model))))
      coefpos <- coeffac[which(is.na(match(coeffac, coefradial)))]
      xvals<-data[,which(names(data)==factorlist.in[i])]
      newX<- sort(unique(xvals))
      newX<- newX[2:length(newX)]
      if(class(newX)=='character'){
        newX<-as.factor(newX)  
      }
      
      partialfit<- coef(model)[c(coefpos)]
      rcoefs<- NULL
      try(rcoefs<- rmvnorm(1000,coef(model), summary(model)$cov.robust), silent=T)
      if(is.null(rcoefs) || length(which(is.na(rcoefs)==T))>0){
        rcoefs<- rmvnorm(1000,coef(model), as.matrix(nearPD(summary(model)$cov.scaled)$mat))}
      
      rpreds<- rcoefs[,c(coefpos)]
      quant.func<- function(x){quantile(x, probs=c(0.025, 0.975))}
      # factor variables have rpreds vector
      if(is.null(dim(rpreds))){
        cis<-quant.func(rpreds)
        if(type=='response'){
          partialfit <- model@family@linkinv(partialfit)
          cis <- model@family@linkinv(cis)
        }
        if(partial.resid==T){
          par.res <- residuals(model, type='partial')
          par.res <- par.res[, c(grep(factorlist.in[i], colnames(par.res)))]
          if(type=='response'){
            par.res<-model@family@linkinv(par.res)
            y.lab = "Partial Fit (response)"
          }else{y.lab="Partial Fit (link)"}
          plot(xvals, par.res, xlab=factorlist.in[i], ylab=y.lab, ylim=range(c(0, cis, par.res)), xaxt="n", cex.lab=1.3, cex.axis=1.3, col='grey', cex=0.5)
          points(newX, partialfit, pch=20)
        }else{
          if(type=='response'){
            y.lab = "Partial Fit (response)"
          }else{y.lab="Partial Fit (link)"}
          plot(newX, partialfit,xlab=factorlist.in[i], ylab=y.lab,  xaxt="n", ylim=c(range(0,cis)), cex.lab=1.3, cex.axis=1.3)
          points(1:length(newX), partialfit, pch=20)
        }
        axis(1, at=newX, labels=colnames(model.matrix(model))[c(coefpos)])
        segments(1:length(newX),cis[1], 1:length(newX) , cis[2], lwd=2)
        abline(h=0, lwd=2)
      }else{
        # cts variables have matrix of rpreds
        cis<- t(apply(rpreds, 2,quant.func))
        if(type=='response'){
          partialfit <- model@family@linkinv(partialfit)
          cis <- model@family@linkinv(cis)
        }
        if(partial.resid==T){
          par.res <- residuals(model, type='partial')
          par.res<- par.res[, c(grep(factorlist.in[i], colnames(par.res)))]
          if(type=='response'){
            par.res<-model@family@linkinv(par.res)
            y.lab = "Partial Fit (response)"
          }else{y.lab="Partial Fit (link)"}
          plot(xvals, par.res, xlab=factorlist.in[i], ylab=y.lab, ylim=range(c(0, cis, par.res)), xaxt="n", cex.lab=1.3, cex.axis=1.3, col='grey', cex=0.5)
          points(1:length(newX), partialfit, pch=20)
        }else{
          if(type=='response'){
            y.lab = "Partial Fit (response)"
          }else{y.lab="Partial Fit (link)"}
          plot(1:length(newX), partialfit, pch=20,
               xlab=factorlist.in[i], ylab=y.lab, lwd=2, xaxt="n", ylim=range(c(0,cis)), cex.lab=1.3, cex.axis=1.3)
        }
        axis(1, at=(1:length(newX)), labels=newX)
        segments(1:length(newX),cis[,1], 1:length(newX) , cis[,2], lwd=2)
        abline(h=0, lwd=2)
      }
      
      if(savedata==T){
        partialdata<-data.frame(newX, partialfit, cis[,1], cis[,2])
        save(partialdata, file=paste('PartialData_', factorlist.in[i], '.RData', sep=''), compress='bzip2' )
      }
      
    }
  } # factor
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Continuous covariates
  if(is.null(varlist.in)==F){
    n<-length(varlist.in)
    
    rcoefs<- NULL
    if (isS4(model)) {
      setClass("vglmMRSea", contains=c("vglm"), slots=c(varshortnames="character", panels="ANY", splineParams="list", data="data.frame", cvfolds="numeric", interactionterm="ANY")) -> vglmMRSea
      sum_mod <- summaryvglm(model)
      covs <- sum_mod@cov.unscaled
      try(rcoefs <- rmvnorm(1000,coef(model), covs), silent=T)
      if(is.null(rcoefs) || length(which(is.na(rcoefs)==T))>0){
        rcoefs<- rmvnorm(1000,coef(model), as.matrix(nearPD(covs)$mat))}
    } else {
      try(rcoefs<- rmvnorm(1000,coef(model), summary(model)$cov.robust), silent=T)
      if(is.null(rcoefs) || length(which(is.na(rcoefs)==T))>0){
        rcoefs<- rmvnorm(1000,coef(model), as.matrix(nearPD(summary(model)$cov.scaled)$mat))}
    }
    
    for(i in 1:n){
      if(includeB0){
        coefpos<- c(grep("Intercept", colnames(model.matrix(model))),grep(varlist.in[i], colnames(model.matrix(model))))  
      }else{
        coefpos<- c(grep(varlist.in[i], colnames(model.matrix(model))))
      }
      
      xvals<-data[,which(names(data)==varlist.in[i])]
      newX<- seq(min(xvals),max(xvals), length=500)
      eval(parse(text=paste(varlist.in[i], "<- newX", sep="")))
      
      #if(dim(initialModel$model[[1]])[2]==2){
      #         response<-
      #        }else{
      response<- rep(1, 500)
      #      }
      
      newBasis<- eval(parse(text=labels(terms(model))[grep(varlist.in[i], labels(terms(model)))]))
      
      if(is.null(dim(newBasis))){
        newBasis<-as.matrix(newBasis)
      }
      
      # reshape or not coefficients for multinomial models
      if (isS4(model)){
        n_y_levels <- dim(model@y)[2] - 1
        coef_shaped <- coef(model)[coefpos]
        coef_shaped <- matrix(coef_shaped, ncol=n_y_levels, byrow=TRUE)
      } else {
        coef_shaped <- coef(model)[coefpos]
      }
      
      if(includeB0){
        partialfit<- cbind(rep(1,500),newBasis)%*%coef_shaped
      }else{
        partialfit<- newBasis%*%coef_shaped
      }
      
      if (isS4(model)) {
        rpreds <- array(rep(NA, n_y_levels*500*1000), c(500, 1000, n_y_levels))
        for (ct in 1:n_y_levels){
          tot_coefs <- length(coefpos)
          coefpos_lev <- coefpos[seq(ct, tot_coefs, by=n_y_levels)]
          if(includeB0){
            rpreds[,,ct] <- cbind(rep(1,500),newBasis)%*%t(rcoefs[,coefpos_lev])
          } else{
            rpreds[,,ct] <- newBasis%*%t(rcoefs[,coefpos_lev])
          }
        }
      } else {
        if(includeB0){
          rpreds<- cbind(rep(1,500),newBasis)%*%t(rcoefs[,coefpos])
        } else{
          rpreds<- newBasis%*%t(rcoefs[,coefpos])
        }
      }
      
      
      quant.func<- function(x){quantile(x, probs=c(0.025, 0.975))}
      if (isS4(model)){
        cis <- array(rep(NA, n_y_levels*500*2), c(500, 2, n_y_levels))
        for (ct in 1:n_y_levels){
          cis[,,ct] <- t(apply(rpreds[,,ct], 1,quant.func))
        }
      } else {
        cis<- t(apply(rpreds, 1,quant.func))
      }
      if(type=='response'){
        if (isS4(model)){
          partialfit <- model@family@linkinv(partialfit, extra=model@extra)
          ci_low <- cis[,1,]
          ci_low <- model@family@linkinv(ci_low, extra=model@extra)
          ci_hi <- cis[,2,]
          ci_hi <- model@family@linkinv(ci_hi, extra=model@extra)
          cis <- abind(ci_low, ci_hi, along=3)
        } else {
          partialfit <- model$family$linkinv(partialfit)
          cis <- model$family$linkinv(cis)
        }
      } else {
        if (isS4(model)){
          ci_low <- cis[,1,]
          ci_hi <- cis[,2,]
          cis <- abind(ci_low, ci_hi, along=3)
        } else {
          cis <- cis
        }
      }
      # need to separate to create function to create graph
      # then need to put in data for each y response
      #if(partial.resid==T){
      # this doesn't work for multinomial
      #par.res <- residuals(model, type='partial')
      #par.res <- par.res[, c(grep(varlist.in[i], colnames(par.res)))]
      #if(type=='response'){
      #par.res<-model$family$linkinv(par.res)
      #y.lab = "Partial Fit (response)"
      #}else{y.lab="Partial Fit (link)"}
      #plot(xvals, par.res, xlab=varlist.in[i], ylab=y.lab, ylim=range(c(cis, par.res)), cex.lab=1.3, cex.axis=1.3, col='grey', cex=0.5)
      #lines(newX, partialfit, lwd=2)
      #}
      
      if (isS4(model)) {
        yvarz <- colnames(model@y)
        if (type=="response") {
          for (yy in 1:length(yvarz)) {
            create_graph_s4(partialfit[, yy], cis[,yy,], varlist.in[i], yvarz[yy], newX, xvals, type, splineParams, showKnots)
          }
        } else {
          nn_levs <- 
            for (yy in 1:n_y_levels) {
              create_graph_s4(partialfit[, yy], cis[,yy,], varlist.in[i], yvarz[yy], newX, xvals, type, splineParams, showKnots)
            }
        }
      } else {
        create_graph_s3(partialfit, cis, varlist.in[i], newX, xvals, type, splineParams, model$family$family, data, showKnots)
      }
      
      
      eval(parse(text=paste("rm(", varlist.in[i], ")", sep="")))
      rm(response)
      
      if(savedata==T){
        partialdata<-data.frame(newX, partialfit,  cis[,1], cis[,2])
        save(partialdata, file=paste('PartialData_', varlist.in[i], '.RData', sep=''), compress='bzip2' )
      }
    }
  } # varlist.in
  
  if(save==T){dev.off()}else{devAskNewPage(ask=FALSE)}
  
  
}


create_graph_s3 <- function(partfit, ciz, var.in, newX, xvals, type, splineParams, family, data, showKnots=FALSE){
  if(type=='response'){
    y.lab = "Partial Fit (response)"
  }else{
    y.lab="Partial Fit (link)"
  }
  plot(newX, partfit, type="l", xlab=var.in, ylab=y.lab, lwd=2, ylim=range(ciz), cex.lab=1.3, cex.axis=1.3)
  #if(family=='binomial'){
  #  if(length(unique(data$response)) ==2 | length(unique(data$successes))==2){
  #    rug(xvals[data$response==0], side=1)
  #    rug(xvals[data$response==1], side=3)
  #  }else{
  #    rug(xvals)
  #  }
  #}else{
  #  rug(xvals)
  #}
  lines(newX,ciz[,1], col="darkred", lty=4)
  lines(newX,ciz[,2], col="darkred", lty=4)
  if(showKnots==TRUE){
    spid<-grep(var.in, sapply(splineParams, '[[','covar'))
    if(length(spid)==0){
      XX<-FALSE
    }else{
      if(length(splineParams[[spid]]$knots)==1 & is.character(splineParams[[spid]]$knots)){
        XX<-FALSE
      }else{XX<-TRUE}  
    }
    if(XX==TRUE){
      abline(v=splineParams[[spid]]$knots, lty=4, col="grey")}
  }
}

create_graph_s4 <- function(partfit, ciz, var.in, yvar, newX, xvals, type, splineParams, showKnots=FALSE){
  if(type=='response'){
    y.lab = "Partial Fit (response)"
  }else{
    y.lab="Partial Fit (link)"
  }
  print(range(ciz))
  plot(newX, partfit, type="l", xlab=paste(var.in, yvar), ylab=y.lab, lwd=2, ylim=range(ciz), cex.lab=1.3, cex.axis=1.3)
  rug(xvals)
  lines(newX,ciz[,1], col="darkred", lty=4)
  lines(newX,ciz[,2], col="darkred", lty=4)
  if(showKnots==TRUE){
    spid<-grep(var.in, sapply(splineParams, '[[','covar'))
    if(length(spid)==0){
      XX<-FALSE
    }else{
      if(length(splineParams[[spid]]$knots)==1 & is.character(splineParams[[spid]]$knots)){
        XX<-FALSE
      }else{XX<-TRUE}  
    }
    if(XX==TRUE){
      abline(v=splineParams[[spid]]$knots, lty=4, col="grey")}
  }
}
