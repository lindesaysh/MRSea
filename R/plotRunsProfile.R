#-----------------------------------------------------------------------------
#' Calculate runs test and plot profile plot.  The output is a plot of runs profiles (with p-value to indicate level of correlation)
#' 
#' @param model Fitted model object (glm or gam)
#' @param varlist Vector of covariate names (continous covariates only)
#' @param label Label printed at the end of the plot name to identify it when \code{save=TRUE}.
#' @param save (\code{default=FALSE}). Logical stating whether plot should be saved into working directory.
#' 
#' @return
#' Runs profile plots are returned for residuals ordered by each covariate in \code{varlist}, predicted value and index of observations (temporally).
#' 
#' The black line is the line of sequences of positive or negative residuals.  The vertical lines are the change between a sequence of positive to negative residuals (or vice versa).
#' 
#' The p-values are from a \code{\link{runsTest}} and indicate whether there is correlation in the residuals (p<0.05) or independence (p>0.05).  The test statistic determines the type of correlation (positive/negative) and the result printed at the bottom of the figure. 
#' 
#' Note: if the covariate is discrete in nature (like the example below), there will be a lot of overplotting of runs.  Some jittering occurs at each discrete value (for covariates with <= 25 unique values).
#' 
#' 
#' @examples 
#' # load data
#' data(ns.data.re)
#' 
#' model<-gamMRSea(birds ~ observationhour + as.factor(floodebb) + as.factor(impact), 
#'             family='quasipoisson', data=ns.data.re)
#' 
#' plotRunsProfile(model, varlist=c('observationhour'))
#' 
#' @export
#' 
plotRunsProfile<- function(model, varlist, label='', save=FALSE){
  
  print("Calculating runs test and plotting profile")
  
  namesOfx<- c(varlist, c("Predicted", "Index"))
  
  if(class(model)[1]=='geeglm' | class(model)[1]=='glm' | class(model)[1]=='gamMRSea'){
    dat<- model$data
  }
  
  if(class(model)[1]=='gam'){
    dat<-model$model 
  }
  
  # find which data columns refer to the used variables
  coefpos<-c()
  for(c in 1:(length(namesOfx)-2)){
    coefpos<- c(coefpos, grep(namesOfx[c], names(dat)))
  }
  
  # make matrix (n x number covar + fitted + orderdata)
  xmatrix<- cbind(dat[,coefpos], fitted(model), ord=1:length(fitted(model)))
  
  for(c in 1:length(namesOfx)){
    rtest<-runsTest(residuals(model, type="pearson")[order(xmatrix[, c])],alternative = c("two.sided"))
    pval<-round(rtest$p.value, 3)
    if(rtest$p.value<0.0001){pval<-'< 0.0001'}
    
    
    b<-xmatrix[,c][order(xmatrix[, c])]
    # if the covariate is very discrete, the x's need to be jittered in order
    if(length(unique(xmatrix[,c]))<25){
      a<-xmatrix[,c][order(xmatrix[, c])]
      tab<-table(a)
      seqrange<- (unique(a)[2] - unique(a)[1])/10
      # jitter discrete variable
      b<-a
      for(j in 1:length(unique(a))){
        b[a==unique(a)[j]]<- unique(a)[j] + seq(-seqrange, seqrange, length=tab[j])
      }
    }
    
    if(save==T){png(paste("RunsProfile_", namesOfx[c], label, ".png", sep=''), height=600, width=700)}else{devAskNewPage(ask=TRUE)}
    plot(stepfun(b, c(0,sign(residuals(model, type="pearson")[order(xmatrix[, c])]))), main=paste("Runs Profile, p=",pval, sep=""), xlab=namesOfx[c], pch=20, ylab="Sign of Pearsons residuals", cex.lab=1.3, cex.axis=1.3)
    #}
    abline(h=0)
    if(rtest$p.value<0.05){
      if(rtest$statistic<0){title(sub=paste("Significant positive correlation; p-value= ",pval, sep="" ), cex.sub=1.3)}
      if(rtest$statistic>0){title(sub=paste("Significant negative correlation; p-value= ",pval, sep="" ), cex.sub=1.3)}
    }
    if(save==T){dev.off()}else{devAskNewPage(ask=FALSE)}
  }
}