#-----------------------------------------------------------------------------
#' functions to create observed vs fitted and fitted vs scaled pearsons residual plots
#' 
#' @param model Fitted model object (glm or gam)
#' @param plotting Plotting options (\code{default='b'}). \code{b}: returns both plots, \code{f}: returns observed vs fitted only and  \code{r}: returns scale pearsons residual plot only.
#' @param save (\code{default=FALSE}). Logical stating whether plot should be saved into working directory. See \code{label} to change directory.
#' @param label Character string indicating an label to be added to the plot when using \code{save = TRUE}. Can also include a pathway to a directory of choice.
#' 
#' @return
#' Two plots:
#' \item{Observed vs Fitted}{Plot of observed vs fitted with concordence correlation and marginal R-squared printed in the plot title.}
#' \item{Fitted vs scaled Pearsons residuals}{The red line is a locally weighted least squares regression line of all of the residuals.}
#' 
#' @examples 
#' # load data
#' data(ns.data.re)
#' 
#' model<-glm(birds ~ observationhour + as.factor(floodebb) + as.factor(impact), 
#'            family='quasipoisson', data=ns.data.re)
#' 
#' runDiagnostics(model)
#'
#' @export
#'

runDiagnostics<-function(model, plotting='b', save=FALSE, label = NULL){
  print("Assessing predictive power")
  
  response<-model$y
  
  if(class(model)[1]=='geeglm'){
    dat<- model$data
  }
  
  if(class(model)[1]=='gam'){
    dat<-model$model 
  }
  
  if(class(model)[1]=='gamMRSea'){
    dat<-model$data
  }
  
  #Assessing predictive power
  
  #r-squared:
  r2<-1-(sum((response-fitted(model))**2)/sum((response-mean(response))**2))
  #concordence correlation
  num<- 2*sum((response-mean(response))*(fitted(model)-mean(fitted(model))))
  den<- sum((response-mean(response))**2)+sum((fitted(model)-mean(fitted(model)))**2)
  rc<-num/den
  
  scaledRes<- residuals(model, type="response")/
    sqrt(family(model)$variance(fitted(model))*as.numeric(summary(model)$dispersion[1]))
  
  df<-data.frame(fits=fitted(model), res=scaledRes, smx=lowess(fitted(model), scaledRes)$x, smy=lowess(fitted(model), scaledRes)$y, response=response )
  
  if(plotting=='b' | plotting=='f'){
    a<- ggplot(df) + geom_point(aes(response, fits), alpha=0.15, size=3) + geom_abline(intercept=0, slope=1) + labs(x='Observed Values', y='Fitted Values', title=paste("Concordence correlation: ", round(rc,4), "\nMarginal R-squared value: ", round(r2,4), sep=""))
    
    a<-a + theme_bw() + theme(panel.grid.major=element_blank(), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), plot.title=element_text(size=15))
    
    if(save==T){ggsave(paste0(label, "FitPlots_fitted.png"), a, height=6, width=8)}else{
      if(plotting=='b'){devAskNewPage(ask=TRUE)}
      plot(a)}
  }
  
  if(plotting=='b' | plotting=='r'){
        
    b<- ggplot(df) + geom_point(aes(fits, res), alpha=0.15, size=3)+ geom_line(aes(smx, smy), col='red') + geom_abline(intercept=0, slope=0) + labs(x='Fitted Values', y='Scaled Pearsons Residuals')
    
    b<- b + theme_bw() + theme(panel.grid.major=element_blank(), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15))
    
    if(save==T){ggsave(paste0(label, "FitPlots_resids.png"), b, height=6, width=8)}else{plot(b)}
    
  }
  
  devAskNewPage(ask=FALSE)
}




# cumResFun_diag<- function(model, xvals, varlist, i, type='response'){
#   maxy<- max(residuals(model),cumsum(residuals(model, type=type)[order(xvals)]))
#   miny<- min(residuals(model),cumsum(residuals(model, type=type)[order(xvals)]))
#   plot(xvals[order(xvals)],residuals(model, type=type)[order(xvals)], ylim=c(miny,maxy), xlab=varlist[i], ylab="Response residuals", pch=20, col="darkgrey", main="Cumulative Residuals", cex.lab=1.3, cex.axis=1.3)
#   lines(xvals[order(xvals)],cumsum(residuals(model, type=type)[order(xvals)]), lwd=2)
#   abline(h=0)
# }
# 
# 
# DiagPlots<- function(model, dat, varlist, i){
#   xvals<-dat[,which(names(dat)==varlist[i])] 
#   cumResFun_diag(model, xvals, varlist, i)
#   rtest<-runs.test(residuals(model, type="pearson")[order(xvals)],alternative = c("two.sided"))
#   if(nrow(dat)<10000){
#     plot(stepfun(xvals[order(xvals)], c(0,sign(residuals(model, type="pearson")[order(xvals)]))), main=paste("Runs Profile, p=", round(rtest$p.value,3), sep=""), xlab=varlist[i], pch=20, ylab="Sign of Pearsons residuals", cex.lab=1.3, cex.axis=1.3)
#   }
#   abline(h=0)
#   if(rtest$p.value<0.05){
#     if(rtest$statistic<0){title(sub=paste("Significant positive correlation; p-value= ",rtest$p.value, sep="" ), cex.sub=1.3)}
#     if(rtest$statistic>0){title(sub=paste("Significant negative correlation; p-value= ",rtest$p.value, sep="" ), cex.sub=1.3)}
#   }
# }
# 
# 
# 
# 
# 
# 
