
#-----------------------------------------------------------------------------
#' run functions to create acf matrix and plot the results
#' 
#' @param block Vector of blocks that identify data points that are correlated
#' @param model Fitted model object (glm or gam)
#' @param store (\code{default=F}). Logical stating whether a list of the matrix of correlations is stored (output from \code{acffunc}.)
#' @param save (\code{default=FALSE}). Logical stating whether plot should be saved into working directory.
#' 
#' @return
#' Plot of lag vs correlation.  Each grey line is the correlation for each individual block in \code{block}.  The red line is the mean values for each lag.
#' 
#' If \code{store=TRUE} then the matrix of correlations (nblocks x length_max_block) is returned and \code{plotacf} may be used to plot the acf.
#' 
#' 
#' @examples 
#' # load data
#' data(ns.data.re)
#' 
#' model<-glm(birds ~ observationhour + as.factor(floodebb) + as.factor(impact), 
#'            family='quasipoisson', data=ns.data.re)
#' 
#' ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear, 
#'                     ns.data.re$DayOfMonth, sep='')
#' ns.data.re$blockid<-as.factor(ns.data.re$blockid)
#' 
#' runACF(ns.data.re$blockid, model)
#' 
#' @export
#' 

runACF<-function(block, model, store=FALSE, save=F, suppress.printout=FALSE){
  acf_result<-acffunc(block, model, suppress.printout)
  if(save==T){png('acfPlot.png', height=500, width=600)}
  plotacf(acf_result$acfmat)
  if(save==T){dev.off()}
  if(store==TRUE){return(acf_result)}
}


#-----------------------------------------------------------------------------
#' calculate correlation for residuals by block
#' 
#' @param block Vector of blocks that identify data points that are correlated
#' @param model Fitted model object (glm or gam)
#' 
acffunc<-function(block, model, suppress.printout=FALSE){
  blocktab<-table(block)
  acfmat<-matrix(NA, length(unique(block)), max(blocktab))
  for(i in 1:length(unique(block))){
    
    if(suppress.printout==FALSE){
      print(i)
    }
    
    corr<-as.vector(acf(residuals(model, type='pearson')[which(block==unique(block)[i])], plot=F,lag.max=max(blocktab))$acf)
    acfmat[i,1:length(corr)]<- corr
  }
  return(list(acfmat=acfmat, blocktab=blocktab))
}


#-----------------------------------------------------------------------------
# plot correlation of residuals by block
#-----------------------------------------------------------------------------
#' run functions to create acf matrix and plot the results
#' @param acfmat Matrix of output from \code{acffunc} (blocks x max block length).
#'  
plotacf<-function(acfmat){
  plot(1:length(na.omit(acfmat[1,])), na.omit(acfmat[1,]), xlim=c(0,ncol(acfmat)), ylim=c(-1,1), type='l', col='grey', xlab='Lag', ylab='Auto correlation', cex.lab=1.3, cex.axis=1.3)
  abline(h=0)
  for(i in 2:nrow(acfmat)){
    lines(1:length(na.omit(acfmat[i,])), na.omit(acfmat[i,]), col='grey')  
  }
  lines(1:ncol(acfmat), apply(acfmat, 2, mean, na.rm=T), col='red', lwd=2)
}
