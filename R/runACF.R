
#-----------------------------------------------------------------------------
#' run functions to create acf matrix and plot the results
#' 
#' @param block Vector of blocks that identify data points that are correlated
#' @param model Fitted model object (glm or gam)
#' @param store (\code{default=FALSE}). Logical stating whether a list of the matrix of correlations is stored (output from \code{acffunc}.)
#' @param save (\code{default=FALSE}). Logical stating whether plot should be saved into working directory.
#' @param suppress.printout (\code{default=TRUE}. Logical stating whether to show a printout of block numbers to assess progress. `FALSE` will show printout.
#' @param  maxlag (\code{default=NULL}). Numeric entry to allow the restriction of the maximum lag on the plots.  If \code{NULL} then the length of the longest panel is used as the maximum plotted lag. 
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
#' model<-gamMRSea(birds ~ observationhour + as.factor(floodebb) + as.factor(impact), 
#'            family='quasipoisson', data=ns.data.re)
#' 
#' ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear, 
#'                     ns.data.re$DayOfMonth, sep='')
#' ns.data.re$blockid<-as.factor(ns.data.re$blockid)
#' 
#' runACF(ns.data.re$blockid, model, suppress.printout=TRUE)
#' 
#' # storing the output and then plotting
#' acfoutput <- runACF(ns.data.re$blockid, model, suppress.printout=TRUE, store=TRUE)
#' plotacf(acfoutput$acfmat)
#' 
#' @author LAS Scott-Hayward, University of St Andrews
#' 
#' @export
#' 

runACF<-function(block, model, store=FALSE, save=F, suppress.printout=TRUE, maxlag=NULL, printplot=TRUE){
  acf_result<-acffunc(block, model, suppress.printout)
  if(save==T){
    png('acfPlot.png', height=500, width=600)
    plotacf(acf_result$acfmat, maxlag)
    dev.off()
  }else{
    if(store==TRUE){
      return(acf_result)
    }else{
      plotacf(acf_result$acfmat, maxlag)
    }
  }
}


#-----------------------------------------------------------------------------
#' calculate correlation for residuals by block
#' 
#' @param block Vector of blocks that identify data points that are correlated
#' @param model Fitted model object (glm or gam)
#' @param suppress.printout (Default: \code{TRUE}. Logical stating whether to show a printout of block numbers to assess progress. `FALSE` will show printout.
#' 
acffunc<-function(block, model, suppress.printout=TRUE){
  blocktab<-table(block)
  acfmat<-matrix(NA, length(unique(block)), max(blocktab))
  
  if(is.list(model)){
    d<-residuals(model, type='pearson')
  }else{
    d<-model
  }
  
  overallacf<-acf(d, lag.max = max(blocktab), plot=F)$acf
  
  
    for(i in 1:length(unique(block))){
    
    if(suppress.printout==FALSE){
      print(i)
    }
    
    corr<-as.vector(acf(d[which(block==unique(block)[i])], plot=F,lag.max=max(blocktab))$acf)
    if(length(which(is.na(corr)))>0)
    {
      corr<-overallacf[1:length(corr)]
    }
    acfmat[i,1:length(corr)]<- corr
  }
  return(list(acfmat=acfmat, blocktab=blocktab))
}


#-----------------------------------------------------------------------------
# plot correlation of residuals by block
#-----------------------------------------------------------------------------
#' run functions to create acf matrix and plot the results
#' @param acfmat Matrix of output from \code{acffunc} (blocks x max block length).
#' @param  maxlag (\code{default=NULL}). Numeric entry to allow the restriction of the maximum lag on the plots.  If \code{NULL} then the length of the longest panel is used as the maximum plotted lag.
#'  
plotacf<-function(acfmat, maxlag=NULL){
  
  suppressWarnings(library(dplyr, quietly = TRUE))
  
  if(is.null(maxlag)){
    maxlag = ncol(acfmat)
  }
  acfdat <- as_tibble(data.frame(t(acfmat))) %>% 
    mutate(Meancor = apply(acfmat, 2, mean, na.rm=T),
           Lag = row_number()-1) %>% 
    tidyr::pivot_longer(names_to = "blocks", values_to = "correlation", cols = -c(Lag)) %>%
    arrange(blocks, Lag)
  
  t <- round(filter(acfdat, Lag==1, blocks!="Meancor") %>% 
    summarise(min=min(correlation), max=max(correlation)),2)
  tmean = round(filter(acfdat, Lag==1, blocks=="Meancor")$correlation,2)
  
  suppressWarnings({
  p <- ggplot() +
    geom_line(data = filter(acfdat, blocks != "Meancor", Lag <= maxlag), 
              aes(x = Lag, y = correlation, group=blocks), 
              colour = "grey", linewidth = 1) +
    theme_bw() + 
    ylab("Auto correlation") +
    geom_hline(aes(yintercept=0)) + 
    geom_line(data = filter(acfdat, blocks == "Meancor", Lag <= maxlag), 
              aes(x = Lag, y = correlation, group=blocks), 
              colour = "red3", linewidth = 1) +
    ggtitle(paste0("Lag 1: min = ", t[1], ", mean = ", tmean, ", max = ", t[2]))
  print(p)
  })
}
