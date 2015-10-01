#-----------------------------------------------------------------------------
#' Identify any significant differences between predicted data before an impact event and predicted data after an impact event
#' 
#' 
#' @param beforePreds Matrix of bootstrap predictions (n x B) to each grid cell before impact (same length and order as \code{afterPreds})
#' @param afterPreds Matrix of bootstrap predictions (n x B) to each grid cell after impact (same length and order as \code{beforePreds})
#' @param quants (\code{default = =c(.025,.975)}) Quantile for significance.
#' 
#' @details
#' This function finds the differences for every predicted grid cell for every bootstrap replicate.  Quantiles are used to determine whether each difference is significantly different from zero and if so, in what direction. 
#' 
#' @return
#' A list is returned consisting of
#' \item{mediandiff}{Vector of the median difference for each grid cell}
#' \item{lowerci}{Vector of the lower 2.5\% difference for each grid cell}
#' \item{upperci}{Vector of the upper 97.5\% difference for each grid cell}
#' \item{significanceMarker}{Vector of significance.  0: not significant, 1: significant and positive, -1: significant and negative}
#' 
#' @examples 
#'  \dontrun{
#' getDifferences(beforePreds, afterPreds)}
#' 
#' @export
#' 
getDifferences<-function(beforePreds, afterPreds, quants=c(.025,.975)){
  
  # check for columns with infinity and remove
  rmid<-unique(which(beforePreds =='Inf', arr.ind=T)[,2])
  if(length(rmid)>0){
    beforePreds<-beforePreds[,-rmid]
    afterPreds<-afterPreds[,-rmid]
    cat('NOTE:', length(rmid), 'bootstrap(s) removed due to infinite values')
  }
  
  
  diff<- matrix(NA, nrow=nrow(beforePreds), ncol=ncol(beforePreds))
  for(i in 1:ncol(beforePreds)){
    diff[,i]<- afterPreds[,i] - beforePreds[,i]
  }
  # find lower and upper quantiles for difference
  difcis<-t(apply(diff, 1,  quantile, probs= quants, na.rm=T ))
  
  # find out if the confidence interval of the difference contains zero
  isin<-as.vector((apply(difcis, 1,  contains)))
  
  # make vector to define what is positive and negative significant differences
  marker<- rep(1, length=length(isin))
  marker[which(isin==0)] <- 0
  marker[which(isin>0 & difcis[,1]>0)] <- 1
  marker[which(isin>0 & difcis[,1]<0)] <- (-1)
  
  mediandiff<-apply(diff,1,median)
  
  return(list(mediandiff=mediandiff, lowerci=difcis[,1], upperci=difcis[,2],significanceMarker=marker))
}


contains<-function(x){
  if(x[1]<=0 & x[2]>=0){z<-0}else{z<-1}
  return(z)
}


