#' Determining the distance bin 
#'
#' For a vector of perpendicular (or radial) distances, this function determines which distance bin it belongs to 
#' (given the input of cut points) and adds the beginning and end points of the respective distance bins in new colunns
#' in \code{dis.data} called "distbegin" and "distend". 
#'
#' @param dis.data A data frame with distance data for which perpendicular (or radial) distances are recorded in the \code{distance} column
#' @param cutpoints A vector of cut points of the intervals (this function is not set up to deal with left-truncation)
#'
#' @details
#' If a value in \code{dis.data$distance} matches a cut point in \code{cutpoints} exactly,  the value of \code{dis.data.re$distance} will be attributed to the bin that is closer to the line/point unless the value of \code{dis.data.re$distance} is 0. 
#' 
#' E.g. if \code{cutpoints=c(0,1,2,3)}, \code{dis.data$distance}=2 will be attributed to interval 2 (and not 3). 
#'
#' @return
#' The \code{dis.data} data frame to which columns "distbegin" and "distend" were added giving the beginning and end cutpoints 
#' of the bin that the respective \code{dis.data$distance} belongs to. 
#' 
#' @export
#' 
which.bin<-function(dis.data,cutpoints){
bins<-matrix(0,length(cutpoints)-1,2)
for (i in 1:(length(cutpoints)-1)){
bins[i,1:2]<-cutpoints[i:(i+1)]
}
y<-array(NA,length(dis.data$distance))
yz<-matrix(NA,length(dis.data$distance),2)
for (i in which(is.na(dis.data$distance)==F)){
      if (dis.data$distance[i]==0) {y[i]<-1}
      else {y[i]<-which(bins[,1]< dis.data$distance[i] & bins[,2]>=dis.data$distance[i])}
      yz[i,]<-bins[y[i],]
}
dis.data[,c("distbegin","distend")]<-yz
return(dis.data)
}

