#' Estimated number of individuals for each detection
#'
#' This function creates a new column in \code{data} which contains the estimated 
#' number of animals for each detection. This is the number of observed individuals 
#' divided by their probability of detection using MCDS methods (size/detection probability). 
#' In the case that no \code{size} column is given in \code{dis.data}, it is assumed that 
#' detections were made of individuals and \code{size} is set to 1 for all detections. The values 
#' for \code{size} and \code{NHAT} are set to zero in case the distance was larger than the 
#' truncation distance \code{w} specified in \code{det.fct.object}.
#' In addition, a new column \code{area} is created which is used as the offset in the 
#' second stage count model (segment length * (truncation distance/1000) * 2). The truncation
#' distance is divided by 1000 to convert it from metres to km. It is assumed that the 
#' segment data represents two-sided surveys. In case the survey was one-sided, this column needs to 
#' be divided by 2 after the call to this function. 
#'
#' @param data distance data object used with \code{det.fct} to estimate probabilities of detection
#' @param ddf.obj detection function object created by \code{ddf}
#'
#' @examples 
#' data(dis.data.re)
#' result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re,method="ds", 
#'      meta.data= list(width=250,binned=FALSE))
#' dis.data<-create.NHAT(dis.data.re,result)
#' 
#' @export

create.NHAT <-
function(data, ddf.obj) {
        within.w<-which(data$distance<=ddf.obj$meta.data$width)     # this ensures that the length of the detections included in the analysis matches the length of ddf.obj$fitted
        data$area<-NA
        data$area<-data$length*2*ddf.obj$meta.data$width/1000       # the size of the covered area (used as the offset in the count model)
        if(length(data$size)==0){                                   # in case the column 'size' is not specified in data it is assumed that the size of each cluster equals one
        data$size<-0
        data$size[within.w]<-1}                
        data$NHAT<-0
        data$NHAT[within.w]<-data$size[within.w]/ddf.obj$fitted     # dividing the size of each cluster by its probability of detection
        data
        }
