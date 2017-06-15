#' Generate a grid of knot locations to run SALSA2D.
#' 
#' @param coordData nx2 matrix or data frame of coordinates representing data locations
#' @param numKnots (\code{default = 300)}).  The user may choose how many legal knot locations are available (only if more than 400 )
#' @param plot (\code{default = TRUE}). Logical stating whether a plot showing the legal knot positions is given.
#' 
#' @details
#' SALSA2D requires a grid of knot locations to determine the best locations.  Illegal knot positions (those not close to data) are kept as a row in the data frame of locations but given c(NA, NA) to avoid a knot considered.
#' 
#' @return
#' A (numKnots x 2) matrix of knot locations.
#' 
#' @examples 
#'  \dontrun{
#' data(dis.data.re)
#' # bootstrap data without stratification
#' dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
#' require(mrds)
#' result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
#'              meta.data=list(width=250))
#' dis.data.re<-create.NHAT(dis.data.re,result)
#' 
#' knotgrid<-getKnotgrid(cbind(dis.data.re$x.pos, dis.data.re$y.pos))}
#' 
#' @export
#' 
getKnotgrid<-function(coordData, numKnots=300, plot=TRUE){
  
  # XGridPoints<- seq(min(coordData[,1]),
  #                   (max(coordData[,1])), length=100)
  # YGridPoints<- seq(min(coordData[,2]),
  #                   (max(coordData[,2])), length=100)
  # Grid<- data.frame(expand.grid(x=XGridPoints,y=YGridPoints))
  # 
  # #find the grid point closest to each candidate knot on the grid
  # d<- c()
  # for(i in 1:nrow(coordData)){
  #   d[i]<- which.min(sqrt((coordData[i,1]-Grid[,1])**2+(coordData[i,2]-Grid[,2])**2))
  # }
  # 
  # RowsToKeep<- rep(0,nrow(Grid))
  # RowsToKeep[d]<-1
  # 
  # GridPosIncludingNAs<- cbind(ifelse(RowsToKeep!=1, NA, RowsToKeep*Grid[,1]), ifelse(RowsToKeep!=1, NA, RowsToKeep*Grid[,2]))
  # Grid <- as.data.frame(GridPosIncludingNAs)
  # knotgrid<-Grid
  
  require(fields)
  
  dupPoints <-paste(coordData[,1], coordData[,2], sep='E')
  duppointid<-which(duplicated(dupPoints)==F)
  
  samplingdata<-coordData[duppointid,]
  
  if(nrow(samplingdata)>1000){
    sampleid<-sample(1:nrow(samplingdata), 1000)
    samplingdata<-samplingdata[sampleid,]
  }
  
  if(nrow(samplingdata)<numKnots){
    numKnots<-nrow(samplingdata)
    spaceid<-1:nrow(samplingdata)
    knotgrid <- samplingdata
  }else{
    spaceid <-cover.design(R = samplingdata, nd = numKnots, nruns = 5)$best.id
    knotgrid <- samplingdata[spaceid,]
  }
  
  if(!exists('sampleid')){
    attr(knotgrid, 'points.selected')<-duppointid[spaceid]
  }else{
    attr(knotgrid, 'points.selected')<-duppointid[sampleid[spaceid]]
  }
  
  # cut down number of knots if there are more than 400 locations available
  #if(nrow(na.omit(knotgrid))>400){
  #   naid<- which(is.na(knotgrid))
  #   spaceid<-cover.design(R = na.omit(knotgrid), nd = numKnots, nruns = 5)$best.id
  #   realid<-(1:nrow(knotgrid))[-naid]
  #   
  #   knotgrid[realid[-spaceid],]<- c(NA, NA)  
  # }
  
  if(plot==TRUE){
    plot(coordData[,1], coordData[,2], pch=20, cex=0.2)
    points(knotgrid, pch=20, cex=0.5, col='red')  
  }
  
  return(knotgrid)
}


