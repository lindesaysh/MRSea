#' Aggregate bootstrapped distance data into count data 
#'
#' This function creates a new data set where \code{dis.data} is aggregated 
#' for each visit to a segment. For bootstrapped data, the column with the ids 
#' for visits to a segment is \code{segment.id2} which is created
#' by \code{create.bootstrap.data} using the default for argument \code{rename}. 
#' The sum of the estimated number of individuals for each segment from 
#' \code{dis.data$NHAT} is given in the column \code{NHAT} in the new data. 
#' All other columns from the observation layer should be discarded. This is achieved 
#' by specifying the columns that should be retained using the argument \code{column.numbers}. 
#' Generally, all columns from the segment and higher levels should be kept. 
#' If the default is used, \code{column.numbers=NULL}, the columns \code{distance}, \code{object}, 
#' \code{size}, \code{distbegin} and \code{distend} from the observation level are automatically 
#' discarded. Note that for those columns from the observation layer that are kept, only the 
#' first recorded value will be transferred. 
#'
#' @param dis.data Data frame containing distance data (one row for each detection). Expects a column \code{NHAT}, i.e. size of detection divided by its probability of detection (see \code{create.NHAT}) and that and that ids in \code{segment.id2} are unique regardless of what resampled transect they belong to.
#' @param column.numbers Optional argument: vector of integers indicating which columns other than \code{NHAT} from \code{dis.data} should be retained in the returned data. 
#' 
#' @return This function returns bootstrapped count data that is suited for second stage count modelling of distance sampling 
#' data. The data includes the columns \code{NHAT} and \code{area} which are the response and
#' the offset required by functions concerned with second stage modelling from this package. 
#' 
#' @examples
#' data(dis.data.re)
#' # bootstrap data without stratification
#' dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
#' result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
#'              meta.data=list(width=250))
#' dis.data.re<-create.NHAT(dis.data.re,result)
#' 
#' bootstrap.data<-create.bootstrap.data(dis.data.re) 
#' 
#' bootcount.data<-create.bootcount.data(bootstrap.data)
#' 
#' @export
#' 

create.bootcount.data<-function(dis.data,column.numbers=NULL){
  print.first<-function(x){paste(x[1])}
  
  if(is.null(column.numbers)==T){
    obs.cols<-match(c("object","distance","size","distbegin","distend","detected","observer","NHAT"),names(dis.data))
    obs.cols<-obs.cols[which(is.na(obs.cols)==F)]
    column.numbers<-c(1:length(names(dis.data)))[-obs.cols]
    }
  # make sure that column 'segment.id2' is transferred over
  if(is.na(match("segment.id2",names(dis.data[,column.numbers])))){column.numbers<-c(column.numbers,which(names(dis.data)=="segment.id2"))}
  
  dis.data$NHAT[which(is.na(dis.data$NHAT)==T)]<-0
  glm.data<-dis.data[0,]
  unique.seg<-sort(unique(dis.data$segment.id2))
  glm.data[1:length(unique.seg),]<-NA
  # making sure that 'area' gets transferred into count.data
  if(is.na(match("area",names(dis.data)[column.numbers]))){column.numbers<-c(column.numbers,match("area",names(dis.data)))}
  for (s in column.numbers){
    glm.data[,s]    <- aggregate(dis.data[,s] ~ dis.data$segment.id2, FUN = print.first)[,2]
  }
  glm.data$NHAT<-aggregate(dis.data$NHAT ~ dis.data$segment.id2,FUN=sum)[,2]
  # converting these back to numbers
  for (n in column.numbers){
    if(is.numeric(dis.data[,n])){glm.data[,n]<-as.numeric(glm.data[,n])}
  }
  
  glm.data<-glm.data[,c(column.numbers,which(names(glm.data)=="NHAT"))]
  glm.data
}


