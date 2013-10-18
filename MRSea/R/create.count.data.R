#' Aggregate distance data into count data
#'
#' This function creates a new data set where \code{dis.data} is aggregated for 
#' each visit to a segment (\code{segment.id}). The sum of the estimated number 
#' of individuals for each segment from \code{dis.data$NHAT} is given in the 
#' column \code{NHAT} in the new data. 
#' Only columns from the segment or higher layers should be carried over into 
#' \code{count.data} from \code{dis.data}. Use argument \code{column.numbers} to 
#' identify these.
#'
#' @return This function returns count data that is suited for second stage count modelling of distance sampling 
#' data. The data includes the columns \code{NHAT} and \code{area} which are the response and
#' the offset required by functions concerned with second stage modelling from this package. 
#' 
#' @param dis.data Data frame containing distance data (one row for each detection). Expects a column \code{NHAT}, i.e. size of detection divided by its probability of detection (see \code{create.NHAT}) and that ids in \code{segment.id} are unique regardless of what transect they belong to 
#' @param column.numbers Optional argument: vector of integers indicating which columns other than \code{NHAT} from \code{dis.data} should be retained in the returned data. Generally all columns from the segment and higher levels should be kept while those from the observation level should be discarded. If the default is used, \code{column.numbers=NULL}, the columns \code{distance}, \code{object}, \code{size}, \code{distbegin} and \code{distend} from the observation level are automatically discarded. Note that for those columns from the observation layer that are kept, only the first recorded value will be transferred. 
#' 
#' @examples
#' data(dis.data.re)
#' dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
#' result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
#'            meta.data=list(width=250))
#' dis.data.re<-create.NHAT(dis.data.re,result)
#' count.data<-create.count.data(dis.data.re)
#' 
#' @export

create.count.data<-function(dis.data,column.numbers=NULL){
print.first<-function(x){paste(x[1])}

if(is.null(column.numbers)==T){
  obs.cols<-match(c("object","distance","size","distbegin","distend","detected","observer","NHAT"),names(dis.data))
  obs.cols<-obs.cols[which(is.na(obs.cols)==F)]
  column.numbers<-c(1:length(names(dis.data)))[-obs.cols]
}
dis.data$NHAT[which(is.na(dis.data$NHAT)==T)]<-0
glm.data<-dis.data[0,]
unique.seg<-sort(unique(dis.data$segment.id))
glm.data[1:length(unique.seg),]<-NA
# glm.data$segment.id<-unique.seg
# making sure that 'area' gets transferred into count.data
if(is.na(match("area",names(dis.data)[column.numbers]))){column.numbers<-c(column.numbers,match("area",names(dis.data)))}
for (s in column.numbers){
  print(names(dis.data)[s])
  glm.data[,s]    <- aggregate(dis.data[,s] ~ dis.data$segment.id, FUN = print.first)[,2]
  }
  glm.data$NHAT<-aggregate(dis.data$NHAT ~ dis.data$segment.id,FUN=sum)[,2]
  # converting these back to numbers
  for (n in column.numbers){
  if(is.numeric(dis.data[,n])){glm.data[,n]<-as.numeric(glm.data[,n])}
  }

glm.data<-glm.data[,c(column.numbers,which(names(glm.data)=="NHAT"))]
glm.data
}
