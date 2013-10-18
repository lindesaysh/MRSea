#' Obtaining a data frame of bootstrapped data using resamples
#' 
#' This function extracts the records corresponding to each resample from the original 
#' distance data and pastes them together in a new data frame which is returned. 
#'
#' @param orig.data Original data to be bootstrapped
#' @param resample Specifies the resampling unit for bootstrapping, default is \code{transect.id}. Must match a column name in \code{orig.data} exactly
#' @param new.resamples String of resampled units from \code{data[,"resample"]}. Created by \code{create.bootstrap.data()}
#' @param resamples.no Length of new.resamples
#'
#' @return
#' Returns bootstrapped data. Internal function called by function \code{create.bootstrap.data}.
#'
#' @examples
#' data(dis.data.re)
#' resample<-"transect.id"
#' samples<-unique(dis.data.re[,resample])
#' resamples.no<-length(samples)
#' new.resamples<-sample(samples,resamples.no,replace=TRUE)
#' bootstrap.data<-bootstrap.orig.data(dis.data.re,resample,new.resamples,resamples.no)
#' 
#' @export

bootstrap.orig.data<-function(orig.data,resample,new.resamples,resamples.no){
  
# giving the resamples a new unique name
new.resamples2<-paste(new.resamples,1:resamples.no,sep="_")

# creating the new column with the renamed resamples
new.name<-paste(resample,2,sep="")
orig.data[,new.name]<-NA
  
### creating the datafile:
# get the data for each new.resamples from the original data
# starting the new data set with the same data frame structure as orig.data
new.data<-orig.data[0,]
# and adding the data for each of the resamples             
# giving the resamplers new names (important for random effects models to keep the number of groups the same as original data)     
for (m in 1:resamples.no){
sample.data<-orig.data[!is.na(match(orig.data[,resample],new.resamples[m])),]
# this is a catch in case \code{orig.data} does not contain any records for any of the resamples. This may be the case if resamples were obtained from count.data containing effort data without detections and distance data is bootstrapped using these resamples.
if(length(sample.data[,new.name])>0){
sample.data[,new.name]<-new.resamples2[m] 
new.data<-rbind(new.data,sample.data)
}
}
new.data
}
        




