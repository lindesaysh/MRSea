#' Create bootstrap data for non-parametric bootstrapping
#'
#' This function creates one realisation of bootstrapped data based on \code{dis.data}. 
#' The default resampling unit is \code{transect.id} which may be modified using the argument \code{resample}. 
#' 
#' @param dis.data Original data to be bootstrapped. Requires a column that matches argument \code{resample} exactly. 
#' @param resample Specifies the resampling unit for bootstrapping, default is \code{transect.id}. Must match a column name in \code{dis.data} exactly
#' @param rename A vector of column names for which a new column needs to be created for the bootstrapped data. This defaults to \code{segment.id} for line transects, however others might be added
#' A new column with new ids will automatically be created for the column listed in \code{resample}
#' @param stratum The column name in \code{dis.data} that identifies the different strata. The default \code{NULL} returns un-stratified bootstrap data. If stratum is specified, this requires a column in \code{dis.data} that matches argument \code{stratum} exactly 
#' 
#' @return
#' Returns one realisation of bootstrapped distance data. Note that a new column 
#' (in addition to those listed under argument \code{rename}) is created. If the default for \code{resample} is used, 
#' a column with new unique ids called \code{transect.id2}. 
#' Note that a new column is created with renamed bootstrap resamples to preserve the number of unique bootstrap resamples. 
#' If the default for \code{resample} is used, i.e. \code{transect.id}, this new column is called \code{transect.id2}. 
#' In addition, a new column \code{segment.id2} is created which is required for other bootstrap functions. 
#' 
#' @examples
#' data(dis.data.re)
#' # run distance analysis to create NHATS
#' dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
#' result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
#'              meta.data=list(width=250))
#' dis.data.re<-create.NHAT(dis.data.re,result)
#' 
#' # bootstrap data without stratification
#' bootstrap.data<-create.bootstrap.data(dis.data.re) 
#' # boostrap data with stratification (here by survey which is composed of 
#' # season and impact)
#' dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
#' bootstrap.data.str<-create.bootstrap.data(dis.data.re, stratum = "survey.id") 
#' 
#' @export

create.bootstrap.data<-function(dis.data,resample = "transect.id",rename = "segment.id", stratum=NULL){
if(is.null(stratum)==T){
  unique.resamples<-unique(dis.data[,resample])
  resamples.no<-length(unique.resamples)   
  # resampling with replacement
  new.resamples<-sample(unique.resamples, size = resamples.no, replace = T)
  # creating the bootstrap data
  boot.data<-bootstrap.orig.data(dis.data,resample=resample,new.resamples=new.resamples,resamples.no=resamples.no)
  # changing the object names to be unique as required by the \code{mrds::ddf} function
  boot.data$object<-NA
  obj.na<-which(is.na(boot.data$distance)==F)
  boot.data$object[obj.na]<-1:length(obj.na) # ddf requires all object ids to be unique for method = "ds"
  # creating new columns listed in \code{rename}
  if(is.null(rename)==F){for (i in 1:length(rename)){
    new.name<-paste(rename[i],2,sep="")
     if(!is.na(match(rename[i],names(dis.data)))){
      boot.data[,new.name]<-paste(boot.data[,paste(resample,2,sep="")],boot.data[,rename[i]],sep="_")}
}
}
}
else{
  uniqueStrata<-unique(dis.data[,stratum])
  stratumNo<-length(uniqueStrata)
  boot.data<-dis.data[0,]
  for (str in 1:stratumNo){
  # the original data containing only rows pertaining to stratum 'str'
  dis.data.str<-dis.data[dis.data[,stratum]==uniqueStrata[str],]
  # the unique resamples from stratum str
  unique.resamples<-unique(dis.data.str[,resample])
  resamples.no<-length(unique.resamples)   
  new.resamples<-sample(unique.resamples,resamples.no,replace=T)
  new.dis.data.str<-bootstrap.orig.data(dis.data.str,resample=resample,new.resamples=new.resamples,resamples.no=resamples.no)
  boot.data<-rbind(boot.data,new.dis.data.str)
  }
  boot.data$object<-NA
  obj.na<-which(is.na(boot.data$distance)==F)
  boot.data$object[obj.na]<-1:length(obj.na) # ddf requires all object ids to be unique for method = "ds"
  if(is.null(rename)==F){for (i in 1:length(rename)){
    new.name<-paste(rename[i],2,sep="")
    if(!is.na(match(rename[i],names(dis.data)))==T){
      boot.data[,new.name]<-paste(boot.data[,paste(resample,2,sep="")],boot.data[,rename[i]],sep="_")}
  }
  }
}
boot.data
}
