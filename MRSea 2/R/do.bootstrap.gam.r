#' Bootstrapping function without model selection using \code{gam} as the second stage count model
#' 
#' This fuction performs a specified number of bootstrapping iterations using gams for fitting the 
#' second stage count model. See below for details. 
#' 
#' @param orig.data The original data. In case \code{ddf.obj} is specified, this should be the original distance data. In case \code{ddf.obj} is \code{NULL}, it should have the format equivalent to \code{count.data} where each record represents the summed up counts at the segments. 
#' @param predict.data The prediction grid data 
#' @param ddf.obj The ddf object created for the best fitting detection model. Defaults to \code{NULL} for nearshore data. 
#' @param model.obj The best fitting \code{gam} model for the original count data
#' @param resample Specifies the resampling unit for bootstrapping, default is \code{transect.id}. Must match a column name in \code{dis.data} exactly
#' @param rename A vector of column names for which a new column needs to be created for the bootstrapped data. This defaults to \code{segment.id} for line transects (which is required for \code{create.bootcount.data}), others might be added. 
#' A new column with new ids will automatically be created for the column listed in \code{resample}. In case of nearshore data, this argument is ignored. 
#' @param stratum The column name in \code{orig.data} that identifies the different strata. The default \code{NULL} returns un-stratified bootstrap data. In case of nearshore data, this argument is ignored.  
#' @param B Number of bootstrap iterations
#' @param name Analysis name. Required to avoid overwriting previous bootstrap results. This name is added at the beginning of "predictionboot.RData" when saving bootstrap predictions. 
#' @param save.data If TRUE, all created bootstrap data will be saved as an RData object in the working directory at each iteration, defaults to FALSE
#' @param nhats (default = FALSE). If you have calculated bootstrap NHATS because there is no simple ddf object then a matrix of these may be fed into the function.  The number of columns of data should >= B.  The rows must be equal to those in \code{orig.data} and \code{d2k} and \emph{must} be in matching order.
#' 
#' @details
#' In case of distance sampling data, the following steps are performed for each iteration: 
#' 
#'   - the original data is bootstrapped
#' 
#'   - a detection function is fitted to the bootstrapped data
#'   
#'   - a count model is fitted to the bootstrapped data
#'   
#'   - coefficients are resampled from a multivariate normal distribution defined by MLE and COV from count model
#'   
#'   - predictions are made to the prediction data using the resampled coefficients 
#'   
#' In case of count data, the following steps are performed for each iteration
#'   
#'   - coefficients are resampled from a multivariate normal distribution defined by MLE and COV from the best fitting count model
#'   
#'   - predictions are made to the prediction data using the resampled coefficients 
#' 
#' @return
#' The function returns a matrix of bootstrap predictions. The number of rows is equal to the number of rows in predict.data.  The number of columns is equal to \code{B}.  The matrix may be very large and so is stored directly into the working directory as a workspace object: '"name"predictionboot.RObj'.  The object inside is called \code{bootPreds}.
#' 
#' @examples
#' # offshore redistribution data
#' data(dis.data.re)
#' data(predict.data.re)
#' dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
#' result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
#'             meta.data=list(width=250))
#' dis.data.re<-create.NHAT(dis.data.re,result)
#' count.data<-create.count.data(dis.data.re)
#' require(mgcv)
#' gam.2<-gam(NHAT~as.factor(impact)+s(x.pos,y.pos,by=as.factor(impact))+offset(log(area)),
#'            data=count.data,family=quasipoisson)
#' do.bootstrap.gam(dis.data.re,predict.data.re,ddf.obj=result,model.obj=gam.2,resample="transect.id",
#'                rename="segment.id",stratum='survey.id',1,name='gam',save.data=FALSE,nhats=NULL)
#' load("gampredictionboot.RData") # loading the predictions into the workspace
#' # look at the first 6 lines of the predictions on the response scale
#' head(bootPreds)
#' 
#' 
#' \dontrun{# nearshore redistribution data
#' data(ns.data.re)
#' data(ns.predict.data.re)
#' require(mgcv)
#' gam.ns2=gam(birds~as.factor(impact)+s(x.pos,y.pos,by=as.factor(impact))+offset(log(area)),
#'          data=ns.data.re,family=quasipoisson)
#' do.bootstrap.gam(ns.data.re,ns.predict.data.re,ddf.obj=NULL,model.obj=gam.ns2,resample=NULL,
#'                 rename=NULL,stratum=NULL,1,name='ns.gam',save.data=FALSE,nhats=NULL)
#' # load the replicate predictions into the workspace               
#' load("ns.gampredictionboot.RData") 
#' # look at the first 6 lines of the predictions on the response scale
#' head(bootPreds)}
#' 
#' @export
#' 
#' 
do.bootstrap.gam<-function(orig.data,predict.data,ddf.obj=NULL,model.obj,resample="transect.id",rename="segment.id",stratum=NULL,B,name=NULL,save.data=FALSE,nhats=NULL){
  require(mvtnorm)
  if(is.null(ddf.obj)==F){
  # fixing the settings for the detection model
  key <- ddf.obj$ds$aux$ddfobj$type
  scale.formula<-ddf.obj$ds$aux$ddfobj$scale$formula
  point<-ddf.obj$ds$aux$point
  binned<-ddf.obj$ds$aux$binned
  breaks<-ddf.obj$ds$aux$breaks
  width<-ddf.obj$ds$aux$width
  my.formula=eval(parse(file="",text=paste("~mcds(key = '",key,"', formula=",scale.formula,")",sep="")))
  #result=ddf(dsmodel = as.formula(my.formula), data = data$dis.object, method = "ds", meta.data = list(point = point, width = data$w, binned = data$binned))
  }
  # object for storing the predictions
  bootPreds<- matrix(NA, nrow=nrow(predict.data), ncol=B)
  
  for (b in 1:B){
  if(is.null(ddf.obj)==F){
  # create bootstrap data
  bootstrap.data<-create.bootstrap.data(orig.data, stratum = stratum)
    
  if (save.data==T){
    filename=paste("bootstrap.data_",b,".RData",sep="")
    save(bootstrap.data,file=filename)}
  
  # fit the same detection model as in ddf.obj to the bootstrap data
  result.boot<-ddf(dsmodel=as.formula(my.formula), data=bootstrap.data, method="ds", meta.data=list(width=width, binned=binned, breaks=breaks))
  
  # create count data
  bootstrap.data<-create.NHAT(bootstrap.data,result.boot)
  boot.count.data<-create.bootcount.data(bootstrap.data)
  }else{boot.count.data<-orig.data}  
  
  # fit count model
  boot.gam<-update(model.obj, .~., data=boot.count.data)

  # sample from multivariate normal
  varcov<- vcov(boot.gam)
  ests<- coefficients(boot.gam)
  samplecoeff<- as.numeric(rmvnorm(1,ests,varcov))
  
  # make predictions
  prediction<- predict(boot.gam, predict.data, type='lpmatrix')%*%samplecoeff +log(predict.data$area)
  
  
  bootPreds[,b]<- exp(prediction)
  } # end of for loop
  
# save predictions  -  and other data? e.g. parameter values?
save(bootPreds, file=paste(name,"predictionboot.RData",sep=""), compress='bzip2')
  
} # end of function
  
  