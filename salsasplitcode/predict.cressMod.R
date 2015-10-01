#' Function for making predictions for a model containing a CReSS basis (two dimensional local smooth).
#' 
#' This function calculates vector of predictions on the scale of the response or link.
#' 
#' @param predict.data Data frame of covariate values to make predictions to
#' @param splineParams spline parameter object that describes the fitting of 2D and 1D splines in the model object
#' @param g2k Matrix of distances between prediction locations and knot locations (n x k). May be Euclidean or geodesic distances.
#' @param model Object from a GEE or GLM model
#' @param type Type of predictions required. (default=`response`, may also use `link`.
#' @param coeff Vector of coefficients (default = NULL). To be used when bootstrapping and sampling coefficients from a distribution e.g. in \code{do.bootstrap.cress}.
#' 
#' @details
#' Calculate predictions for a model whilst centering the CReSS bases in the same way as the fitted model. Note, if there is an offset in the model it must be called 'area'.
#' 
#' @return
#' Returns a vector of predictions on either the response or link scale
#' 
#' @examples
#' 
#' # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # offshore redistribution data
#' # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' data(dis.data.re)
#' data(predict.data.re)
#' data(knotgrid.off)

#' # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # distance sampling
#' dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
#' result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
#'         meta.data=list(width=250))
#' dis.data.re<-create.NHAT(dis.data.re,result)
#' count.data<-create.count.data(dis.data.re)
#' 
#' # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # spatial modelling
#' splineParams<-makesplineParams(data=count.data, varlist=c('depth'))
#' #set some input info for SALSA
#' count.data$response<- count.data$NHAT
#' # make distance matrices for datatoknots and knottoknots
#' distMats<-makeDists(cbind(count.data$x.pos, count.data$y.pos), na.omit(knotgrid.off))
#' # choose sequence of radii
#' r_seq<-getRadiiChoices(8,distMats$dataDist)
#' # set initial model without the spatial term
#' initialModel<- glm(response ~ as.factor(season) + as.factor(impact) + offset(log(area)),  
#'                 family='quasipoisson', data=count.data)
#' # make parameter set for running salsa2d
#' salsa2dlist<-list(fitnessMeasure = 'QICb', knotgrid = knotgrid.off, knotdim=c(26,14), startKnots=4, minKnots=4, 
#'                  maxKnots=20, r_seq=r_seq, gap=4000, interactionTerm="as.factor(impact)")
#' salsa2dOutput_k6<-runSALSA2D(initialModel, salsa2dlist, d2k=distMats$dataDist, 
#'                    k2k=distMats$knotDist, splineParams=splineParams) 
#' 
#' splineParams<-salsa2dOutput_k6$splineParams
#' # specify parameters for local radial function:
#' radiusIndices <- splineParams[[1]]$radiusIndices
#' dists <- splineParams[[1]]$dist
#' radii <- splineParams[[1]]$radii
#' aR <- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]
#' count.data$blockid<-paste(count.data$transect.id, count.data$season, count.data$impact, sep='')
#' # Re-fit the chosen model as a GEE (based on SALSA knot placement) and GEE p-values
#' geeModel<- geeglm(formula(salsa2dOutput_k6$bestModel), data=count.data, family=poisson, id=blockid)
#' dists<-makeDists(cbind(predict.data.re$x.pos, predict.data.re$y.pos), na.omit(knotgrid.off), 
#'        knotmat=FALSE)$dataDist
#'        
#' # make predictions on response scale
#' preds<-predict.cress(predict.data.re, splineParams, dists, geeModel)
#' 
#' @export
#' 

predict.cress<-function(predict.data, splineParams, g2k, model, type='response', coeff=NULL, modelav=TRUE){
  
  attributes(model$formula)$.Environment<-environment()
  radii<-splineParams[[1]]$radii
  radiusIndices<-splineParams[[1]]$radiusIndices
  
  if(modelav==TRUE){
  aR<- splineParams[[1]]$knotPos
  }
  if(modelav==FALSE){
    aR<- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]
      }
  
  #print(aR)
  dists<- g2k
  x2<- data.frame(response=rpois(nrow(predict.data),lambda = 5), predict.data)
  #fakemodel<- eval(parse(text=paste("glm(model$formula, data=x2, family=", model$family$family,"(link='", model$family$link, "'))", sep='')))
  fakemodel<- glm(model$formula, data=x2, family=model$family)
  modmat<- model.matrix(fakemodel)
  
#   
#   # get centering of fakeglm above
#   bad_cntr<-attributes(fakemodel$model$'LocalRadialFunction(radiusIndices, dists, radii, aR)')$`scaled:center`
#   
#   # get centering of original model
#   cntr<- attributes(model$model$'LocalRadialFunction(radiusIndices, dists, radii, aR)')$`scaled:center`
#   
#   # find columns with interaction in them
#   intcol<- grep(':', attributes(modmat)$dimnames[[2]])
#   
#   # find columns with radial in them but not interaction (if there)
#   if(length(intcol)==0){
#     radcol<- grep('Local', attributes(modmat)$dimnames[[2]])
#   }else{
#     radcol<- grep('Local', attributes(modmat)$dimnames[[2]][-intcol])  
#   }
#   
#   # take local radial model matrix and add bad centre back on
#   
#   bad_cntr.mat<-matrix(rep(bad_cntr, each = nrow(modmat)), ncol=length(bad_cntr), byrow=F)
#   #modmat[,radcol]<-modmat[,radcol]+bad_cntr.mat
#   #head(modmat[,6:9])
#   # now re-centre with centering from original fitted model
#   cntr.mat<-matrix(rep(cntr, each = nrow(modmat)), ncol=length(cntr), byrow=F)
#   modmat[,radcol]<-modmat[,radcol]+bad_cntr.mat-cntr.mat
#   #head(modmat[,6:9])
# 
#   # do the same for the interaction if present
#   if(length(intcol)!=0){
#     modmat[,intcol]<-modmat[,intcol]+bad_cntr.mat-cntr.mat
#   }
#   
  if(is.null(coeff)){
    modcoef<- as.vector(model$coefficients)  
  }else{
    modcoef<-coeff
  }
  if(type=='response'){
    if(length(model$offset)>0 & sum(model$offset)!=0){
      preds<- model$family$linkinv(modmat%*%modcoef)* predict.data$area  
    }else{
      preds<- model$family$linkinv(modmat%*%modcoef)
    }
  }
  
  if(type=='link'){
    preds<- modmat%*%modcoef
    print('warning: no offset included as link response specified')
  }
  
  return(preds)
}
