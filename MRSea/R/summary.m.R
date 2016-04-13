#' Summary function for MRSea models to shorten names.
#' 
#' @param model Fitted model object (glm or geeglm)
#' @param varlist (\code{default=NULL}). Vector of names of continuous variables.
#' @param facls (\code{default=NULL}). Vector or names of factor variables.
#' 
#' @return
#' summary object as for \code{summary.glm} but with the variable names shortened for ease of viewing. 
#' 
#' @examples
#' # load data
#' data(ns.data.re)
#' 
#' splineParams<-makesplineParams(data=ns.data.re, varlist=c('observationhour'))
#' 
#' # make column with foldid for cross validation calculation
#' ns.data.re$foldid<-getCVids(ns.data.re, folds=5)
#' 
#' #' # set initial model without the spline terms in there 
#' # (so all other non-spline terms)
#' ns.data.re$response<- ns.data.re$birds
#' initialModel<- glm(response ~ as.factor(floodebb) + as.factor(impact) + offset(log(area)), family='quasipoisson',data=ns.data.re)
#' 
#' #set some input info for SALSA
#' salsa1dlist<-list(fitnessMeasure = 'QBIC', minKnots_1d=c(1), maxKnots_1d = c(5), startKnots_1d = c(2), degree=c(2), maxIterations = 10, gaps=c(0))
#' 
#' # run SALSA
#' salsa1dOutput<-runSALSA1D_withremoval(initialModel, salsa1dlist, varlist=c('observationhour'), factorlist=c('floodebb', 'impact'), ns.predict.data.re, splineParams=splineParams, datain=ns.data.re)
#' 
#' summary.m(salsa1dOutput$bestModel,varlist=c('observationhour'), facls=c('floodebb', 'impact'))
#' 
#' # compare output with:
#' summary(salsa1dOutput$bestModel)
#' 
#' @export
#' 
summary.m<-function(model, varlist=NULL, facls=NULL){
  bob<-attr(model$coefficients, 'names')
  for(i in 1:length(varlist)){
    id<-grep(varlist[i], bob)
    if(length(id)>1){
    for(j in 1:length(id)){
      bob[id][j]<-paste('s(', varlist[i],')', j, sep='') 
    }
    }else{
      bob[id]<-paste(varlist[i], sep='')
    }
  }
  # check for interaction terms
  localid<-grep('LocalRadial', bob)
  if(length(localid>1)){
    intid<-grep(':', bob)
    # change smooth term (without interaction)
    smoothid<- localid[which(is.na(match(localid, intid)))]
    for(k in 1:length(smoothid)){
      bob[smoothid][k]<-paste('s(x.pos, y.pos)b', k, sep='')
    }
    # change interaction terms:
    for(i in 1:length(facls)){
      a<-grep(facls[i], bob[intid])  
      if(length(a)>0){
        a<-i
        break
      }
    }
    
    counter<-1
    for(k in 1:length(smoothid)){
        for(i in 1:(length(intid)/length(smoothid))){
          #print(c(b, k, i))
          textin<-paste('as.factor(',facls[a],')', i, ':s(x.pos, y.pos)b', k, sep='')
          bob[intid[counter]]<-textin
         #print(textin)
          counter<-counter+1
        }  
      }
  }
  attr(model$coefficients, 'names')<-bob
  if(class(model)[1]=='geeglm'){
    attr(model$geese$beta, 'names') <- bob
  }
  summary(model)
}






