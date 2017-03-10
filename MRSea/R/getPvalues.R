#-----------------------------------------------------------------------------
#' Calculate marginal p-values from a \code{model}. 
#' 
#' An ANOVA is fitted repeatedly with each covariate being the last so that the output is marginal.  \code{varlist} and \code{factorlist} are optional and shorten the variable names in the output.
#' 
#' @param model Fitted model object of class \code{gee}.
#' @param varlist (default =\code{NULL}). Vector of covariate names (continous covariates only) used to make the output table names shorter.  Useful if spline parameters are specified in the model.
#' @param factorlist (default =\code{NULL}). Vector of covariate names (factor covariates only) used to make the output table names shorter. Useful if spline parameters are specified in the model.
#' 
#' @return
#' Print out table of each variable and its associated marginal p-value.
#' 
#' @examples 
#' 
#' # load data
#' data(ns.data.re)
#' 
#' # make blocking structure
#' ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear, 
#'                     ns.data.re$DayOfMonth, sep='')
#' ns.data.re$blockid<-as.factor(ns.data.re$blockid)
#' 
#' initialModel<- geeglm(birds ~ as.factor(floodebb) + as.factor(impact) + observationhour + x.pos + 
#'               y.pos + offset(log(area)), family='poisson',data=ns.data.re, id=blockid)
#' 
#' getPvalues(initialModel, varlist=c('observationhour', 'x.pos', 'y.pos'), 
#'             factorlist=c('floodebb', 'impact'))
#' 
#' getPvalues(initialModel)
#' 
#' @export
#' 
getPvalues<-function(model, varlist=NULL, factorlist=NULL){
  
  if(class(model)[1]!='geeglm'){stop('Class of model not geeglm. Fit model as GEE or use an appropriate ANOVA function (e.g. Anova or anova)')}
  
  print("Getting marginal p-values")
  
  # make list of terms if varlist and factorlist are spcified
  if(is.null(varlist)==FALSE){
    termlist<-rep(NA, length=length(labels(terms(model))))
    cressterm<-grep('LocalRadial', labels(terms(model)))
    
    if(length(varlist)+length(factorlist) + length(cressterm)!=length(termlist)) stop('not all one dimensional terms specified in varlist or factorlist')
    
    if(length(varlist)>0){
      for(v in 1:length(varlist)){
        termlist[grep(varlist[v], labels(terms(model)))]<- varlist[v]  
      }
    }
    
    if(length(factorlist)>0){
      for(f in 1:length(factorlist)){
        termlist[grep(factorlist[f], labels(terms(model)))]<- factorlist[f]  
      }
    }
    
    if(length(cressterm)>0){
      for(a in 1:length(cressterm)){
        if(is.na(termlist[cressterm[a]])==FALSE){
          termlist[cressterm[a]]<- paste('s(x.pos, y.pos):', termlist[cressterm[a]], sep='')
        }else{termlist[cressterm[a]]<- 's(x.pos, y.pos)'  }
      }  
    }
    
  }
  
  
  #get marginal p-values for each term:
  store<- matrix(0, ncol=2, nrow=length(labels(terms(model))))
  for(i in 1:length(labels(terms(model)))){
    covariate<- labels(terms(model))[i]  
    
    text1<-paste("anova(update(model, . ~ . - ",  covariate, "+",  covariate, "))", sep="")
    test<-eval(parse(text=text1))
    
    if(is.null(varlist)){
      store[i,1] <- covariate
    }else{      
      store[i,1] <- termlist[i]
    }
    
    store[i,2]<- round(test[which(labels(test)[[1]]==covariate),3],6)
    if(as.numeric(store[i,2])<0.0001){store[i,2]<-'<0.0001'}
  }
  
  # check for zero p-values
  which(store[,2]==0)
  #print results with p-values
  store<- as.data.frame(store)
  names(store)<-c("Variable", "p-value")
  print(store)
  
}
