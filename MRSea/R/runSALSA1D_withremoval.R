
#' Running SALSA for continuous one-dimensional covariates.
#' 
#' This function finds spatially adaptive knot locations for one or more continuous one-dimensional covariates.  It differs to \code{\link{runSALSA1D}} in that if the CV score of a model does not improve with the addition of a covariate in \code{varlist} then that term is either reduced to linear or removed from the model.
#' 
#' @param initialModel The best fitting \code{CReSS} model with no continuous covariates specified.  This must be a model of class \code{glm}.
#' @param salsa1dlist Vector of objects required for \code{runSALSA1D}: \code{fitnessMeasure}, \code{minKnots_1d}, \code{maxKnots_1d}, \code{startKnots_1d} \code{degree}, \code{maxIterations} \code{gap}. 
#' @param varlist Vector of variable names for the covariates required for knot selection
#' @param factorlist vector of factor variables specified in \code{initialModel}.  Specified so that a check can be made that there are non-zero counts in all levels of each factor. Uses the function \code{checkfactorlevelcounts}. Default setting is NULL.
#' @param predictionData The data for which predictions are to be made. Column names must correspond to the data in \code{initialModel}. If predictionData is not specified (\code{NULL}), then the range of the data is used to create the smooth terms.
#' @param varlist_cyclicSplines Vector of variable names for covariates to be modelled with cyclic cubic splines.  This must be a subset of \code{varlist}.The default is \code{NULL}
#' @param splineParams List object containing information for fitting splines to the covariates in \code{varlist}. If not specified (\code{NULL}) this object is created and returned. See \code{\link{makesplineParams}} for details.
#' @param datain Data used to fit the initial Model.
#' @param removal (Default: \code{TRUE}). Logical stating whether a selection procedure should be done to choose smooth, linear or removal of covariates.  If \code{FALSE} all covariates are returned and smooth.
#' 
#' @details
#' There must be columns called \code{response} (response variable) and \code{foldid} (for cross-validation calculation) in the data used in the initial model to be fitted. If the data is proportion, then there should be two columns called \code{successess} and \code{failures}.
#' 
#' The object \code{salsa1dlist} contains parameters for the \code{runSALSA1D} function.
#'   
#' \code{fitnessMeasure}. The criterion for selecting the `best' model.  Available options: AIC, AIC_c, BIC, QIC_b.
#' 
#'    \code{minKnots_1d}. Minimum number of knots to be tried.
#'    
#'    \code{maxKnots_1d}. Maximum number of knots to be tried.
#'    
#'    \code{startKnots_1d}. Starting number of knots (spaced at quantiles of the data).
#'    
#'    \code{degree}. The degree of the B-spline. Does not need to be specified if \code{splineParams} is a parameter in \code{runSALSA1D}.
#'    
#'    \code{maxIterations}.The exchange/improve steps will terminate after maxIterations if still running.
#'    
#'    \code{gaps}. The minimum gap between knots (in unit of measurement of explanatory).
#'    
#'    
#' \code{minKnots_1d}, \code{maxKnots_1d}, \code{startKnots_1d} and \code{gaps} are vectors the same length as \code{varlist}.  This enables differing values of these parameters for each covariate.
#'
#'  The initial model contains all the factor level covariates and any covariates of interest that are not specified in the \code{varlist} argument of \code{runSALSA1D} 
#' 
#' \emph{Note:} The algorithm may remove variables in \code{varlist} but not the variables in \code{factorlist}.  If there is no better model than with a knot at the mean, the output will include that covariate with a knot at the mean.  The best model with a given smooth term is tested both against a model with the term as linear or removed. Cross-Validation is used in the selection process.
#' 
#' @return
#' A list object is returned containing 4 elements:
#' 
#' \item{bestModel}{A glm model object from the best model fitted}
#' \item{modelFits1D}{A list object with an element for each new term fitted to the model.  The first element is a model fitted with a knot at the mean for each of the covariates (startmodel) in \code{varlist}.  Within the first element, the current fit and formula of the start model.  
#' 
#' The second element is the result of SALSA on the first term in \code{varlist}.  Within this element:
#' \itemize{
#' \item \code{term}: term of interest
#' \item \code{kept}: Statement of whether the term is kept in the model (yes- initial knots, yes - new knots, yes -linear or no)
#' \item \code{basemodelformula}: the resulting model formula.  If \code{kept=yes} or \code{kept=linear} then the term of interest is included in the model otherwise it is removed.
#' \item \code{knotSelected}: the knots chosen for the term of interest (NA if term removed or linear)
#' \item \code{baseModelFits}: fit statistics for the resulting formula
#' \item \code{modelfits}: fit statistics for the model with the term included (same as resulting formula if \code{kept=yes})
#' } 
#' 
#' This continues till all covariates in \code{varlist} have been through SALSA.}
#' \item{splineParams}{The updated spline parameter object, with the new (if chosen) knot locations for each covariate in \code{varlist}}
#' \item{fitstat}{The final fit statistic of \code{bestModel}.  The type of statistic was specified in \code{salsa1dlist}.}
#' \item{keptvarlist}{The covariates from \code{varlist} that have been retained in the model}
#' 
#' 
#' @references Walker, C.; M. Mackenzie, C. Donovan and M. O'Sullivan. SALSA - a Spatially Adaptive Local Smoothing Algorithm. Journal of Statistical Computation and Simulation, 81(2):179-191, 2010
#' 
#' @examples
#' # load data
#' data(ns.data.re)
#' # load prediction data
#' data(ns.predict.data.re)
#' 
#' splineParams<-makesplineParams(data=ns.data.re, varlist=c('observationhour', 'DayOfMonth'))
#' 
#' # make column with foldid for cross validation calculation
#' ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear, ns.data.re$DayOfMonth, sep='')
#' ns.data.re$blockid<-as.factor(ns.data.re$blockid)
#' ns.data.re$foldid<-getCVids(ns.data.re, folds=5, block='blockid')
#' 
#' #' # set initial model without the spline terms in there 
#' # (so all other non-spline terms)
#' ns.data.re$response<- ns.data.re$birds
#' initialModel<- glm(response ~ as.factor(floodebb) + as.factor(impact) + offset(log(area)), 
#'                     family='quasipoisson',data=ns.data.re)
#' 
#' #set some input info for SALSA
#' salsa1dlist<-list(fitnessMeasure = 'QBIC', minKnots_1d=c(2,2), maxKnots_1d = c(5, 5), 
#'                   startKnots_1d = c(2,2), degree=c(2,2), maxIterations = 10, gaps=c(1,1))
#' 
#' # run SALSA
#' salsa1dOutput<-runSALSA1D_withremoval(initialModel, salsa1dlist, varlist=c('observationhour', 'DayOfMonth'), 
#'                  factorlist=c('floodebb', 'impact'), ns.predict.data.re, splineParams=splineParams, datain=ns.data.re)
#' 
#' @export
#' 
runSALSA1D_withremoval<-function(initialModel, salsa1dlist, varlist, factorlist=NULL, predictionData=NULL, varlist_cyclicSplines=NULL, splineParams=NULL, datain, removal=TRUE, panelid=NULL){
  
  require(splines)
  require(fields)
  
  if(class(initialModel)[1]!='glm') stop('Class of model not supported.  Please use glm')
  
  # check for response variable
  if(nrow(initialModel$data)!=length(initialModel$model[[1]])){
    fam<-'BinProp'
  }else{
    fam<-'other'
  }
  
  if(is.null(factorlist)==F & fam=='other'){
    # check factor level counts:
    checkfactorlevelcounts(factorlist, initialModel$data, initialModel$y)
  }
  
  # ~~~~~~~~~~~~ SET UP ~~~~~~~~~~~~~~~~
  # set parameters for SALSA
  winHalfWidth = 0
  family<-initialModel$family$family
  link<-initialModel$family$link
  data<-initialModel$data
  if(sum(abs(dim(data)-dim(datain)))>0) stop('Data dimensions do not match the data in initialModel')
  
  attributes(initialModel$formula)$.Environment<-environment()

  
  if(fam=='other'){
    if(is.null(data$response)) stop('data does not contain response column')  
  }else{
    if(is.null(data$successes)) stop('data does not contain successes column')
    if(is.null(data$failures)) stop('data does not contain failures column')
  }
  
  
  # check parameters in salsa1dlist are same length as varlist
  if(length(varlist)!=length(salsa1dlist$minKnots_1d)) stop('salsa1dlist$minKnots_1d not same length as varlist')
  if(length(varlist)!=length(salsa1dlist$maxKnots_1d)) stop('salsa1dlist$maxKnots_1d not same length as varlist')
  if(length(varlist)!=length(salsa1dlist$startKnots_1d)) stop('salsa1dlist$startKnots_1d not same length as varlist')
  if(length(varlist)!=length(salsa1dlist$degree)) stop('salsa1dlist$degree not same length as varlist')
  if(length(varlist)!=length(salsa1dlist$gaps)) stop('salsa1dlist$gaps not same length as varlist')
  
  if(is.null(data$foldid)) stop('no column called "foldid" in data')
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(splineParams)){
    splineParams<-makesplineParams(data, varlist, predictionData, salsa1dlist$degree)
    varID<-(1:length(varlist))+1
  }else{
    # check whats in varlist and splineParams to see if they match
    if((length(splineParams)-1)!=length(varlist)){
      varID<-vector(length=length(varlist))
      for(i in 1:length(varlist)){
        varID[i]<-grep(varlist[i], splineParams)
      }
    }else{varID<-(1:length(varlist))+1}
  }
  
  terms1D <- list(length(varlist))
  
  counter<-1 # needed to loop through cyclics if more than one
    for(i in 2:(length(varlist)+1)){
    if(varlist[varID[(i-1)]-1]%in%varlist_cyclicSplines){
      require(mgcv)
      terms1D[[(i-1)]]<- paste("smooth.construct(s(", varlist[varID[(i-1)]-1], ", bs='cc', k=(length(splineParams[[", varID[(i-1)], "]]$knots))+2), knots = list(",varlist[varID[(i-1)]-1], "=as.numeric(c(splineParams[[",varID[(i-1)], "]]$bd[1], splineParams[[", varID[(i-1)], "]]$knots, splineParams[[",varID[(i-1)], "]]$bd[2]))), data=data.frame(",varlist[varID[(i-1)]-1],"))$X[,-1]", sep="")
      splineParams[[i]]$knots<-eval(parse(text=paste('quantile(data$', varlist_cyclicSplines[counter], ', probs = c(0.25, 0.5, 0.75))', sep='')))
      counter<-counter+1
    }else{
      terms1D[[(i-1)]]<- paste("bs(", varlist[(i-1)], ", knots = splineParams[[", varID[(i-1)], "]]$knots, degree=splineParams[[", varID[(i-1)], "]]$degree, Boundary.knots=splineParams[[",varID[(i-1)], "]]$bd)", sep='')
    }
  }
  
  splineParams<<-splineParams
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # initial model is all variables in varlist with knots at locations in splineParams and NO 2D smooth
  require(mgcv)
  if(fam=='BinProp'){
    baseModel <- eval(parse(text=paste("glm(cbind(successes,failures) ~ ", paste(formula(initialModel)[3],sep=""), "+", paste(terms1D, collapse="+"),", family =", family,"(link=", link,"), data = data)", sep='')))
  }else{
    baseModel <- eval(parse(text=paste("glm(response ~ ", paste(formula(initialModel)[3],sep=""), "+", paste(terms1D, collapse="+"),", family =", family,"(link=", link,"), data = data)", sep='')))
  }
  
  if(removal==TRUE){
    cv_initial <- getCV_CReSS(data, baseModel, splineParams=splineParams)  
  }else{
    cv_initial=NULL
  }
  
  # re-fit base model with max knots allowed to get best estimate of dispersion parameter
  if(salsa1dlist$fitnessMeasure=='QAIC' | salsa1dlist$fitnessMeasure=='QAICc' | salsa1dlist$fitnessMeasure=='QBIC'){
    splineParams_base<-splineParams
    for(i in 2:(length(varlist)+1)){
      eval(parse(text=paste("splineParams[[", i, "]]$knots<-seq(min(splineParams[[", i,"]]$explanatory),max(splineParams[[", i,"]]$explanatory),length=(salsa1dlist$maxKnots_1d[",(i-1),"])+2)[2:(salsa1dlist$maxKnots_1d[",(i-1),"]+1)]", sep='')))
    }
    splineParams<<-splineParams
    dispersion_Model<-update(baseModel, .~.)
    salsa1dlist$fitnessMeasure<-c(salsa1dlist$fitnessMeasure, summary(dispersion_Model)$dispersion)
    
    # return the splineParams object back to the original
    splineParams<<-splineParams_base
    if(dispersion_Model$conv==FALSE) stop('Model to get dispersion parameter, with max knots evenly spaced did not converge.  Use fewer max knots or change fitness measure.')
    
    # update the modelling to use the poisson family. The dispersion has been calculated and will be used by the get.measure function.
    if(family=='quasipoisson'){family='poisson'}
    if(family=='quasibinomial'){family='binomial'}
    baseModel<-eval(parse(text=paste("update(baseModel, .~., family=",substitute(family), "(link=", substitute(link),"))", sep='')))
  }
  
  initDisp<-getDispersion(baseModel)
  fitStat<-get.measure(salsa1dlist$fitnessMeasure,'NA',baseModel, initDisp)$fitStat
  
  modelFits1D <- list((length(varlist)+1))
  modelFits1D[[1]] <- list(term = 'startmodel', kept=NULL, basemodelformula = baseModel$call, knotsSelected = NULL, tempfits = c(CV = cv_initial, fitStat=fitStat))
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~ loop through 1D covar ~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  starttime = proc.time()
  timings<- vector(length=length(varlist))
  knots=NULL
  varkeepid<-NULL
  for (i in 2:(length(varlist)+1)){
    explanatory <- splineParams[[varID[(i-1)]]]$explanatory
    if(fam=='BinProp'){
      response<-cbind(data$successes, data$failures)
    }else{
      response<-data$response  
    }
    bd <- as.numeric(splineParams[[varID[(i-1)]]]$bd)   # i is the location of covar in varid +1 (2d has 1st entry in spline params)
    gap <- (salsa1dlist$gaps[(i-1)])
    term<- terms1D[[(i-1)]]
    interactionTerm<-NULL  #(salsa1dlist$interactionTerm[(i-1)])
    baseModel <- eval(parse(text=paste("update(baseModel, .~. -", term, ")", sep="")))
    
    if(removal==TRUE){
      cv_without<-getCV_CReSS(data, baseModel, splineParams=splineParams)
      fitStat_without<-get.measure(salsa1dlist$fitnessMeasure,'NA', baseModel, initDisp)$fitStat  
    }
    
    if(length(grep(varlist[(i-1)], baseModel$formula))>0){stop(paste('Multiple instances of covariate in model. Remove ',splineParams[[varID[(i-1)]]]$covar , ' before proceding', sep=''))}
    
    if(varlist[(i-1)]%in%varlist_cyclicSplines){spl<- "cc"}else{spl="bs"}
    
    if(spl == "cc"){salsa1dlist$minKnots_1d[(i-1)] <- 3; salsa1dlist$startKnots_1d[(i-1)]<-3}
    sttime<- proc.time()[3]
    output <- return.reg.spline.fit(response,explanatory,splineParams[[varID[(i-1)]]]$degree,salsa1dlist$minKnots_1d[(i-1)],salsa1dlist$maxKnots_1d[(i-1)],salsa1dlist$startKnots_1d[(i-1)], gap, winHalfWidth, salsa1dlist$fitnessMeasure, maxIterations=100, baseModel=baseModel, bd=bd, spl=spl, interactionTerm=interactionTerm)
    
    timings[(i-1)]<- proc.time()[3] - sttime
    
    thisFit <- output$output[2] 
    # update spline parameters if fit statistic for updated knots is better than for one knot at the mean.
    #if (thisFit < fitStat) {
    tempinitialknots<-splineParams[[varID[(i-1)]]]$knots
    splineParams[[varID[(i-1)]]]$knots= sort(output$aR)
    #}
    # update best model to have new knot locations and covariate back in model
    tempModel<- eval(parse(text=paste("update(baseModel, .~. +", term, ")", sep="")))
    # calculate a cv score here too
    if(removal==TRUE){
      cv_with<- getCV_CReSS(data=data, tempModel, splineParams)  
    }
    models<<-output$models
    knotSites<<-output$knotSites
    
    if(removal=='TRUE'){
    cat('Fitting Linear Model...')
    tempModel_lin<- eval(parse(text=paste("update(baseModel, . ~. +",varlist[(i-1)] , ")", sep="")))
    cv_linear<- getCV_CReSS(data=data, tempModel_lin, splineParams)
    
    cvid<-which(c(cv_initial, cv_with, cv_without, cv_linear)==min(na.omit(c(cv_initial, cv_with, cv_without, cv_linear))))[1]
    
    cat('Choosing smooth vs linear model...')
    if(cvid==1){
      # initial model is best - keep term but with original knots
      fitStat = fitStat
      splineParams[[varID[(i-1)]]]$knots<- tempinitialknots
      baseModel<-update(tempModel, .~.)
      cv_initial<-cv_initial
      varkeepid<-c(varkeepid, i)
      kept='YES - initial'
    }
    if(cvid==2){
      # model with covariate with new knots is best
      fitStat = thisFit
      splineParams[[varID[(i-1)]]]$knots= sort(output$aR)
      baseModel<-update(tempModel, .~.)
      cv_initial<-cv_with
      varkeepid<-c(varkeepid, i)
      kept='YES - new knots'
    }
      
    if(cvid==3){
        # model with parameter removed is best
        fitStat = fitStat_without
        splineParams[[varID[(i-1)]]]$knots<- 'NA'
        baseModel<-update(baseModel, .~.)
        cv_initial<-cv_without
        kept='NO'
    }
      
    if(cvid==4){
        # model with parameter linear is best
        splineParams[[varID[(i-1)]]]$knots<- 'NA'
        baseModel<-update(tempModel_lin, .~.)
        fitStat = get.measure(salsa1dlist$fitnessMeasure,'NA', baseModel, initDisp)$fitStat
        cv_initial<-cv_linear
        varkeepid<-c(varkeepid, i)
        kept='YES - linear'
    }
    }else{
#       cvid<-which(c(cv_initial, cv_with)==min(cv_initial, cv_with))
#       if(cvid==1){
#         # initial model is best - keep term but with original knots
#         fitStat = fitStat
#         splineParams[[varID[(i-1)]]]$knots<- tempinitialknots
#         baseModel<-update(tempModel, .~.)
#         cv_initial<-cv_initial
#         varkeepid<-c(varkeepid, i)
#         kept='YES - initial'
#       }
#      if(cvid==2){
        # model with covariate with new knots is best
        fitStat = thisFit
        splineParams[[varID[(i-1)]]]$knots= sort(output$aR)
        baseModel<-update(tempModel, .~.)
        cv_initial=cv_with=NULL
        varkeepid<-c(varkeepid, i)
        kept='YES - new knots'
#      }      
    }
    
    
    modelFits1D[[i]] <- list(term = term, kept=kept, basemodelformula = baseModel$call, knotsSelected = splineParams[[varID[(i-1)]]]$knots, baseModelFits = c(CV = cv_initial, fitStat = fitStat), modelfits = c(CV = cv_with, fitStat = thisFit))
    splineParams<<-splineParams
  }
  
  # turn the model data back into what came in
  
  outTerms <- list(length(varlist_cyclicSplines))
  counter<-1
  origfamily<-initialModel$family$family
  
  outModel<-eval(parse(text=paste("update(baseModel, .~., family=",substitute(origfamily), "(link=", substitute(link),"))",  sep='')))
  eval(parse(text=paste(substitute(datain),"<-data", sep="" )))
  #  for(i in 2:(length(varlist)+1)){
  #if(varlist[varID[(i-1)]-1]%in%varlist_cyclicSplines){
  # outTerms[[counter]]<-paste("as.matrix(data.frame(gam(response ~ s(", varlist[varID[(i-1)]-1], ", bs='cc', k=(length(splineParams[[", varID[(i-1)], "]]$knots) +2)), knots = list(",varlist[varID[(i-1)]-1], "=c(splineParams[[",varID[(i-1)], "]]$bd[1], splineParams[[", varID[(i-1)], "]]$knots, splineParams[[",varID[(i-1)], "]]$bd[2])), data=",substitute(datain),", fit=F)$X[,-1]))", sep="")  
  outModel<-eval(parse(text=paste("update(outModel, ~ ., data=", substitute(datain),")", sep="")))
  #    counter<-counter+1
  #   }
  class(outModel)<-c('gammrsea', class(outModel))
  outModel$varshortnames<-varlist[(varkeepid-1)]
  outModel$panels<-panelid
  
  attributes(outModel$formula)$.Environment<-.GlobalEnv
  
  if(salsa1dlist$fitnessMeasure[1]=='QAIC' | salsa1dlist$fitnessMeasure[1]=='QAICc' | salsa1dlist$fitnessMeasure[1]=='QBIC'){
    fitStatlist<-list(fitStat=fitStat, CV = cv_initial, chat=salsa1dlist$fitnessMeasure[2])  
  }else{
    fitStatlist<-list(fitStat=fitStat, CV = cv_initial)  
  }
  #save.image("Test.RData")
  return(list(bestModel=outModel, modelFits1D=modelFits1D, splineParams=splineParams, fitStat=fitStatlist, keptvarlist = varlist[(varkeepid-1)]))
}
