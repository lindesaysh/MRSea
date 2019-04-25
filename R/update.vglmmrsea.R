update.vglmMRSea<-function (object, formula., ..., evaluate = TRUE, panels=NULL)
{
  #rm('splineParams', envir = globalenv())
  ###splineParams<<-object@splineParams
  data <- object@data
  if (is.null(call <- getCall(object)))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.))
    call$formula <- update.formula(formula(object), formula.)
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate) {
    #if(is.null(splineParams)){
      #newmodel<-eval(call, parent.frame()) 
      newmodel<-eval(call) 
      class(newmodel)<-class(object)
    #}else{
    #  newmodel<-eval(call, environment())
    #}
    ###rm('splineParams', envir = globalenv())
    if(is.null(panels)){
      newmodel@panels<-object@panels  
    }else{
      newmodel@panels<-panels
    }
    newmodel@cvfolds<-object@cvfolds
    newmodel@varshortnames<-object@varshortnames
    newmodel@splineParams<-object@splineParams
    newmodel@data<-object@data
    newmodel
  }
  else{
    rm('splineParams', envir = globalenv())
    call
  }
}
