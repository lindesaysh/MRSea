update.gamMRSea<-function (object, formula., ..., evaluate = TRUE)
{
  #rm('splineParams', envir = globalenv())
  splineParams<<-object$splineParams
  
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
      newmodel<-eval(call, parent.frame())  
    #}else{
    #  newmodel<-eval(call, environment())
    #}
    rm('splineParams', envir = globalenv())
    newmodel$panels<-object$panels
    newmodel$varshortnames<-object$varshortnames
    newmodel$splineParams<-object$splineParams
    class(newmodel)<-class(object)
    newmodel
  }
  else{
    rm('splineParams', envir = globalenv())
    call
  }
}
