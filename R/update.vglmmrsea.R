update.vglmMRSea <- function (object, formula., ..., evaluate = TRUE, panels=NULL) 
{
  splineParams<<-object@splineParams
  data <- object@data
  if (is.null(call <- getCall(object))) 
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) {
    call$formula <- update_formula(formula(object), formula.)
  }
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate) {
    newmodel<-eval(call) 
    rm('splineParams', envir = globalenv())
    if(is.null(panels)){
      newmodel@panels<-object@panels  
    }else{
      newmodel@panels<-panels
    }
    newmodel@cvfolds<-object@cvfolds
    newmodel@varshortnames<-object@varshortnames
    newmodel@splineParams<-object@splineParams
    newmodel@data <- data
    class(newmodel)<-class(object)
    newmodel
  }
  else{
    rm('splineParams', envir = globalenv())
    call
  }
}


# 
# update.vglmMRSea<-function (object, formula., ..., evaluate = TRUE, panels=NULL)
# {
#   #rm('splineParams', envir = globalenv())
#   ###splineParams<<-object@splineParams
#   data <- object@data
#   #constraints <- object@constraints
#   constraints <- NULL
#   #etastart <- object@etastart
#   #mustart <- object@mustart
#   #coefstart <- object@coefstart
#   
#   if (is.null(call <- getCall(object)))
#     stop("need an object with call component")
#   extras <- match.call(expand.dots = FALSE)$...
#   if (!missing(formula.))
#     call$formula <- update.formula(formula(object), formula.)
#   extras$constraints <- constraints
#   #extras$etastart <- etastart
#   #extras$mustart <- mustart
#   #extras$coefstart <- coefstart
#   if (length(extras)) {
#     existing <- !is.na(match(names(extras), names(call)))
#     for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
#     if (any(!existing)) {
#       call <- c(as.list(call), extras[!existing])
#       call <- as.call(call)
#     }
#   }
#   if (evaluate) {
#     #if(is.null(splineParams)){
#       #newmodel<-eval(call, parent.frame())
#       newmodel<-eval(call)
#       print("class")
#       print(class(object))
#       class(newmodel)<-class(object)
#     #}else{
#     #  newmodel<-eval(call, environment())
#     #}
#     ###rm('splineParams', envir = globalenv())
#     if(is.null(panels)){
#       newmodel@panels<-object@panels  
#     }else{
#       newmodel@panels<-panels
#     }
#     newmodel@cvfolds<-object@cvfolds
#     newmodel@varshortnames<-object@varshortnames
#     newmodel@splineParams<-object@splineParams
#     newmodel@data<-object@data
#     newmodel@interactionterm <- object@interactionterm
#     newmodel
#   }
#   else{
#     rm('splineParams', envir = globalenv())
#     call
#   }
# }
