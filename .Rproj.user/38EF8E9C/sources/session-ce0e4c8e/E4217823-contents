#' Shortening names in summary object
#' 
#' @param object a model object
#' @param varshortnames a list of short variable names
#' 
#' @export
#' 
summaryshortnames<-function(object, varshortnames){
  

#bob <- attr(object$coefficients, "names")

#factorlist<-names(which(attr(terms(object), 'dataClasses')=='factor'))

if(length(varshortnames)>0){
  bob2<-attr(terms(object), 'term.labels')
 
  varshortnames.temp<-varshortnames#[varshortnames%in%bob2]
  
  if(length(varshortnames.temp)==0){
    v<-NULL
    for(i in 1:length(varshortnames)){
      v<-c(v, ifelse(length(grep(varshortnames[i], bob2))>0, 1, 0))
    }
    varshortnames<-varshortnames[which(v==1)]
  }else{
    varshortnames<-varshortnames.temp
  }
}

#varshortnames <- fit.int$varshortnames
b <- attributes(object$terms)$term.labels
namingtable <- NULL

for(i in 1:length(b)){
  newname <- b[i]
  # does it have colon
  intcheck <- strsplit(b[i], ":")[[1]]
  if(length(intcheck)>1){
    newname <- vector(length=length(intcheck))
    for(d in 1:length(intcheck)){
      if(length(grep("splineParams", intcheck[d]))>0){
        nameid <- which(stringr::str_detect(intcheck[d], varshortnames, negate = FALSE)==TRUE)
        newname[d] <- paste0("s(", varshortnames[nameid],")")
      }else{
        newname[d] <- intcheck[d]
      }
    } # d loop 
    nn <- cbind(intcheck, newname)
  }else{
    # is it in varshortnames
    nameid <- which(stringr::str_detect(b[i], varshortnames, negate = FALSE)==TRUE)
    if(length(nameid)>0){
      newname <- paste0("s(", varshortnames[nameid], ")")
    }else{
      if(length(grep("LRF", b[i]))>0){
        newname <- "s(x,y)"
      }else{
      newname <- b[i]
    }}
    
    nn <- c(b[i], newname)
  }
  
  namingtable <- rbind(namingtable, nn)
  
}

namingtable <- namingtable[!duplicated(namingtable),]

coefnames <- attr(object$coefficients, "names")

if(is.null(nrow(namingtable))){
  coefnames <- stringr::str_replace(coefnames, fixed(namingtable[1]), namingtable[2])
}else{
  for(nam in 1:nrow(namingtable)){
    coefnames <- stringr::str_replace(coefnames, fixed(namingtable[nam,1]), namingtable[nam,2])
  }
}


attr(object$coefficients, "names") <- coefnames

return(object)

}


