#' Shortening names in summary object
#' 
#' @param object a model object
#' @param varshortnames a list of short variable names
#' 
#' @export
#' 
summaryshortnames<-function(object, varshortnames){
  

bob <- attr(object$coefficients, "names")

factorlist<-names(which(attr(terms(object), 'dataClasses')=='factor'))

if(is.null(varshortnames)==FALSE | length(varshortnames)==0){
  bob2<-attr(terms(object), 'term.labels')
  v<-NULL
  for(i in 1:length(varshortnames)){
    v<-c(v, ifelse(length(grep(varshortnames[i], bob2))>0, 1, 0))
  }
  varshortnames<-varshortnames[which(v==1)]
}

# interaction terms id
id_int <- grep(':', bob)

# main effect id
if(length(id_int)==0){
  idmaineffect<-1:length(bob)
}else{
  idmaineffect <- (1:length(bob))[-id_int]  
}


# remove factors from shortnames
finvid<-NULL
if(length(factorlist)>0){
  for(i in 1:length(factorlist)){
    finvid<-c(finvid, grep(factorlist[i], varshortnames))
  }
  if(length(finvid)>0 & length(varshortnames)>0){
    varshortnames<-varshortnames[-finvid]
  }
}

# remove linear terms from shortnames

if(length(varshortnames)>0){
  for (i in 1:length(varshortnames)) {
    idvar <- grep(varshortnames[i], bob)
    idvarmain <- na.omit(match(idvar, idmaineffect))
    if (length(idvarmain) > 1) {
      for (j in 1:length(idvarmain)) {
        bob[idvarmain][j] <- paste("s(", varshortnames[i], ")", j, sep = "")
      }
    }else {bob[idvarmain] <- paste(varshortnames[i], sep = "")}
  }
}

# general interaction terms

if(length(id_int)>0){
  if(length(varshortnames)>0){
    for(v in 1:length(varshortnames)){
      idvar<-grep(varshortnames[v], bob[id_int])
      idvar_main<-grep(varshortnames[v], bob[idmaineffect])
      if (length(idvar) > 1) {
        counter<-1
        for (j in 1:(length(idvar)/length(idvar_main))) {
          for (i in 1:length(idvar_main)){
            bob[id_int[idvar][counter]] <- paste(strsplit(bob[id_int[idvar][counter]], ':')[[1]][1], ":s(", varshortnames[v], ")", i, sep = "")
            counter<-counter+1  
          }
        }
      }
      if(length(idvar)==1){bob[id_int[idvar]] <- paste(varshortnames[i], sep = "")}
    }  
  }
}


# local radial terms
localid <- grep("LRF", bob)
localint<-id_int[na.omit(match(localid, id_int))]
if (length(localid > 1)) {
  #intid <- grep(":", bob)
  smoothid <- localid[which(is.na(match(localid, localint)))]
  for (k in 1:length(smoothid)) {
    bob[smoothid][k] <- paste("s(x.pos, y.pos)b", k,
                              sep = "")
  }
  
  # local radial interactions
  
  counter <- 1
  for (k in 1:(length(localint)/length(smoothid))) {
    for (i in 1:length(smoothid)) {
      if(i==1){
        nonlocalid<-c(2,1)[grep(pattern = 'LRF', strsplit(bob[localint[counter]], ":")[[1]])]
      }
    textin <- paste(strsplit(bob[localint[counter]], ':')[[1]][nonlocalid], ":s(x.pos, y.pos)b", i, sep = "")
    #print(textin)
       bob[localint[counter]] <- textin
       counter <- counter + 1
    }
  }
}
attr(object$coefficients, "names") <- bob

return(object)

}