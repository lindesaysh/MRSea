
anova.gamMRSea<-function(object, varshortnames=NULL, panelid=NULL){

  if(!is.null(object$varshortnames)){
    varshortnames=object$varshortnames
  }


  x<-model.matrix(object)
  varlist<-attr(object$terms, 'variables')
  varseq<-attr(x, 'assign')
  nvars<-max(0, varseq)

  beta<-object$coefficients

  panelid <- if(is.null(panelid) & is.null(object$panels)) 1:nrow(object$data) else object$panels
  vbeta<-sandcov(object, id=panelid)
  ncoefs<-table(varseq)
  X2Vec<-NULL

  for(j in 1:nvars){
    beta0<-rep(0, length(beta))
    beta0[varseq==j]<-1
    #beta0<-rev(beta0)
    zeroidx<-beta0==1
    X2<-t(beta[zeroidx]%*% solve(vbeta[zeroidx, zeroidx, drop=FALSE]) %*% beta[zeroidx])
    X2Vec<-c(X2Vec, X2)
  }

  hasIntercept <- (length(grep("(Intercept)", names(beta))) != 0)

  if(hasIntercept){
    dfVec<-as.vector(table(varseq[-1]))
  }

  resdf<-dfVec
  resdev<-X2Vec
  table<-data.frame(resdf, resdev, 1-pchisq(resdev, resdf))
  tl<-attr(object$terms, 'term.labels')

  if (length(varshortnames) > 0) {
    for (v in 1:length(tl)) {
      tl[grep(varshortnames[v], tl)]<-paste('s(',varshortnames[v],')', sep='')
    }
  }

  dimnames(table)<-list(c(tl), c("Df", "X2", "P(>|Chi|)"))
  title <- paste("Analysis of 'Wald statistic' Table", "\nModel: ",
                 object$family$family, ", link: ", object$family$link,
                 "\nResponse: ", as.character(varlist[-1])[1], "\nMarginal Testing\n", "Max Panel Size = ",max(table(panelid)),"; Number of panels = ", max(panelid),"\n", sep = "")
  structure(table, heading = title, class = c("anova", "data.frame"))
}
