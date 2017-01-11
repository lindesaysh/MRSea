a=
anova.gamMRSea<-function(object, varshortnames=NULL, panelid=NULL, test='Wald'){

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
  
  if(test=='Wald'){
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
    
  }
  
  if(test=='F'){
    which.nms <- function(name) which(varseq == which(names == 
                                                      name))
    df.res <- df.residual(object)
    error.SS <- sum(residuals(object, "pearson")^2, na.rm = TRUE)
    
    fac <- attr(terms(object), "factors")
    names <- labels(terms(object))
    n.terms <- length(names)
    y <- object$y
    if (is.null(y)) 
      y <- model.response(model.frame(object), "numeric")
    wt <- object$prior.weights
    if (is.null(wt)) 
      wt <- rep(1, length(y))
    p <- df <- f <- SS <- rep(0, n.terms + 1)
    f[n.terms + 1] <- p[n.terms + 1] <- NA
    df[n.terms + 1] <- df.res
    SS[n.terms + 1] <- error.SS
    dispersion <- error.SS/df.res
  
    for (term in 1:n.terms) {
      rels <- names[relatives(names[term], names, fac)]
      exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), 
                                           which.nms)))
      mod.1 <- glm.fit(x[, -exclude.1, drop = FALSE], y, wt, 
                       offset = object$offset, family = object$family, control = object$control)
      dev.1 <- deviance(mod.1)
      mod.2 <- if (length(rels) == 0) 
        object
      else {
        exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
        glm.fit(x[, -exclude.2, drop = FALSE], y, wt, offset = object$offset, 
                family = object$family, control = object$control)
      }
      dev.2 <- deviance(mod.2)
      df[term] <- df.residual(mod.1) - df.residual(mod.2)
      if (df[term] == 0) 
        SS[term] <- f[term] <- p[term] <- NA
      else {
        SS[term] <- dev.1 - dev.2
        f[term] <- SS[term]/(dispersion * df[term])
        p[term] <- pf(f[term], df[term], df.res, lower.tail = FALSE)
      }
    }
    result <- data.frame(SS, df, f, p)
    }
  
  
  tl<-attr(object$terms, 'term.labels')

  if (length(varshortnames) > 0) {
    for (v in 1:length(tl)) {
      tl[grep(varshortnames[v], tl)]<-paste('s(',varshortnames[v],')', sep='')
    }
  }

  lr<-grep('LocalRadialFunction', tl)
  if(length(lr>0)){
    if(length(lr)>1){
      intid<-grep(':', tl)
      tl[lr[-which(lr==intid)]]<-'s(x.pos, y.pos)'
      
      splitint<-strsplit(tl[lr[2]], ':')
      id2<-which(c(length(grep('LocalRadialFunction', splitint[[1]][1])), length(grep('LocalRadialFunction', splitint[[1]][2])))==1)
      tl[lr[which(lr==intid)]]<-paste('s(x.pos, y.pos):', splitint[[1]][-id2], sep='')
    }else{
      tl[lr]<-'s(x.pos, y.pos)'
    }
  }
  
  if(test=='Wald'){
    dimnames(table)<-list(c(tl), c("Df", "X2", "P(>|Chi|)"))
    title <- paste("Analysis of 'Wald statistic' Table", "\nModel: ",
                   object$family$family, ", link: ", object$family$link,
                   "\nResponse: ", as.character(varlist[-1])[1], "\nMarginal Testing\n", "Max Panel Size = ",max(table(panelid)),"; Number of panels = ", max(panelid),"\n", sep = "")
    result<-structure(table, heading = title, class = c("anova", "data.frame"))
    
  }
  
  if(test=='F'){
  row.names(result) <- c(tl, "Residuals")
  names(result) <- c("SS", "Df", "F", "Pr(>F)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n", 
                               paste("Response:", as.character(varlist[-1])[1]), paste("Error estimate based on", "Pearson residuals", "\n"))
  }
  
  result
}




relatives<-function (term, names, factors) 
{
  is.relative <- function(term1, term2) {
    all(!(factors[, term1] & (!factors[, term2])))
  }
  if (length(names) == 1) 
    return(NULL)
  which.term <- which(term == names)
  (1:length(names))[-which.term][sapply(names[-which.term], 
                                        function(term2) is.relative(term, term2))]
}