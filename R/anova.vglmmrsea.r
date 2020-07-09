anova.vglmMRSea<-function(object, varshortnames=NULL, panelid=NULL, test='Wald'){
  
  x<-model.matrix.vlm(object)
  varlist<-attr(object@terms, 'variables')
  varseq<-attr(x, 'assign')
  varseq_adj <- c()
  for (vs in 1:length(varseq)){
    varseq_adj <- c(varseq_adj, rep(vs-1, length(varseq[vs][[1]])))
  }
  nvars<-max(0, varseq_adj)
  
  beta<-object@coefficients
  
  panelid <- if(is.null(panelid) & is.null(object@panels)) 1:nrow(object@data) else object@panels
  
  #warning message not to use F when correlation present
  if(length(unique(panelid))<nrow(object@data) && test=='F'){
    stop("If panels are provided (number of panels < number of observations), please use type = 'Wald' to make use of the robust standard errors.")
  }
  
  sumobj <- summaryvglm(object)
  
  if(max(table(panelid))==1){
    vbeta<-sumobj@cov.unscaled
  }else{
    vbeta<-sandcovS4(object, id=panelid)  
  }
  
  if(length(object@varshortnames)>0){
    varshortnames=object@varshortnames
    object<-summaryshortnamesS4(object, varshortnames)
  }
  
  ncoefs<-table(varseq_adj)
  
  if(test=='Wald'){
    X2Vec<-NULL
    
    for(j in 0:nvars){
      beta0<-rep(0, length(beta))
      beta0[varseq_adj==j]<-1
      #beta0<-rev(beta0)
      zeroidx<-beta0==1
      X2<-t(beta[zeroidx]%*% solve(vbeta[zeroidx, zeroidx, drop=FALSE]) %*% beta[zeroidx])
      X2Vec<-c(X2Vec, X2)
    }
    
    hasIntercept <- (length(grep("(Intercept)", names(beta))) != 0)
    
    if(hasIntercept){
      dfVec<-as.vector(table(varseq_adj[-1]))
    }
    
    resdf<-dfVec
    resdev<-X2Vec
    table<-data.frame(resdf, resdev, 1-pchisq(resdev, resdf))
    
  }
  
  if(test=='F'){
    which.nms <- function(name) which(varseq_adj == which(names == name))
    df.res <- df.residual(object)
    error.SS <- sum(residuals(object, "pearson")^2, na.rm = TRUE)
    fac <- attr(terms(object), "factors")
    names <- labels(terms(object))
    n.terms <- length(names)
    y <- object@y
    if (is.null(y)) 
      y <- model.response(model.frame(object), "numeric")
    wt <- object@prior.weights
    if (is.null(wt)) 
      wt <- rep(1, length(y))
    p <- df <- f <- SS <- rep(0, n.terms + 1)
    f[n.terms + 1] <- p[n.terms + 1] <- NA
    df[n.terms + 1] <- df.res
    SS[n.terms + 1] <- error.SS
    dispersion <- error.SS/df.res
    nyy <- ncol(y) - 1
    BasesFromMod <- get_all_bases(object)
    colnams <- colnames(BasesFromMod)
    varshortnam <- c(varshortnames, "LRF")
    for (term in 1:n.terms) {
      colmask <- c()
      varnam <- varshortnam[term]
      lenvarnam <- nchar(varnam)
      for (cn in 1:ncol(BasesFromMod)) {
        colshrt <- substr(colnams[cn],1,lenvarnam)
        if (colshrt==varnam) {
          colmask <- c(colmask, FALSE)
        } else {
          colmask <- c(colmask, TRUE)
        }
      }
      BasesNoVar <- BasesFromMod[,colmask]
      newdf <- as.data.frame(cbind(y, BasesNoVar))
      # deal with intercept only model
      if (sum(colmask) > 1) {
        text_mod <- paste0("mod.1<-vglm(cbind(", paste(colnames(y),collapse=','), ")~", 
                           paste(colnames(BasesNoVar)[2:ncol(BasesNoVar)],collapse='+'),
                           ",family='multinomial',data=newdf)")
      } else {
        text_mod <- paste0("mod.1<-vglm(cbind(", paste(colnames(y),collapse=','), 
                           ")~1,family='multinomial',data=newdf)")
      }
      eval(parse(text=text_mod))
      dev.1 <- deviance(mod.1)
      mod.2 <- object
      dev.2 <- deviance(mod.2)
      df[term] <- df.residual(mod.1) - df.residual(mod.2)
      if (df[term] == 0) {
        SS[term] <- f[term] <- p[term] <- NA
      } else {
        SS[term] <- dev.1 - dev.2
        f[term] <- SS[term]/(dispersion * df[term])
        p[term] <- pf(abs(f[term]), df[term], df.res, lower.tail = FALSE)
      }
    }
    result <- data.frame(SS, df, f, p)
    
  }
  
  tl<-attr(object@terms$terms, 'term.labels')
  
  if (length(varshortnames) > 0) {
    for (v in 1:length(tl)) {
      tl[grep(varshortnames[v], tl)]<-paste('s(',varshortnames[v],')', sep='')
    }
  }
  
  lr <-grep('LRF', tl)
  
  if(length(lr>0)){
    if(length(lr)>1){
      intid<-grep(':', tl)
      tl[lr[-which(lr==intid)]]<-'s(x.pos, y.pos)'
      splitint<-strsplit(tl[lr[2]], ':')
      id2<-which(c(length(grep('LRF', splitint[[1]][1])), length(grep('LRF', splitint[[1]][2])))==1)
      tl[lr[which(lr==intid)]]<-paste('s(x.pos, y.pos):', splitint[[1]][-id2], sep='')
    }else{
      tl[lr]<-'s(x.pos, y.pos)'
    }
  }
  
  if(test=='Wald'){
    if(max(table(panelid))==1){
      maxpanels <- paste(max(table(panelid)), ' (independence assumed)', sep='')
    }else{
      maxpanels<-max(table(panelid))
    }
    tl <- c("Intercept", tl)
    dimnames(table) <-list(c(tl), c("Df", "X2", "P(>|Chi|)"))
    title <- paste("Analysis of 'Wald statistic' Table", "\nModel: ",
                   object@family@vfamily[1], ", link: ", object@misc$link,
                   "\nResponse: ", as.character(varlist[-1])[1], "\nMarginal Testing\n", "Max Panel Size = ",maxpanels,"; Number of panels = ", length(unique(panelid)),"\n", sep = "")
    result<-structure(table, heading = title, class = c("anova", "data.frame"))
    
  }
  
  if(test=='F'){
    row.names(result) <- c(tl, "Residuals")
    names(result) <- c("SS", "Df", "F", "Pr(>F)")
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)", "Marginal Testing\n", 
                                 paste("Response:", as.character(varlist[-1])[1]), paste("Error estimate based on", "Pearson residuals", "\n"))
  }
  
  result
}


get_all_bases <- function(modin) {
  data <- modin@data
  splineParams <- modin@splineParams
  aR <- splineParams[[1]]$knotPos
  dists <- splineParams[[1]]$dist
  radii <- splineParams[[1]]$radii
  radiusIndices <- splineParams[[1]]$radiusIndices
  nbspl <- length(labels(terms(modin)))
  base_mod_cols <- modin@varshortnames
  if (length(base_mod_cols) > 0){
    for (nb in 1:length(base_mod_cols)){
      createcoltext <- paste0(base_mod_cols[nb],"<-data$", base_mod_cols[nb])
      eval(parse(text=createcoltext))
    }
  }
  BasesInMod <- eval(parse(text=labels(terms(modin))[1]))
  if (nbspl >= 2) {
    for (nb in 2:nbspl){
      newBas <- eval(parse(text=labels(terms(modin))[nb]))
      BasesInMod <- cbind(BasesInMod, newBas)
    }
  }
  # add intercept
  BasesInMod <- cbind(rep(1, nrow(BasesInMod)),BasesInMod)
  newnam <- create_temp_col_names(modin)
  colnames(BasesInMod) <- c("Intercept", newnam$columnz)
  
  return(BasesInMod)
}


create_temp_col_names <- function(mod_in){
  varz <- mod_in@varshortnames
  n_bases <- c()
  new_col_names <- c()
  if (length(varz)>0){
    for (vv in 1:length(varz)){
      n_bas <- length(mod_in@splineParams[[vv+1]]$knots) + 2
      n_bases <- c(n_bases, n_bas)
      new_col_names <- c(new_col_names, paste0(rep(varz[vv], n_bas), seq(n_bas)))
    }
  }
  n_smooth <- length(mod_in@splineParams[[1]]$knotPos)
  n_bases <- c(n_bases, n_smooth)
  new_col_names <- c(new_col_names, paste0(rep("LRF", n_smooth), seq(n_smooth)))
  return(list(columnz=new_col_names, nbases=n_bases))
}

