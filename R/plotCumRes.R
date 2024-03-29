#-----------------------------------------------------------------------------
#' Calculate cumulative residuals and plot.  
#' 
#' The output is plots of cumulative residuals.  
#' 
#' @param model Fitted model object (glm or gam)
#' @param varlist Vector of covariate names (continous covariates only)
#' @param label Label printed at the end of the plot name to identify it if \code{save=TRUE}.
#' @param save (\code{default=FALSE}). Logical stating whether plot should be saved into working directory.
#'
#' @return
#' Cumulative residual plots are returned for residuals ordered by each covariate in \code{varlist}, predicted value and index of observations (temporally).
#' The blue dots are the residuals
#' The black line is the line of cumulative residual.
#' On the covariate plots (those in \code{varlist}) the grey line indicates what we would expect from a well fitted covariate. i.e. one that is fitted with excessive knots.
#' 
#' Note: if the covariate is discrete in nature (like the example below), there will be a lot of overplotting of residuals.
#'
#' @examples 
#' # load data
#' data(ns.data.re)
#' 
#' model<-gamMRSea(birds ~ observationhour + as.factor(floodebb) + as.factor(impact), 
#'            family='quasipoisson', data=ns.data.re)
#' 
#' plotCumRes(model, varlist=c('observationhour'))
#' 
#' @export
#' 
plotCumRes<- function(model, varlist=NULL, label='', save=FALSE, variableonly = FALSE){
  
  require(splines)
  
  print("Calculating cumulative residuals")
  
  if(is.null(varlist)){
    namesOfx<- c(c("Predicted", "Index"))  
  }else{
    namesOfx<- c(varlist, c("Predicted", "Index"))
  }
  
  if(class(model)[1]=='geeglm' | class(model)[1]=='glm'){
    dat<- data.frame(model$data)
  }
  
  if(class(model)[1]=='gam'){
    dat<-data.frame(model$model)
  }
  
  if(class(model)[1]=='gamMRSea'){
    dat<-data.frame(model$data)
  }
  
  if(!is.null(varlist)){
    # find which data columns refer to the used variables
    coefpos<-c()
    for(z in 1:(length(namesOfx)-2)){
      coefpos<- c(coefpos, grep(namesOfx[z], names(dat)))
    }
  }
  
  # make matrix (n x number covar + fitted + orderdata)
  #xmatrix<- cbind(dat[,coefpos], fitted(model), ord=1:length(fitted(model)))
  
  dat <- dat %>% 
          mutate(Predicted = fitted(model),
                 Index = 1:n())
  
  #attributes(model$formula)$.Environment<-environment()

  if(variableonly == FALSE){
    plotvar <- namesOfx
  }else{
    plotvar <- varlist
  }
    
    
  for(z in 1:length(plotvar)){
    
    type="response"
    yl <- "Response Residuals"
    if(plotvar[z] == "Predicted"){
      type = "response"
      yl = "Response Residuals"
    }
    
    
    
    
    plotdat <- dat %>% 
      mutate(m.resids = residuals(model, type=type))
    
    plotdat <- eval(parse(text = paste0("arrange(plotdat, ", plotvar[z], ")")))
    
    plotdat <- plotdat %>% 
      mutate(m.cumsum = cumsum(m.resids))
    
    if(z < (length(namesOfx)-1)){
      covardat<- dat[model$y>0, coefpos[z]]
      newknots<-as.vector(quantile(covardat, seq(0.05, 0.95, length=7)))
      newknots<-unique(newknots)
      #dists<-d2k
      
      # removal term
      #term<-attributes(terms(model))$term.labels[grep(namesOfx[i], attributes(terms(model))$term.labels)]
      term<-labels(terms(model))[grep(namesOfx[z], labels(terms(model)))]
      #newterm<- paste('bs(', namesOfx[c], ', knots= c(newknots)', sep='')
      newterm<- paste('bs(', namesOfx[z], ', knots= c(', paste(newknots, sep=" ", collapse=','), '))', sep='')
      eval(parse(text=paste('covarModelUpdate<-update(model, .~. -',term, ' +', newterm,', data=plotdat)', sep='')))
      #xmatrix_test<- cbind(dat[,coefpos], fitted(covarModelUpdate), ord=1:length(fitted(covarModelUpdate)))
      plotdat <- plotdat %>% 
        mutate(new.fits = fitted(covarModelUpdate),
               new.resids = residuals(covarModelUpdate, type=type),
               new.cumsum = cumsum(new.resids))
    }
    
    
    maxy <- max(plotdat$m.resids, plotdat$m.cumsum)
    miny <- min(plotdat$m.resids, plotdat$m.cumsum)
    
    
   p <- ggplot(plotdat) + 
      geom_point(aes(x=pull(plotdat, plotvar[z]), y=m.resids), colour = "turquoise4") +
      xlab(plotvar[z]) + ylab(yl) +
      geom_line(aes(x=pull(plotdat, plotvar[z]), y= m.cumsum)) +
      geom_hline(yintercept = 0) +
     theme_bw()
    
   if(z < (length(namesOfx)-1)){  
    p <- p + geom_line(aes(x=pull(plotdat, plotvar[z]), y=new.cumsum), colour = "grey") +
      geom_point(aes(x=pull(plotdat, plotvar[z]), y=new.resids), colour = "grey", alpha=1/5)
   }
    #maxy<- max(residuals(model),cumsum(residuals(model, type=type)[order(xmatrix[,z])]))
    #miny<- min(residuals(model),cumsum(residuals(model, type=type)[order(xmatrix[,z])]))
    
    # make plot
    if(save==T){png(paste("CumRes_", namesOfx[z],label, ".png", sep=''), height=600, width=700)}
    else{devAskNewPage(ask=TRUE)}
   print(p)
    #plot(xmatrix[,z][order(xmatrix[,z])],residuals(model, type=type)[order(xmatrix[,z])], ylim=c(miny,maxy), xlab=namesOfx[z], ylab=yl, pch='.',cex=0.1, main="Cumulative Residuals", cex.lab=1.3, cex.axis=1.3)
    
    # if(z < (length(namesOfx)-1)){
    #   lines(xmatrix_test[,z][order(xmatrix_test[,z])],cumsum(residuals(covarModelUpdate, type=type)[order(xmatrix_test[,z])]), lwd=2, col='grey')  
    # }
    #points(xmatrix[,z][order(xmatrix[,z])],residuals(model, type=type)[order(xmatrix[,z])], pch=20, col='turquoise4', cex=0.5)
    #lines(xmatrix[,z][order(xmatrix[,z])],cumsum(residuals(model, type=type)[order(xmatrix[,z])]), lwd=2)
    #abline(h=0)
    if(save==T){dev.off()}else{devAskNewPage(ask=FALSE)}
  }
}
