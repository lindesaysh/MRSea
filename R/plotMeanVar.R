#-----------------------------------------------------------------------------
#' Functions to create a Mean-Variance plot for checking the distribution assumptions of the mean and the variance.  Distributions available are Gaussian, Poisson, QuasiPoisson, Gamma and Tweedie. 
#' 
#' @param model Fitted model object (glm or gam)
#' @param cut.prob.by Numerical input to state the increment for the sequence of cut probabilities. 
#' @param save (\code{default=FALSE}). Logical stating whether plot should be saved into working directory. See \code{label} to change directory.
#' @param label Character string indicating an label to be added to the plot when using \code{save = TRUE}. Can also include a pathway to a directory of choice.
#' @param print Logical stating whether or not to print the plot. If FALSE then the plot object is returned. 
#' 
#' @return
#' A plot showing the observed mean and variance (cutting the fitted values into bins and finding the mean fitted value and the variance for each bin) and the assumed relationship under various distributions depending on the model fitted (lines on the plot). 
#' 
#' @examples 
#' # load data
#' data(ns.data.re)
#' 
#' model<-glm(birds ~ observationhour + as.factor(floodebb) + as.factor(impact), 
#'            family='quasipoisson', data=ns.data.re)
#' 
#' plotMeanVar(model)
#'
#' @export
#'

plotMeanVar<-function(model, cut.bins = 20, save=FALSE, label = NULL, print=TRUE){
  
  cutpts<-unique(quantile(fitted(model), prob=seq(0,1, length=cut.bins)))
  mycuts<-cut(fitted(model), breaks= cutpts)
  meanfits<-tapply(fitted(model), mycuts, mean)
  varresid<-tapply(residuals(model, type='response'), mycuts, var)

  if(length(which(is.na(meanfits)))>0){
    stop("Error in number of bins, try reducing using cut.bins=.  The default is 20.")
  }
     
  p <- ggplot() +
    geom_point(aes(meanfits, varresid), size=2) +
    xlab('Fitted Values (Mean)') +  ylab('V(residuals)') +
    theme_bw() +
    theme(panel.grid.major=element_blank(), 
          axis.text.x=element_text(size=15), 
          axis.text.y=element_text(size=15), 
          axis.title.x=element_text(size=15), 
          axis.title.y=element_text(size=15), 
          plot.title=element_text(size=15),
          legend.position = "top")
  
  
  if(model$family[[1]] == "poisson"){
    pdat <- tibble(mu = seq(min(meanfits), max(meanfits), length=100)) %>% 
      mutate(
        V.p = mu
      )  
    
  p <- p + 
     geom_line(data=pdat, aes(x=mu, y=V.p), linewidth=1, linetype=2, colour="darkgrey") #+
      # scale_colour_manual("Distritubion",
      #                     values=c("darkgrey"),
      #                     breaks = c("V.p"),
      #                     labels=c("Poisson"))
  }
  
  
  if(model$family[[1]] == "quasipoisson"){
    pdat <- tibble(mu = seq(min(meanfits), max(meanfits), length=100)) %>% 
      mutate(
        V.p = mu, 
        V.qp = mu * summary(model)$dispersion
      )  
    
    pdat <- tidyr::pivot_longer(pdat, 
                                cols = V.p:V.qp, 
                                names_to = "Distribution", 
                                values_to = "EstVariance")
    
    p <- p + 
      geom_line(data=pdat, aes(x=mu, y=EstVariance, group=Distribution, colour=Distribution, linetype=Distribution), linewidth=1) +
      scale_colour_manual("Distribution", 
                          values=c("darkgrey", "firebrick3"),
                          labels=c("Poisson", "QuasiPoisson")) +
      scale_linetype_manual("Distribution", 
                            values = c(2,1),
                            labels=c("Poisson", "QuasiPoisson"))
    }
  
  if(model$family[[1]] == "Tweedie"){
    pdat <- tibble(mu = seq(min(meanfits), max(meanfits), length=100)) %>% 
      mutate(
      V.mu = model$family$variance(mu),
      V.p = mu, 
      V.qp = mu * summary(model)$dispersion,
      V.tw = V.mu * summary(model)$dispersion
    )  
    
    pdat <- tidyr::pivot_longer(pdat, 
                                cols = V.p:V.tw, 
                                names_to = "Distribution", 
                                values_to = "EstVariance")
    
    p <- p + 
      geom_line(data=pdat, aes(x=mu, y=EstVariance, group=Distribution, colour=Distribution, linetype=Distribution), linewidth=1) +
      scale_colour_manual("Distribution",
                          values=c("darkgrey", "firebrick3", "deepskyblue2"),
                          labels=c("Poisson", "QuasiPoisson", "Tweedie"))+
      scale_linetype_manual("Distribution", 
                            values = c(2,1,3),
                            labels=c("Poisson", "QuasiPoisson", "Tweedie"))
  }
  
  if(model$family[[1]] == "Gamma"){
    pdat <- tibble(mu = seq(min(meanfits), max(meanfits), length=100)) %>% 
      mutate(
        V.p = mu, 
        V.g = model$family$variance(mu) * summary(model)$dispersion
      )  
    
    pdat <- tidyr::pivot_longer(pdat, 
                                cols = V.p:V.g, 
                                names_to = "Distribution", 
                                values_to = "EstVariance")
    
    p <- p + 
      geom_line(data=pdat, aes(x=mu, y=EstVariance, group=Distribution, colour=Distribution, linetype=Distribution), linewidth=1) +
      scale_colour_manual("Distribution",
                          values=c("darkgrey", "firebrick3"),
                          labels=c("Poisson", "Gamma")) +
      scale_linetype_manual("Distribution", 
                            values = c(2,1),
                            labels=c("Poisson", "Gamma"))
  }
  
  if(model$family[[1]] == "gaussian"){
    pdat <- tibble(mu = seq(min(meanfits), max(meanfits), length=100)) %>% 
      mutate(
        V.p = mu, 
        V.g = summary(model)$dispersion
      )  
    
    pdat <- tidyr::pivot_longer(pdat, 
                                cols = V.p:V.g, 
                                names_to = "Distribution", 
                                values_to = "EstVariance")
    
    p <- p + 
      geom_line(data=pdat, aes(x=mu, y=EstVariance, group=Distribution, colour=Distribution, linetype=Distribution), linewidth=1) +
      scale_colour_manual("Distribution",
                          values=c("darkgrey", "firebrick3"),
                          labels=c("Poisson", "Gaussian")) +
      scale_linetype_manual("Distribution", 
                            values = c(2,1),
                            labels=c("Poisson", "Gaussian"))
  }
  
  # if(model$family[[1]] == "binomial"){
  #   pdat <- tibble(mu = seq(min(meanfits), max(meanfits), length=100)) %>% 
  #     mutate(
  #       V.b = model$family$variance(mu),
  #       V.b2 = (mu) * (1-mu)
  #     )  
  #   
  #   p <- p + 
  #     geom_line(data=pdat, aes(x=mu, y=V.b2), linewidth=1) +
  #     scale_colour_manual(values=c("firebrick3"),
  #                         breaks=c("V.b2"),
  #                         labels=c("Binomial"))
  # }
  
  if(save==T){ggsave(paste0(label, "MeanVarplot.png"), a, height=6, width=8)
  }
  if(print==FALSE){
    return(p)
  }else{
    print(p)
  }
}
