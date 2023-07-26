#' Function to generate the empirical distribution of the runs test statistic, given some data and a model.
#' 
#' @param n.sim The number of simulated sets of data.
#' @param simData A matrix, where each column is a set of data simulated under independence, with rows the same length as the data used for the model.
#' @param model a glm or gamMRSea model object
#' @param data data set used to fit the model.
#' @param plot logical flag. If TRUE, a plot is made showing the 5\% critical values for the empirical distribution vs the N(0,1) distribution. Default is 'FALSE'
#' @param returnDist logical flag.  Do you want the distribution of test statistics returned ('TRUE') or just the 5\% critical values ('FALSE')
#' @param dots (Default: \code{TRUE}. Logical stating whether to show a printout of block numebers to assess progress. 'TRUE' will print dots into the workspace. 
#' 
#'
#' @examples 
#' data(ns.data.re)
#' 
#' model<-gamMRSea(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),  
#'               family='poisson', data=ns.data.re)
#' 
#' simData<-generateNoise(n=500, response=fitted(model), family='poisson')
#' empdist<-getEmpDistribution(500, simData, model, data=ns.data.re, plot=FALSE, 
#'                 returnDist=TRUE,dots=FALSE)
#' 
#' @author LAS Scott-Hayward, University of St Andrews
#' 
#' @export
#' 
getEmpDistribution<-function(n.sim, simData, model, data, plot=FALSE, returnDist=TRUE,dots=TRUE){
  runsstatH0=betas=vector(length=n.sim)
  model$panels<-1:nrow(data) # make sure that independent panels specified for null distribution
for(i in 1:n.sim){
  if(dots==TRUE){if((i/500)%%1 == 0){cat(i, '\n')}else{cat('.')}}

  splineParams = model$splineParams
  data$response=simData[,i]
  sim_glm<- update(model, response ~ ., data=data, splineParams=splineParams)
  runs<-runsTest(residuals(sim_glm, type='pearson'))
  runsstatH0[i]<-runs$statistic
  # if(imp.beta.coverage){
  #   betas[i]<-coefficients(model)[length(coefficients(model))]
  # }
}

  critvals<-quantile(runsstatH0, probs = c(0.025, 0.975))

  if(plot){
    hist(runsstatH0, main='Empirical Distribution: Runs Test Statistic', xlim=range(c(qnorm(0.025), qnorm(0.975), runsstatH0)), xlab='Test Statistic')
    abline(v=critvals, col='red', lwd=2)
    abline(v=c(qnorm(0.025), qnorm(0.975)), col='blue', lwd=2)
  }



  # if(returnDist){
  #   output<-list=c(empirical.distribution.null=runsstatH0, beta.coverage.null = )
  # }else{
  #   output<-list=c(empirical.critical=as.vector(critvals), beta.coverage.null = )
  #   }
  if(returnDist){
    output<-runsstatH0
  }else{
    output<-as.vector(critvals)
  }
  return(output)
}
