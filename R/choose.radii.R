#' Function to choose the radii for the CReSS local radial basis function
#' 
#' @export
#' 

choose.radii <- function(currentFit,indices,radiusIndices,radii,out.lm,dists,
                         aR,baseModel,fitnessMeasure,response,models, interactionTerm, data, initDisp, fit.opts, basis, printout) {
  #print("+++++++++++++++++++++++++++")
  #print("Fitting Radii")
  #print("+++++++++++++++++++++++++++")
  iterations = 0
  bestRadii=radiusIndices
  currentModel=out.lm
  improving = 1
  last = rep(0,length(radiusIndices))
  if(length(radii)> 1){
    #print("Fitting Radii")
    #cat('Start: ', radiusIndices, '\n')
    while (improving) {
      iterations = iterations + 1
      #print(c("iteration: ",iterations))
      improving = 0
      for (i in indices) {
        tempRadii = radiusIndices
        if ((tempRadii[i] > 1)&(last[i] <= 0)) {
          thisImprove = 1
          while (thisImprove) {
            thisImprove = 0
            tempRadii[i] = tempRadii[i] - 1
            output<- fit.thinPlate_2d(fitnessMeasure, dists,aR,radii,baseModel,tempRadii,models, currentFit, interactionTerm, data, initDisp, fit.opts, basis, printout)
            tempModel<-output$currentModel
            models<-output$models
#            tempFit <- get.measure_2d(fitnessMeasure,currentFit,tempModel, data, dists,aR,radii, tempRadii, initDisp)$fitStat
            tempFit <- output$fitStat
            #print(c(i,"Down",tempFit))
            if (tempFit < currentFit) {
              #print("UPDATING")
              bestRadii = tempRadii
              currentFit = tempFit
              currentModel = tempModel
              improving = 1
              if (tempRadii[i] > 1) thisImprove = 1
              last = rep(0,length(radiusIndices))
              last[i] = -1
            }
          }
        }
        tempRadii = radiusIndices
        if ((tempRadii[i]<length(radii))&(last[i] >= 0)) {
          thisImprove = 1
          while (thisImprove) {
            thisImprove = 0
            tempRadii[i] = tempRadii[i] + 1
            output<- fit.thinPlate_2d(fitnessMeasure, dists,aR,radii,baseModel,tempRadii,models, currentFit, interactionTerm, data, initDisp, fit.opts, basis, printout)
            tempModel<-output$currentModel
            models<-output$models
            tempFit <- output$fitStat
            #tempFit <- get.measure_2d(fitnessMeasure,currentFit,tempModel, data,  dists,aR, radii, tempRadii, initDisp)$fitStat
            #print(c(i,"Up",tempFit))
            if (tempFit < currentFit) {
              # print("UPDATING")
              bestRadii = tempRadii
              currentFit = tempFit
              currentModel = tempModel
              improving = 1
              if (tempRadii[i]<length(radii)) thisImprove = 1
              last = rep(0,length(radiusIndices))
              last[i] = 1
            }
          }
        }
      }
      #print(c("Current fit:",currentFit))
      radiusIndices = bestRadii
      if (isS4(currentModel)){
        currentModel@splineParams[[1]]$radiusIndices<-bestRadii
      } else {
        currentModel$splineParams[[1]]$radiusIndices<-bestRadii
      }
    }
  }# if length radii > 1 loop
  # print("+++++++++++++++++++++++++++")
  #print("Fitted Radii")
  #print("+++++++++++++++++++++++++++")
  list(BIC=currentFit,out.lm=currentModel,radiusIndices=radiusIndices,models=models)
}