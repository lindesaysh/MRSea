# devtools::install_github("lindesaysh/MRSea")
devtools::load_all("C:/Users/kryzi/OneDrive - University of St Andrews/PhD/Code/MRsea")
# library(MRSea)
require(tidyverse)
require(patchwork)
require(knitr)
require(fields)
require(MuMIn)

count.data <- read.csv("1.count.data.csv")
data <- count.data
data$response <- round(data$NHAT)
data$blockid <- paste(data$transect.id, data$season, data$impact,sep = "")
attach(data)

data("nysted.predictdata")  # contains predict.data
# This is a spatial grid for making predictions.  All covariates in 
# final model must be in this data frame and the naming must be the 
# same as for the data
predictData <- nysted.predictdata

# Not including impact or season in the model so included data for only one impact and season, impact = 1 and season = 1 has highest count
predictData <- predictData[predictData$impact==1&predictData$season==1,]

initialModel <- glm(response ~ 1 + offset(log(area)), 
                    family = "quasipoisson", data = data)

salsa1dlist <- list(
  fitnessMeasure = "cv.gamMRSea", 
  minKnots_1d = 2,
  maxKnots_1d = 5, 
  startKnots_1d = 1, 
  degree = 2, 
  maxIterations = 10, 
  gaps = c(0), 
  cv.opts=list(cv.gamMRSea.seed=1, K=10)
)

set.seed(1234)
salsa1dOutput <- runSALSA1D(initialModel, salsa1dlist, c("depth"), predictionData=predictData, datain=data,removal=TRUE)

# How many knots were chosen for depth?
depth_knots <- salsa1dOutput$splineParams[[2]]$knots
# ~~~~~~~~~~~~~~~~~~~~~~~

cat("Number of knots selected for depth:", length(depth_knots), "this should be 1. \n\n")
cat("Location of knot selected for depth:", round(depth_knots, 4), "this should be -12.212. \n\n")

knotgrid <- read.csv("4.knotgrid.csv")

cat("Sum of knotgrid:", sum(knotgrid), "this should be 2018630. \n\n")


# make distance matrices for datatoknots and knottoknots
distMats <- makeDists(cbind(data$x.pos, data$y.pos), knotgrid)

# ~~~~~~~~~~~~~~~~~~~~~~~

cat("Sum of datadist:", sum(distMats$dataDist[1,]), "this should be 8224.849. \n\n")
cat("Sum of knotdist:", sum(distMats$knotDist), "this should be 1775357. \n\n")

# make parameter set for running salsa2d
salsa2dlist<-list(
  fitnessMeasure = 'cv.gamMRSea', 
  knotgrid = knotgrid, 
  startKnots=5, 
  minKnots=4, 
  maxKnots=12, 
  gap=0, 
  cv.opts=list(cv.gamMRSea.seed=1, K=10)
)

start_time <- Sys.time()
salsa2dOutput<-runSALSA2D(salsa1dOutput$bestModel, salsa2dlist, d2k=distMats$dataDist, k2k=distMats$knotDist)
total_time <- Sys.time() - start_time





