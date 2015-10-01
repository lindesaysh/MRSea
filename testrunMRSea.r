library(MRSea)

setwd('Work/MarineScotland/packagetesting/')

data(dis.data.no)
dis.data<-dis.data.no
result <- ddf(dsmodel=~mcds(key="hn", formula=~1),
              data = dis.data, method="ds", meta.data=list(width=250))
#result.hr <- ddf(dsmodel=~mcds(key="hr", formula=~1),
#                 data = dis.data, method="ds", meta.data=list(width=250))
# Using a half-normal function with impact as a covariate
#result.imp <- ddf(dsmodel=~mcds(key="hn", formula=~1+impact),
#                  data = dis.data, method="ds", meta.data=list(width=250))
dis.data <- create.NHAT(dis.data,result)
count.data <- create.count.data(dis.data)
data <- count.data
data$response<-data$NHAT
attach(data)
fullModel <- glm(NHAT ~ as.factor(season) + as.factor(impact) +
                   depth + x.pos + y.pos, family = poisson, data = data)
vif(fullModel)
checkfactorlevelcounts(c("season", "impact"), data, data$response)
knots <- mean(depth) # must be specified as an object
fullModel <- glm(NHAT ~ as.factor(season) + as.factor(impact) +
                   bs(depth, knots = knots) + x.pos + y.pos, family = quasipoisson,
                 data = data)
runs.test(residuals(fullModel, type = "pearson"), alternative = c("two.sided"))
data$blockid <- paste(data$transect.id, data$season, data$impact,
                      sep = "")
runACF(data$blockid, fullModel, store = F)
data$response <- data$NHAT
salsa1dlist <- list(fitnessMeasure = "QICb", minKnots_1d = 2,
                    maxKnots_1d = 20, startKnots_1d = 2, degree = 2, maxIterations = 10,
                    gaps = c(1))

data(predict.data.no)  # contains predict.data
predictData <- predict.data.no
initialModel <- glm(response ~ as.factor(season) + as.factor(impact) +offset(area), family = "quasipoisson", data = data)
# run SALSA
salsa1dOutput <- runSALSA1D(initialModel, salsa1dlist, c("depth"),
                            predictionData=predictData)
salsa1dOutput$splineParams[[2]]$knots



#data(knotgrid.off) # knotgrid <- read.csv("data/Danishknotgrid_salsa.csv")
#knotgrid <- knotgrid.off
require(splancs)
XGridPoints<- seq(min(data$x.pos),
                  (max(data$x.pos)), length=100)
YGridPoints<- seq(min(data$y.pos),
                  (max(data$y.pos)), length=100)
Grid<- as.points(expand.grid(x=XGridPoints,y=YGridPoints))
dataXandY<- as.points(x.pos=data$x.pos,y.pos=data$y.pos)


# ~~~~~~~~~~~~~~~~~~~~~~

#visualise the candidate knots
par(mfrow=c(1,1))
plot(data$x.pos, data$y.pos, col="grey", pch=16)
points(Grid, pch=20, col=2)
#points(data$x.pos, data$y.pos, col="grey", pch=16)
dim(Grid)

# ~~~~~~~~~~~~~~~~~~~~~~
# remove knots inside land and way outside boundary
# ~~~~~~~~~~~~~~~~~~~~~~
#find the grid point closest to each candidate knot on the grid
d<- c()
for(i in 1:nrow(data)){
  d[i]<- which.min(sqrt((data$x.pos[i]-Grid[,1])**2+(data$y.pos[i]-Grid[,2])**2))
}

RowsToKeep<- rep(0,nrow(Grid))
RowsToKeep[d]<-1

GridPosIncludingNAs<- cbind(ifelse(RowsToKeep!=1, NA, RowsToKeep*Grid[,1]), ifelse(RowsToKeep!=1, NA, RowsToKeep*Grid[,2]))
Grid <- as.data.frame(GridPosIncludingNAs)

knotgrid<-Grid

points(knotgrid, pch=20)

dim(knotgrid)
dim(na.omit(knotgrid))

require(fields)
naid<- which(is.na(knotgrid))
spaceid<-cover.design(R = na.omit(knotgrid), nd = 400, nruns = 5)$best.id

realid<-(1:nrow(knotgrid))[-naid]

points(knotgrid[realid[spaceid],], pch=20, col='purple')

knotgrid[realid[-spaceid],]<- c(NA, NA)

plot(data$x.pos, data$y.pos, pch=20, cex=0.2)
points(knotgrid, pch=20, cex=0.5, col='red')

#write.csv(knotgrid, file='Data/knotgrid_fullanalysis.csv', row.names=F)

# make distance matrices for datatoknots and knottoknots
distMats <- makeDists(cbind(data$x.pos, data$y.pos), na.omit(knotgrid))

#r_seq <- getRadiiChoices(8, distMats$dataDist)

mindist<-mean(apply(distMats$dataDist, 2, min))
maxdist<-mean(apply(distMats$dataDist, 2, max))

#r_seq=exp(seq(log(1/mindist), log(1/maxdist),length=50))[10:50]
r_seq=1/(seq((mindist), (maxdist),length=50))


salsa2dlist <- list(fitnessMeasure = "QICb", knotgrid = knotgrid,startKnots = 6, minKnots = 4, maxKnots = 20, r_seq = r_seq,gap = 1, interactionTerm = "as.factor(impact)")

splineParams <- salsa1dOutput$splineParams
salsa2dOutput_k6 <- runSALSA2D(salsa1dOutput$bestModel, salsa2dlist,
                               distMats$dataDist, distMats$knotDist, splineParams = splineParams)
data$foldid <- getCVids(data, folds = 5, block = "blockid")
cv1 <- getCV_CReSS(data, salsa1dOutput$bestModel, salsa1dOutput$splineParams)
baseModel <- salsa2dOutput_k6$bestModel
# update spline parameter object
splineParams <- salsa2dOutput_k6$splineParams

radiusIndices <- splineParams[[1]]$radiusIndices
dists <- splineParams[[1]]$dist
radii <- splineParams[[1]]$radii
aR <- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]
baseModel <- update(baseModel, . ~ .)
vif(baseModel)

geeModel <- geeglm(formula(baseModel), data = data, family = poisson,id = blockid)
getPvalues(geeModel, varlist = c("depth"), factorlist = c("season","impact"))
model <- update(geeModel, . ~ . - as.factor(impact))

# reshow p-values
getPvalues(model, varlist = c("depth"), factorlist = c("season","impact"))
# ~~~~~~~~~~~~~~~~~~~~~~~
runPartialPlots(geeModel, data, factorlist = c("season", "impact"),
                varlist = c("depth"))
# ~~~~~~~~~~~~~~~~~~~~~~~
runDiagnostics(geeModel)
# ~~~~~~~~~~~~~~~~~~~~~~~
runsTest(geeModel, c("observationhour"), dists, splineParams,
         label = "geemodel")
# ~~~~~~~~~~~~~~~~~~~~~~~
timeInfluenceCheck(geeModel, data$blockid, dists, splineParams)
# influence plots (covratio and press statistics)
influence<-runInfluence(geeModel, data$blockid, dists, splineParams) 

# ~~~~~~~~~~~~~~~~~~~~~~~
resids<-fitted(geeModel) - data$response
dims<-getPlotdimensions(data$x.pos, data$y.pos, 1000, 1000)

par(mfrow=c(1,2), mar=c(3,3,3,5))
quilt.plot(data$x.pos[data$impact==0], data$y.pos[data$impact==0], resids[data$impact==0], asp=1, ncol=dims[2], nrow=dims[1], zlim=c(-2.2, 2.2))
quilt.plot(data$x.pos[data$impact==1], data$y.pos[data$impact==1], resids[data$impact==1], asp=1, ncol=dims[2], nrow=dims[1], zlim=c(-2.2, 2.2))

# ~~~~~~~~~~~~~~~~~~~~~~~
dists<-makeDists(cbind(predictData$x.pos, predictData$y.pos), na.omit(knotgrid), knotmat=FALSE)$dataDist
predslink<- predict(baseModel, predictData, type='link')
preds<-exp(predslink)

# ~~~~~~~~~~~~~~~~~~~~~~~
try1<- aggregate(preds, by = list(predictData$segment.id, predictData$impact), mean)
maxlim<-max(abs(min(try1$x)), max(try1$x))

# prediction plot
par(mfrow=c(1,2), mar=c(3,3,3,5))
quilt.plot(predictData$x.pos[predictData$impact==0], predictData$y.pos[predictData$impact==0], preds[predictData$impact==0], asp=1, nrow=104, ncol=55, zlim=c(0, maxlim))
points(0,0,pch=20, col='grey', cex=3)
points(data$x.pos[which(data$GridCode=='e7')[1]], data$y.pos[which(data$GridCode=='e7')[1]],cex=2, pch='*', lwd=2, col='grey')
points(knotgrid[splineParams[[1]]$knotPos,], pch=2, col='darkgrey', cex=1.5)
quilt.plot(predictData$x.pos[predictData$impact==1], predictData$y.pos[predictData$impact==1], preds[predictData$impact==1], asp=1, nrow=104, ncol=55, zlim=c(0, maxlim))
points(knotgrid[splineParams[[1]]$knotPos,], pch=2, col='darkgrey', cex=1.5)


# ~~~~~~~~~~~~~~~~~~~~~~~
do.bootstrap.cress(data, predictData, ddf.obj=NULL, baseModel, splineParams, dists, B=250)
# read in bootstrap predictions
load('predictionboot.RData')
cis<-makeBootCIs(bootPreds)


# ~~~~~~~~~~~~~~~~~~~~~~~
try1<- aggregate(cis[,1], by = list(predictData$segment.id, predictData$impact), mean)
lowermax<-max(try1$x); lowermin<-min(try1$x)
try2<- aggregate(cis[,2], by = list(predictData$segment.id, predictData$impact), mean)
uppermax<-max(try2$x); uppermin<-min(try2$x)
lims<-c(0, max(lowermax, uppermax))
# make plot and save
par(mfrow=c(1,2), mar=c(3,3,3,5))
quilt.plot(predictData$x.pos[predictData$impact==0], predictData$y.pos[predictData$impact==0], cis[predictData$impact==0, 1], nrow=104, ncol=55, zlim=c(0,uppermax), asp=1)
quilt.plot(predictData$x.pos[predictData$impact==1], predictData$y.pos[predictData$impact==1], cis[predictData$impact==1, 1], nrow=104, ncol=55, zlim=c(0,uppermax), asp=1)

par(mfrow=c(1,2), mar=c(3,3,3,5))
quilt.plot(predictData$x.pos[predictData$impact==0], predictData$y.pos[predictData$impact==0], cis[predictData$impact==0, 2], nrow=104, ncol=55, zlim=c(0,uppermax), asp=1)
quilt.plot(predictData$x.pos[predictData$impact==1], predictData$y.pos[predictData$impact==1], cis[predictData$impact==1, 2], nrow=104, ncol=55, zlim=c(0,uppermax), asp=1)

# ~~~~~~~~~~~~~~~~~~~~~~~
differences<-getDifferences(beforePreds=bootPreds[predictData$impact==0,], afterPreds=bootPreds[predictData$impact==1,])
meandiff<-differences$meandiff
marker<-differences$significanceMarker
# make plot
par(mfrow=c(1,1))
quilt.plot(predictData$x.pos[predictData$impact==0], predictData$y.pos[predictData$impact==0], meandiff, asp=1, nrow=104, ncol=55)
# add + or - depending on significance of cells.  Just requires one significance out of all to be allocated
points(predictData$x.pos[predictData$impact==0][marker==1], predictData$y.pos[predictData$impact==0][marker==1], pch='+', col='darkgrey', cex=2)
points(predictData$x.pos[predictData$impact==0][marker==(-1)], predictData$y.pos[predictData$impact==0][marker==(-1)], col='darkgrey', cex=3)
points(0,0,pch=20, col='grey', cex=3)
points(data$x.pos[which(data$GridCode=='e7')[1]], data$y.pos[which(data$GridCode=='e7')[1]],cex=2, pch='*', lwd=2, col='grey')