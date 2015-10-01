rm(list=ls())

require(MRSea)
# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~

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
# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~

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
#runACF(data$blockid, fullModel, store = F)
data$response <- data$NHAT
salsa1dlist <- list(fitnessMeasure = "QICb", minKnots_1d = 2,
                    maxKnots_1d = 20, startKnots_1d = 1, degree = 2, maxIterations = 10,
                    gaps = c(5))

data(predict.data.no)  # contains predict.data
predictData <- predict.data.no
initialModel <- glm(response ~ as.factor(season) + as.factor(impact) +offset(log(area)), family = "quasipoisson", data = data)

# run SALSA
salsa1dOutput <- runSALSA1D(initialModel, salsa1dlist, c("depth"),
                            predictionData=predictData)
salsa1dOutput$splineParams[[2]]$knots
splineParams<-salsa1dOutput$splineParams
# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~
#data(knotgrid.off) # knotgrid <- read.csv("data/Danishknotgrid_salsa.csv")
#knotgrid <- knotgrid.off
require(splancs)
knotgrid<-read.csv('knotgrid_fullanalysis.csv')
# XGridPoints<- seq(min(data$x.pos),
#                   (max(data$x.pos)), length=100)
# YGridPoints<- seq(min(data$y.pos),
#                   (max(data$y.pos)), length=100)
# Grid<- as.points(expand.grid(x=XGridPoints,y=YGridPoints))
# dataXandY<- as.points(x.pos=data$x.pos,y.pos=data$y.pos)
# 
# #visualise the candidate knots
# par(mfrow=c(1,1))
# plot(data$x.pos, data$y.pos, col="grey", pch=16)
# points(Grid, pch=20, col=2)
# #points(data$x.pos, data$y.pos, col="grey", pch=16)
# dim(Grid)
# #find the grid point closest to each candidate knot on the grid
# d<- c()
# for(i in 1:nrow(data)){
#   d[i]<- which.min(sqrt((data$x.pos[i]-Grid[,1])**2+(data$y.pos[i]-Grid[,2])**2))
# }
# 
# RowsToKeep<- rep(0,nrow(Grid))
# RowsToKeep[d]<-1
# 
# GridPosIncludingNAs<- cbind(ifelse(RowsToKeep!=1, NA, RowsToKeep*Grid[,1]), ifelse(RowsToKeep!=1, NA, RowsToKeep*Grid[,2]))
# Grid <- as.data.frame(GridPosIncludingNAs)
# knotgrid<-Grid
# 
# points(knotgrid, pch=20)
# dim(knotgrid)
# dim(na.omit(knotgrid))
# 
# require(fields)
# naid<- which(is.na(knotgrid))
# spaceid<-cover.design(R = na.omit(knotgrid), nd = 400, nruns = 5)$best.id
# realid<-(1:nrow(knotgrid))[-naid]
# 
# points(knotgrid[realid[spaceid],], pch=20, col='purple')
# 
# knotgrid[realid[-spaceid],]<- c(NA, NA)
# 
# plot(data$x.pos, data$y.pos, pch=20, cex=0.2)
# points(knotgrid, pch=20, cex=0.5, col='red')
# 
# write.csv(knotgrid, file='knotgrid_fullanalysis.csv', row.names=F)
# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~

# make distance matrices for datatoknots and knottoknots
distMats <- makeDists(cbind(data$x.pos, data$y.pos), na.omit(knotgrid))

#r_seq <- getRadiiChoices(8, distMats$dataDist)

mindist<-mean(apply(distMats$dataDist, 2, min))
maxdist<-mean(apply(distMats$dataDist, 2, max))

#r_seq=exp(seq(log(1/mindist), log(1/maxdist),length=50))[10:50]
r_seq=1/(seq((mindist), (maxdist),length=50))
r_seq
# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~

# make parameter set for running salsa2d
salsa2dlist<-list(fitnessMeasure = 'QICb', knotgrid = knotgrid, knotdim=c(100,100), startKnots=15, minKnots=4, maxKnots=45, r_seq=r_seq, gap=4000, interactionTerm="as.factor(impact)")
salsa2dOutput_k6<-runSALSA2D(salsa1dOutput$bestModel, salsa2dlist, d2k=distMats$dataDist,
                             k2k=distMats$knotDist, splineParams=splineParams)

splineParams<-salsa2dOutput_k6$splineParams
# specify parameters for local radial function:
radiusIndices <- splineParams[[1]]$radiusIndices
dists <- splineParams[[1]]$dist
radii <- splineParams[[1]]$radii
aR <- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]
data$blockid<-paste(data$transect.id, data$season, data$impact, sep='')
# Re-fit the chosen model as a GEE (based on SALSA knot placement) and GEE p-values
geeModel<- geeglm(formula(salsa2dOutput_k6$bestModel), data=data, family=poisson, id=blockid)
dists<-makeDists(cbind(predict.data.re$x.pos, predict.data.re$y.pos), na.omit(knotgrid),
                 knotmat=FALSE)$dataDist

# make predictions on response scale
preds<-predict.cress(predict.data.re, splineParams, dists, geeModel, modelav = FALSE)

#preds<-predict(geeModel, predict.data.re, type='response')

par(mfrow=c(2,2))
quilt.plot(predict.data.re$x.pos[predict.data.re$impact==0], predict.data.re$y.pos[predict.data.re$impact==0], preds[predict.data.re$impact==0], nrow=104, ncol=55)

quilt.plot(predict.data.re$x.pos[predict.data.re$impact==1], predict.data.re$y.pos[predict.data.re$impact==1], preds[predict.data.re$impact==1], nrow=104, ncol=55)
points(na.omit(knotgrid)[aR,], pch=20, cex=2)

sum((preds-predict.data.re$truth)^2)
length(aR)

quilt.plot(predict.data.re$x.pos[predict.data.re$impact==0], predict.data.re$y.pos[predict.data.re$impact==0], predict.data.re$truth[predict.data.re$impact==0] - preds[predict.data.re$impact==0], nrow=104, ncol=55)

quilt.plot(predict.data.re$x.pos[predict.data.re$impact==1], predict.data.re$y.pos[predict.data.re$impact==1], predict.data.re$truth[predict.data.re$impact==1] - preds[predict.data.re$impact==1], nrow=104, ncol=55)
points(knotgrid[splineParams[[1]]$knotPos,], pch=20, cex=2)

data$foldid<-getCVids(data = data, folds = 5, block='blockid')
getCV_CReSS(data, geeModel, splineParams = splineParams)


  b<- LocalRadialFunction(splineParams[[1]]$radiusIndices, dists, radii,aR)
# 
# b<- LocalRadialFunction(rep(25, length(aR)), dists, r_seq, 
#                         splineParams[[1]]$invInd[splineParams[[1]]$knotPos])
# 
# 
par(mfrow=c(2,2))
for(i in 1:(length(coef(salsa2dOutput_k6$bestModel))-3)){
  print(i)
  quilt.plot(predict.data.re$x.pos, predict.data.re$y.pos, b[,i], asp=1)
  points(knotgrid[splineParams[[1]]$knotPos,], pch=20, cex=2)
}

# # 
# dists<-makeDists(cbind(count.data$x.pos, count.data$y.pos), na.omit(knotgrid.off))$dataDist
# quilt.plot(count.data$x.pos, count.data$y.pos, fitted(geeModel))
