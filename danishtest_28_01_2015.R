rm(list=ls())

require(MRSea)
# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~

data(dis.data.no)
dis.data<-dis.data.no
require(mrds)
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
require(MRSea)

data <- count.data
data$response<-data$NHAT
attach(data)
fullModel <- glm(round(NHAT) ~ as.factor(season) + as.factor(impact) +
                   depth + x.pos + y.pos, family = poisson, data = data)
require(car)
vif(fullModel)
checkfactorlevelcounts(c("season", "impact"), data, data$response)
knots <- mean(depth) # must be specified as an object
require(splines)
fullModel <- glm(round(NHAT) ~ as.factor(season) + as.factor(impact) +
                   bs(depth, knots = knots) + x.pos + y.pos, family = poisson,
                 data = data)

fullModel <- glm(NHAT ~ as.factor(season) + as.factor(impact) +
                   bs(depth, knots = knots) + x.pos + y.pos, family = quasipoisson,
                 data = data)

# for correlated data:
data$blockid <- paste(data$transect.id, data$season, data$impact,sep = "")
data$foldid <- getCVids(data = data, folds = 5, block = 'blockid')
# for uncorrelated data@
data$foldid<- getCVids(data=data, folds=5)

data$response <- round(data$NHAT)
salsa1dlist <- list(fitnessMeasure = "QAIC", minKnots_1d = 2,
                    maxKnots_1d = 5, startKnots_1d = 1, degree = 2, maxIterations = 10,
                    gaps = c(1))

data(predict.data.no)  # contains predict.data
# This is a spatial grid for making predictions.  All covariates in final model must be in this data frame and the naming must be the same as for the data
predictData <- predict.data.no
initialModel <- glm(response ~ as.factor(season) + as.factor(impact) +offset(log(area)), family = "quasipoisson", data = data)

# run SALSA
salsa1dOutput <- runSALSA1D_withremoval(initialModel, salsa1dlist, c("depth"), predictionData=predictData, datain=data, removal=TRUE)
salsa1dOutput$splineParams[[2]]$knots
splineParams<-salsa1dOutput$splineParams
# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~

knotgrid<- getKnotgrid(coordData = cbind(data$x.pos, data$y.pos))
# this may take while and could be different every time you run it so I suggest saving the knotgrid as a file.
#
# write.csv(knotgrid, file='knotgrid_fullanalysis.csv', row.names=F)
# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~

# make distance matrices for datatoknots and knottoknots
distMats <- makeDists(cbind(data$x.pos, data$y.pos), na.omit(knotgrid))

r_seq <- getRadiiChoices(10, distMats$dataDist)

# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~

# make parameter set for running salsa2d
salsa2dlist<-list(fitnessMeasure = 'QICb', knotgrid = knotgrid, knotdim=c(100,100), startKnots=6, minKnots=4, maxKnots=10, r_seq=r_seq, gap=4000, interactionTerm="as.factor(impact)")
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
require(geepack)
geeModel<- geeglm(formula(salsa2dOutput_k6$bestModel), data=data, family=poisson, id=blockid)
dists<-makeDists(cbind(predict.data.re$x.pos, predict.data.re$y.pos), na.omit(knotgrid),knotmat=FALSE)$dataDist

# make predictions on response scale
preds<-predict.cress(predict.data.re, splineParams, dists, geeModel)

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

getCV_CReSS(data, geeModel, splineParams = splineParams)


  b<- LocalRadialFunction(splineParams[[1]]$radiusIndices, dists, radii,aR)
# 
# b<- LocalRadialFunction(rep(25, length(aR)), dists, r_seq, 
#                         splineParams[[1]]$invInd[splineParams[[1]]$knotPos])
# 
# 
par(mfrow=c(2,2))
for(i in 1:(length(splineParams[[1]]$knotPos)){
  print(i)
  quilt.plot(predict.data.re$x.pos, predict.data.re$y.pos, b[,i], asp=1)
  points(knotgrid[splineParams[[1]]$knotPos,], pch=20, cex=2)
}

distMats <- makeDists(cbind(data$x.pos, data$y.pos), na.omit(knotgrid))

mindist<-mean(apply(distMats$dataDist, 2, min))
maxdist<-mean(apply(distMats$dataDist, 2, max))
r_seq=exp(seq(log(1/mindist)*1.5, log(1/maxdist)*2,length=10))
r_seq
r_seq<-getRadiiChoices(10, distMats$dataDist)
b<- LocalRadialFunction(1:length(r_seq), dists, r_seq,rep(aR[1], length(r_seq)))
summary(b)
par(mfrow=c(2,2))
for(i in 1:(length(r_seq))){
  print(i)
  quilt.plot(predict.data.re$x.pos, predict.data.re$y.pos, b[,i], asp=1, zlim=c(0,1))
  points((knotgrid)[splineParams[[1]]$knotPos,], pch=20, cex=2)
}

# # 
# dists<-makeDists(cbind(count.data$x.pos, count.data$y.pos), na.omit(knotgrid.off))$dataDist
# quilt.plot(count.data$x.pos, count.data$y.pos, fitted(geeModel))
