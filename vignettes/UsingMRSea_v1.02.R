## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
require(knitcitations)
cleanbib()
biblio <- read.bibtex("newref.bib")
cite_options(citation_format = 'pandoc', cite.style = 'authoryear', max.names = 1, longnamesfirst=FALSE)
knitr::opts_chunk$set(fig=TRUE, warning=FALSE, message=FALSE, eval=TRUE, comment = '')

## ----message=FALSE------------------------------------------------------------
#devtools::load_all(path='../../MRSea/')
require(MRSea)
# we will use the dataset with a known re-distribution of animals
data(dis.data.re)
dis.data<-dis.data.re
require(mrds) # distance sampling package
result <- ddf(dsmodel=~mcds(key="hn", formula=~1),
              data = dis.data, method="ds", 
              meta.data=list(width=250))

## ----dist, cache=TRUE, results='hide', warning=FALSE, message=FALSE-----------
# create.NHAT and create.count.data are MRSea functions to adjust the 
# sightings for the detection function estimated above.
dis.data <- create.NHAT(dis.data,result)
count.data <- create.count.data(dis.data)

## ----message=FALSE, warning=FALSE---------------------------------------------
data <- count.data
data$response <- round(data$NHAT)
attach(data)
fullModel <- glm(response ~ as.factor(season) + as.factor(impact) +
                   depth + x.pos + y.pos, family = poisson, data = data)

## ----message=FALSE------------------------------------------------------------
require(splines)
fullModel <- glm(response ~ as.factor(season) + as.factor(impact) +
                   bs(depth, knots = mean(depth)) + x.pos + y.pos, 
                 family = poisson,data = data)

## -----------------------------------------------------------------------------
data$blockid <- paste(data$transect.id, data$season, data$impact,sep = "")

## -----------------------------------------------------------------------------
data("nysted.predictdata")  # contains predict.data
# This is a spatial grid for making predictions.  All covariates in 
# final model must be in this data frame and the naming must be the 
# same as for the data
predictData <- nysted.predictdata
range(data$depth)
range(predictData$depth)

## ----message=FALSE------------------------------------------------------------
initialModel <- glm(response ~ as.factor(season) + as.factor(impact) 
                    + offset(log(area)), family = "quasipoisson", 
                    data = data)

## -----------------------------------------------------------------------------
salsa1dlist <- list(fitnessMeasure = "cv.gamMRSea", 
                    minKnots_1d = 2,
                    maxKnots_1d = 5, 
                    startKnots_1d = 1, 
                    degree = 2, 
                    maxIterations = 10,
                    gaps = c(0), 
                    cv.opts=list(cv.gamMRSea.seed=1, K=10))

## ----message=FALSE, warning=FALSE, echo=FALSE, results='hide'-----------------
# run SALSA
require(MuMIn)
salsa1dOutput <- runSALSA1D(initialModel, salsa1dlist, c("depth"),
                      predictionData=predictData, datain=data, removal=TRUE, 
                      panelid = data$blockid)

## ----eval=FALSE---------------------------------------------------------------
#  # run SALSA
#  salsa1dOutput <- runSALSA1D(initialModel, salsa1dlist, c("depth"),
#                        predictionData=predictData, datain=data, removal=TRUE,
#                        panelid = data$blockid)

## -----------------------------------------------------------------------------
summary(salsa1dOutput$bestModel)

## -----------------------------------------------------------------------------
# How many knots were chosen for depth?
salsa1dOutput$splineParams[[2]]$knots
# ~~~~~~~~~~~~~~~~~~~~~~~

## ----knotgrid, message=FALSE, fig=TRUE, fig.align='center', fig.width=9, fig.height=6, cache=TRUE----
knotgrid<- getKnotgrid(coordData = cbind(data$x.pos, data$y.pos), numKnots = 300)
#
# write.csv(knotgrid, file='knotgrid_fullanalysis.csv', row.names=F)
# ~~~~~~~~~~~~~~~~~~~~~~~

## -----------------------------------------------------------------------------
# make distance matrices for datatoknots and knottoknots
distMats <- makeDists(cbind(data$x.pos, data$y.pos), knotgrid)

# ~~~~~~~~~~~~~~~~~~~~~~~

# make parameter set for running salsa2d
salsa2dlist<-list(fitnessMeasure = 'cv.gamMRSea', 
                  knotgrid = knotgrid,
                  startKnots=10, 
                  minKnots=4, 
                  maxKnots=15, 
                  gap=0, 
                  interactionTerm="as.factor(impact)", 
                  cv.opts=list(cv.gamMRSea.seed=1, K=10))

## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'-----------------
salsa2dOutput<-runSALSA2D(salsa1dOutput$bestModel, 
                          salsa2dlist, 
                          d2k=distMats$dataDist,
                          k2k=distMats$knotDist)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  salsa2dOutput<-runSALSA2D(salsa1dOutput$bestModel,
#                            salsa2dlist,
#                            d2k=distMats$dataDist,
#                            k2k=distMats$knotDist)

## -----------------------------------------------------------------------------
require(ggplot2)
ggplot() + geom_point(data=data, aes(x=x.pos, y=y.pos), colour='grey') +
  theme_bw() + 
  geom_point(data=data.frame(knotgrid), aes(X1, X2), col='blue') + 
  geom_point(data=data.frame(knotgrid)[salsa2dOutput$aR[[1]],], aes(X1, X2), col='darkgreen', size=4) 


## ----acfplot, fig.cap='ACF plot showing correlation in each block (grey lines), and the mean correlation by lag across blocks (red line).'----
runACF(block = data$blockid, model = salsa2dOutput$bestModel, suppress.printout=TRUE)

## -----------------------------------------------------------------------------
simData<-generateNoise(n=500, response=fitted(salsa2dOutput$bestModel), family='poisson', d=summary(salsa2dOutput$bestModel)$dispersion)
empdist<-getEmpDistribution(500, simData, salsa2dOutput$bestModel, data=data,dots=FALSE)
runsTest(residuals(salsa2dOutput$bestModel, type='pearson'),emp.distribution=empdist)

## -----------------------------------------------------------------------------
set.seed(1)
# 2D model
cv.gamMRSea(data=data, modelobject = salsa2dOutput$bestModel, K=10)$delta[2]
# salsa2dOutput$fitStat # alternatively

# 1D model
cv.gamMRSea(data=data, modelobject = salsa1dOutput$bestModel, K=10)$delta[2]

# intitial model
cv.gamMRSea(data=data, modelobject = initialModel, K=10)$delta[2]

# full glm model with simple depth smooth
cv.gamMRSea(data=data, modelobject = fullModel, K=10)$delta[2]


## -----------------------------------------------------------------------------
anova(salsa2dOutput$bestModel)

## ----fig=TRUE, fig.align='center', fig.width=6, fig.height=4, message=FALSE----
par(mfrow=c(2,2))
runPartialPlots(model = salsa2dOutput$bestModel, data = data, 
                factorlist.in = c('season', 'impact'), 
                varlist.in = 'depth', showKnots = T, 
                includeB0 = TRUE)

## ----fig=TRUE, fig.align='center', fig.width=6, fig.height=4, message=FALSE----
par(mfrow=c(2,2))
runPartialPlots(model = salsa2dOutput$bestModel, data = data, 
                factorlist = c('season', 'impact'), varlist = 'depth', 
                showKnots = T, type='link', 
                includeB0 = TRUE)

## -----------------------------------------------------------------------------
preddist<-makeDists(cbind(predictData$x.pos, predictData$y.pos), 
                 knotgrid,knotmat=FALSE)$dataDist


# make predictions on response scale
preds<-predict.gamMRSea(newdata = predictData, g2k = preddist, object = salsa2dOutput$bestModel)

## ----fig=TRUE, fig.align='center', fig.width=9, fig.height=4------------------
require(RColorBrewer)
predictData$preds<-preds[,1]
ggplot() + geom_tile(data=predictData, aes(x.pos, y.pos, fill=preds), height=0.5, width=0.5) +
  facet_wrap(~impact + season, nrow=2) + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette = "Spectral",name="Animal Counts")

## ----boots, warning=FALSE, message=FALSE, results='hide'----------------------
dis.data$seasonimpact <- paste(dis.data$season, dis.data$impact)

bootPreds<-do.bootstrap.cress.robust(model.obj = salsa2dOutput$bestModel, 
                                     predictionGrid = predictData,
                                     g2k=preddist,
                                     B = 100, 
                                     robust=TRUE)

## -----------------------------------------------------------------------------
#load('predictionboot.RData')
cis <- makeBootCIs(bootPreds)

## -----------------------------------------------------------------------------
differences <- getDifferences(beforePreds = 
                      bootPreds[predictData$impact == 0, ],
                      afterPreds = bootPreds[predictData$impact == 1, ])

## ----fig=TRUE, fig.align='center', fig.width=9, fig.height=6------------------
mediandiff <- differences$mediandiff
# The marker for each after - before difference:
# positive ('1') and negative ('-') significant differences
marker <- differences$significanceMarker
par(mfrow = c(1, 1))
quilt.plot(predictData$x.pos[predictData$impact == 0], 
           predictData$y.pos[predictData$impact == 0],
           mediandiff, asp = 1, nrow = 104, ncol = 55)
# add + or - depending on significance of cells. Just
# requires one significance out of all to be allocated
points(predictData$x.pos[predictData$impact == 0][marker == 1],
       predictData$y.pos[predictData$impact == 0][marker == 1],
       pch = "+", col = "darkgrey", cex = 0.75)
points(predictData$x.pos[predictData$impact == 0][marker == (-1)],
       predictData$y.pos[predictData$impact == 0][marker == (-1)],
       col = "darkgrey", cex = 0.75)
points(681417.3/1000, 6046910/1000, cex = 3, pch = "*", lwd = 1, col = "grey")

## ----fig.width=9, fig.height=6------------------------------------------------
require(dplyr)
diffdata<-data.frame(predictData[predictData$impact==0,], mediandiff, marker)
diffdata_s1<-filter(diffdata, season==1)

wf<-data.frame(x=(681417.3/1000), y= (6046910/1000))

ggplot() + geom_tile(data=diffdata_s1, aes(x=x.pos, y=y.pos, fill=mediandiff), 
                     height=0.5, width=0.5) + 
  geom_point(data=filter(diffdata_s1, marker == 1), aes(x=x.pos, y=y.pos), 
             shape=3, colour='darkgrey', size=1) +
  geom_point(data=filter(diffdata_s1, marker == -1), aes(x=x.pos, y=y.pos), 
             shape=1, colour='darkgrey', size=1.5) +
  theme_bw() + coord_equal() +  
  scale_fill_distiller(palette = "Spectral",name="Difference") + 
  geom_point(data=wf, aes(x, y), shape=8, size=4)
  

## ---- fig.width=9, fig.height=7-----------------------------------------------
ggplot() + geom_tile(data=diffdata, aes(x=x.pos, y=y.pos, fill=mediandiff), 
                     height=0.5, width=0.5) + 
  geom_point(data=filter(diffdata, marker == 1), aes(x=x.pos, y=y.pos), 
             shape=3, colour='darkgrey', size=1) +
  geom_point(data=filter(diffdata, marker == -1), aes(x=x.pos, y=y.pos), 
             shape=1, colour='darkgrey', size=1.5) +
  theme_bw() + coord_equal() + facet_wrap(~season) + 
  scale_fill_distiller(palette = "Spectral",name="Difference") + 
  geom_point(data=wf, aes(x, y), shape=8, size=4)

## -----------------------------------------------------------------------------
finalmod<-salsa2dOutput$bestModel

## ---- fig.width=9, fig.height=6-----------------------------------------------
runDiagnostics(finalmod)

## -----------------------------------------------------------------------------
timeInfluenceCheck(finalmod, id = data$blockid)

## ---- fig.width=9, fig.height=6-----------------------------------------------
runInfluence(finalmod, id = data$blockid)

## ----fig.width=9, fig.height=7------------------------------------------------
plotCumRes(model = finalmod, varlist = 'depth')

## -----------------------------------------------------------------------------
fullModel.gamMRSea <- make.gamMRSea(fullModel, panelid = data$blockid, gamMRSea = TRUE)
summary(fullModel.gamMRSea)

## ---- eval=FALSE--------------------------------------------------------------
#  m1.robustse <- make.gamMRSea(m1, panelid =  data$blockid)

