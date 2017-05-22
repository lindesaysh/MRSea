## ----echo=FALSE, message=FALSE, warning=FALSE----------------------------
require(knitcitations)
cleanbib()
biblio <- read.bibtex("newref.bib")
cite_options(citation_format = 'pandoc', cite.style = 'authoryear', max.names = 1, longnamesfirst=FALSE)

## ----message=FALSE-------------------------------------------------------
devtools::load_all(pkg='../../MRSea')
# we will use the dataset with a known re-distribution of animals
data(dis.data.re)
dis.data<-dis.data.re
require(mrds) # distance sampling package
result <- ddf(dsmodel=~mcds(key="hn", formula=~1),
              data = dis.data, method="ds", 
              meta.data=list(width=250))

## ----dist, cache=TRUE, results='hide', warning=FALSE, message=FALSE------
# create.NHAT and create.count.data are MRSea functions to adjust the 
# sightings for the detection function estimated above.
dis.data <- create.NHAT(dis.data,result)
count.data <- create.count.data(dis.data)

## ------------------------------------------------------------------------
data <- count.data
data$response <- round(data$NHAT)
attach(data)
fullModel <- glm(response ~ as.factor(season) + as.factor(impact) +
                   depth + x.pos + y.pos, family = poisson, data = data)

## ----message=FALSE-------------------------------------------------------
require(splines)
fullModel <- glm(response ~ as.factor(season) + as.factor(impact) +
                   bs(depth, knots = mean(depth)) + x.pos + y.pos, 
                 family = poisson,data = data)

## ------------------------------------------------------------------------
# for correlated data:
data$blockid <- paste(data$transect.id, data$season, data$impact,sep = "")
data$foldid <- getCVids(data = data, folds = 5, block = 'blockid')

## ----eval=FALSE----------------------------------------------------------
#  # for uncorrelated data```
#  data$foldid<- getCVids(data=data, folds=5)

## ------------------------------------------------------------------------
salsa1dlist <- list(fitnessMeasure = "cv.gamMRSea", minKnots_1d = 2,maxKnots_1d = 5, 
                    startKnots_1d = 1, degree = 2, maxIterations = 10,
                    gaps = c(0), seed.in=1)

## ------------------------------------------------------------------------
data("nysted.predictdata")  # contains predict.data
# This is a spatial grid for making predictions.  All covariates in 
# final model must be in this data frame and the naming must be the 
# same as for the data
predictData <- nysted.predictdata
range(data$depth)
range(predictData$depth)

## ----message=FALSE-------------------------------------------------------
initialModel <- glm(response ~ as.factor(season) + as.factor(impact) 
                    + offset(log(area)), family = "quasipoisson", 
                    data = data)

## ----message=FALSE, warning=FALSE, echo=FALSE, results='hide'------------
# run SALSA
require(MuMIn)
salsa1dOutput <- runSALSA1D(initialModel, salsa1dlist, c("depth"),
                      predictionData=predictData, datain=data, removal=TRUE, panelid = data$blockid)

## ----eval=FALSE----------------------------------------------------------
#  # run SALSA
#  salsa1dOutput <- runSALSA1D(initialModel, salsa1dlist, c("depth"),
#                        predictionData=predictData, datain=data, removal=TRUE, panelid = data$blockid)

## ------------------------------------------------------------------------
summary(salsa1dOutput$bestModel)

## ------------------------------------------------------------------------
# How many knots were chosen for depth?
salsa1dOutput$splineParams[[2]]$knots
# ~~~~~~~~~~~~~~~~~~~~~~~

## ----knotgrid, message=FALSE, fig=TRUE, fig.align='center', fig.width=9, fig.height=6, cache=TRUE----
knotgrid<- getKnotgrid(coordData = cbind(data$x.pos, data$y.pos), numKnots = 300)
#
# write.csv(knotgrid, file='knotgrid_fullanalysis.csv', row.names=F)
# ~~~~~~~~~~~~~~~~~~~~~~~

## ------------------------------------------------------------------------
# make distance matrices for datatoknots and knottoknots
distMats <- makeDists(cbind(data$x.pos, data$y.pos), knotgrid)

# ~~~~~~~~~~~~~~~~~~~~~~~

# make parameter set for running salsa2d
salsa2dlist<-list(fitnessMeasure = 'cv.gamMRSea', knotgrid = knotgrid, 
                 startKnots=5, minKnots=4, maxKnots=12, gap=0, 
                 interactionTerm="as.factor(impact)", cv.gamMRSea.seed=1)

## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'------------
salsa2dOutput<-runSALSA2D(salsa1dOutput$bestModel, salsa2dlist, 
                      d2k=distMats$dataDist,k2k=distMats$knotDist)

## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  salsa2dOutput<-runSALSA2D(salsa1dOutput$bestModel, salsa2dlist,
#                               d2k=distMats$dataDist, k2k=distMats$knotDist)

## ------------------------------------------------------------------------
plot(data$x.pos, data$y.pos, col="grey", pch=16,
    xlab="X", ylab="Y", asp=1)
points(knotgrid, pch=16, col=4)
points(knotgrid[salsa2dOutput$aR[[1]],], 
       col="darkgreen", pch=16, cex=2)  

## ----acfplot, fig.cap='ACF plot showing correlation in each block (grey lines), and the mean correlation by lag across blocks (red line).'----
runACF(block = data$blockid, model = salsa2dOutput$bestModel, suppress.printout=TRUE)

