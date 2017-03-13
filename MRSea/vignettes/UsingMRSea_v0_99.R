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
salsa1dlist <- list(fitnessMeasure = "QBIC", minKnots_1d = 2,maxKnots_1d = 5, 
                    startKnots_1d = 1, degree = 2, maxIterations = 10,
                    gaps = c(0))

## ------------------------------------------------------------------------
data(predict.data.re)  # contains predict.data
# This is a spatial grid for making predictions.  All covariates in 
# final model must be in this data frame and the naming must be the 
# same as for the data
predictData <- predict.data.re
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

## ------------------------------------------------------------------------
# make distance matrices for datatoknots and knottoknots
distMats <- makeDists(cbind(data$x.pos, data$y.pos), na.omit(knotgrid))

r_seq <- getRadiiChoices(numberofradii=10, distMatrix=distMats$dataDist)

# ~~~~~~~~~~~~~~~~~~~~~~~

# make parameter set for running salsa2d
salsa2dlist<-list(fitnessMeasure = 'QBIC', knotgrid = knotgrid, 
                  knotdim=c(100,100), startKnots=10, minKnots=4,
                  maxKnots=12, r_seq=r_seq, gap=0, 
                  interactionTerm="as.factor(impact)")

