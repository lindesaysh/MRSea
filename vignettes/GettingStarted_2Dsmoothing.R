## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
knitr::opts_chunk$set(fig=TRUE, warning=FALSE, message=FALSE, 
                      eval=TRUE, cache=FALSE,
                      comment = '#>', collapse=TRUE, dev='png')

## ----message=FALSE, warning=FALSE---------------------------------------------
library(dplyr)
library(MRSea)
library(ggplot2)

## -----------------------------------------------------------------------------
# load the data
data("nysted.analysisdata")
wfdata <- filter(nysted.analysisdata, impact==0, season==1)
# load the prediction grid
data("nysted.predictdata")
preddata <- filter(nysted.predictdata, impact==0, season==1)

## ----message=FALSE, warning=FALSE---------------------------------------------
ggplot(wfdata) + geom_tile(aes(x=x.pos, y=y.pos, fill=response, height=sqrt(area), width=sqrt(area))) +
  scale_fill_distiller(palette = "Spectral",name="No. Birds") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()

## ----message=FALSE------------------------------------------------------------
initialModel <- glm(response ~ 1 + offset(log(area)), family = "quasipoisson", 
                    data = wfdata)

## ----knotgrid, fig=TRUE, fig.align='center', fig.width=9, fig.height=6, fig.cap="Plot showing the candidate knot locations in red and the raw data locations in black."----
set.seed(123)
knotgrid<- getKnotgrid(coordData = cbind(wfdata$x.pos, wfdata$y.pos),
                       numKnots = 300,
                       plot = TRUE)

## -----------------------------------------------------------------------------
distMats <- makeDists(cbind(wfdata$x.pos, wfdata$y.pos), knotgrid)

## -----------------------------------------------------------------------------
# make parameter set for running salsa2d
salsa2dlist<-list(fitnessMeasure = 'QBIC',
                  knotgrid = knotgrid,
                  startKnots=10,
                  minKnots=4,
                  maxKnots=15,
                  gap=0)

## ----echo=TRUE, message=FALSE, warning=FALSE, results='hide'------------------
salsa2dOutput<-runSALSA2D(model = initialModel,
                          salsa2dlist = salsa2dlist,
                          d2k=distMats$dataDist,
                          k2k=distMats$knotDist,
                          suppress.printout = TRUE)

## ----eval=TRUE----------------------------------------------------------------
summary(salsa2dOutput$bestModel)

## ----eval=TRUE----------------------------------------------------------------
# How many knots were chosen for depth?
salsa2dOutput$bestModel$splineParams[[1]]$knotPos

## ----message=FALSE, warning=FALSE---------------------------------------------
ggplot(wfdata) + geom_tile(aes(x=x.pos, y=y.pos, fill=fitted(salsa2dOutput$bestModel), height=sqrt(area), width=sqrt(area))) +
  scale_fill_distiller(palette = "Spectral",name="No. Birds") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()

## -----------------------------------------------------------------------------
preddist<-makeDists(cbind(preddata$x.pos, preddata$y.pos),
                 knotgrid, knotmat=FALSE)$dataDist


# make predictions on response scale
preds<-predict(newdata = preddata,
               g2k = preddist,
               object = salsa2dOutput$bestModel)

## ----message=FALSE, warning=FALSE---------------------------------------------
ggplot(preddata) + geom_tile(aes(x=x.pos, y=y.pos, fill=preds, height=sqrt(area), width=sqrt(area))) +
  scale_fill_distiller(palette = "Spectral",name="No. Birds") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()

## ----eval=FALSE---------------------------------------------------------------
#  salsa2dlist<-list(fitnessMeasure = 'QBIC',
#                    knotgrid = knotgrid,
#                    startKnots=10,
#                    minKnots=4,
#                    maxKnots=15,
#                    gap=0,
#                    interactionTerm="as.factor(impact)")

