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

## ----message=FALSE, warning=FALSE---------------------------------------------
ggplot() + geom_point(data=wfdata, aes(x=depth, y=response), alpha=1/5) +
  xlab("Sea Depth (m)") + ylab("Number of Birds")

## ----message=FALSE------------------------------------------------------------
initialModel <- glm(response ~ 1 + offset(log(area)), family = "quasipoisson", 
                    data = wfdata)

## -----------------------------------------------------------------------------
salsa1dlist <- list(fitnessMeasure = "QBIC", 
                    minKnots_1d = 1,
                    maxKnots_1d = 3, 
                    startKnots_1d = 1, 
                    degree = 2,
                    gaps = c(0))

## ----message=FALSE, warning=FALSE, echo=TRUE, results='hide'------------------
salsa1dOutput <- runSALSA1D(initialModel = initialModel, 
                            salsa1dlist = salsa1dlist,
                            varlist = c("depth"),
                            predictionData = nysted.predictdata, 
                            datain = wfdata,
                            suppress.printout = TRUE)

## ----eval=TRUE----------------------------------------------------------------
summary(salsa1dOutput$bestModel)

## ----eval=TRUE----------------------------------------------------------------
# How many knots were chosen for depth?
salsa1dOutput$bestModel$splineParams[[2]]$knots

## ----eval=TRUE, fig=TRUE, fig.align='center', fig.width=6, fig.height=4, message=FALSE----
runPartialPlots(model = salsa1dOutput$bestModel, data = wfdata, 
                varlist = 'depth', 
                showKnots = TRUE, 
                type='link', 
                includeB0 = TRUE)

## -----------------------------------------------------------------------------
runPartialPlots(model = salsa1dOutput$bestModel, data = wfdata, 
                varlist = 'depth', 
                showKnots = TRUE, type='response', 
                includeB0 = TRUE)

## ----message=FALSE------------------------------------------------------------
initialModel <- glm(response ~ 1 + as.factor(season) + offset(log(area)), family = "quasipoisson", 
                    data = nysted.analysisdata)

## -----------------------------------------------------------------------------
salsa1dlist <- list(fitnessMeasure = "QBIC", 
                    minKnots_1d = 1,
                    maxKnots_1d = 3, 
                    startKnots_1d = 1, 
                    degree = 2,
                    gaps = c(0))

salsa1dOutput.f <- runSALSA1D(initialModel = initialModel, 
                            salsa1dlist = salsa1dlist,
                            varlist = c("depth"),
                            predictionData = nysted.predictdata, 
                            datain = nysted.analysisdata,
                            suppress.printout = TRUE)

## ----eval=TRUE, fig=TRUE, fig.align='center', fig.width=6, fig.height=4, message=FALSE----
runPartialPlots(model = salsa1dOutput.f$bestModel, 
                data = nysted.analysisdata, 
                varlist.in = 'depth', 
                factorlist.in = "season",
                showKnots = TRUE, 
                type='link', 
                includeB0 = TRUE)

## ----message=FALSE------------------------------------------------------------
data(ns.data.re)
vpdata <- filter(ns.data.re, impact==0) %>%
  mutate(response = birds)
head(vpdata)

## -----------------------------------------------------------------------------
initialModel <- glm(response ~ 1 + as.factor(floodebb) + offset(log(area)), family = "quasipoisson", 
                    data = vpdata)

## -----------------------------------------------------------------------------
varlist <- c("observationhour", "x.pos")

salsa1dlist <- list(fitnessMeasure = "QBIC", 
                    minKnots_1d = rep(1, length(varlist)),
                    maxKnots_1d = rep(1, length(varlist)), 
                    startKnots_1d = rep(1, length(varlist)), 
                    degree = rep(2, length(varlist)),
                    gaps = rep(0, length(varlist)))

salsa1dOutput.multi <- runSALSA1D(initialModel = initialModel, 
                            salsa1dlist = salsa1dlist,
                            varlist = varlist, 
                            datain = vpdata,
                            suppress.printout = TRUE)

## ----eval=TRUE, fig=TRUE, fig.align='center', fig.width=6, fig.height=4, message=FALSE----
runPartialPlots(model = salsa1dOutput.multi$bestModel, 
                data = vpdata, 
                varlist.in = varlist,
                factorlist.in = "floodebb",
                showKnots = TRUE, 
                type='link', 
                includeB0 = TRUE)

## -----------------------------------------------------------------------------
initialModel <- glm(response ~ 1 + floodebb + offset(log(area)), family = "quasipoisson", 
                    data = vpdata)

varlist <- c("observationhour", "x.pos", "y.pos", "MonthOfYear")

salsa1dlist <- list(fitnessMeasure = "QBIC", 
                    minKnots_1d = rep(1, length(varlist)),
                    maxKnots_1d = rep(1, length(varlist)), 
                    startKnots_1d = rep(1, length(varlist)), 
                    degree = rep(2, length(varlist)),
                    gaps = rep(0, length(varlist)))

salsa1dOutput.multi.rm <- runSALSA1D(initialModel = initialModel, 
                            salsa1dlist = salsa1dlist,
                            varlist = varlist, 
                            datain = vpdata,
                            removal = TRUE, ##
                            suppress.printout = TRUE)

## -----------------------------------------------------------------------------
salsa1dOutput.multi.rm$modelFits[[5]]

## ----eval=TRUE, fig=TRUE, fig.align='center', fig.width=6, fig.height=4, message=FALSE----
runPartialPlots(model = salsa1dOutput.multi.rm$bestModel, 
                data = vpdata, 
                varlist.in = salsa1dOutput.multi.rm$keptvarlist,  ##
                factorlist.in = "floodebb",
                showKnots = TRUE, 
                type='link', 
                includeB0 = TRUE)

## -----------------------------------------------------------------------------
fit_rmfloodeb <- update(salsa1dOutput.multi.rm$bestModel, . ~ . - floodebb)

set.seed(123)
cv.gamMRSea(modelobject = salsa1dOutput.multi.rm$bestModel, 
            data = vpdata, 
            K = 10)$delta[2]
set.seed(123)
cv.gamMRSea(modelobject = fit_rmfloodeb, 
            data = vpdata, 
            K = 10)$delta[2]


## -----------------------------------------------------------------------------
anova(salsa1dOutput.multi.rm$bestModel)

