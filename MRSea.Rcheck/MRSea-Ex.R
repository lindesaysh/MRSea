pkgname <- "MRSea"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('MRSea')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("LocalRadialFunction")
### * LocalRadialFunction

flush(stderr()); flush(stdout())

### Name: LocalRadialFunction
### Title: Function for creating an exponential basis function for a
###   spatial smooth using the CReSS method.
### Aliases: LocalRadialFunction

### ** Examples

# load data
data(ns.data.re)
# load knot grid data
data(knotgrid.ns)

splineParams<-makesplineParams(data=ns.data.re, varlist=c('observationhour'))

#set some input info for SALSA
ns.data.re$response<- ns.data.re$birds

# make distance matrices for datatoknots and knottoknots
distMats<-makeDists(cbind(ns.data.re$x.pos, ns.data.re$y.pos), na.omit(knotgrid.ns), knotmat=FALSE)

# choose sequence of radii
r_seq<-getRadiiChoices(8, distMats$dataDist)

# using the fourth radius and picking 5 knots
basis<-LocalRadialFunction(radiusIndices=rep(4, 5), dists=distMats$dataDist, radii = r_seq,
        aR=c(3, 10, 15, 28, 31))



cleanEx()
nameEx("bootstrap.orig.data")
### * bootstrap.orig.data

flush(stderr()); flush(stdout())

### Name: bootstrap.orig.data
### Title: Obtaining a data frame of bootstrapped data using resamples
### Aliases: bootstrap.orig.data

### ** Examples

data(dis.data.re)
resample<-"transect.id"
samples<-unique(dis.data.re[,resample])
resamples.no<-length(samples)
new.resamples<-sample(samples,resamples.no,replace=TRUE)
bootstrap.data<-bootstrap.orig.data(dis.data.re,resample,new.resamples,resamples.no)



cleanEx()
nameEx("checkfactorlevelcounts")
### * checkfactorlevelcounts

flush(stderr()); flush(stdout())

### Name: checkfactorlevelcounts
### Title: Factor level response check
### Aliases: checkfactorlevelcounts

### ** Examples

# load data
data(ns.data.re)

checkfactorlevelcounts(factorlist=c('floodebb', 'impact'), ns.data.re,
     ns.data.re$birds)



cleanEx()
nameEx("create.NHAT")
### * create.NHAT

flush(stderr()); flush(stdout())

### Name: create.NHAT
### Title: Estimated number of individuals for each detection
### Aliases: create.NHAT

### ** Examples

data(dis.data.re)
result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re,method="ds",
     meta.data= list(width=250,binned=FALSE))
dis.data<-create.NHAT(dis.data.re,result)



cleanEx()
nameEx("create.bootcount.data")
### * create.bootcount.data

flush(stderr()); flush(stdout())

### Name: create.bootcount.data
### Title: Aggregate bootstrapped distance data into count data
### Aliases: create.bootcount.data

### ** Examples

data(dis.data.re)
# bootstrap data without stratification
dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
             meta.data=list(width=250))
dis.data.re<-create.NHAT(dis.data.re,result)

bootstrap.data<-create.bootstrap.data(dis.data.re)

bootcount.data<-create.bootcount.data(bootstrap.data)



cleanEx()
nameEx("create.bootstrap.data")
### * create.bootstrap.data

flush(stderr()); flush(stdout())

### Name: create.bootstrap.data
### Title: Create bootstrap data for non-parametric bootstrapping
### Aliases: create.bootstrap.data

### ** Examples

data(dis.data.re)
# run distance analysis to create NHATS
dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
             meta.data=list(width=250))
dis.data.re<-create.NHAT(dis.data.re,result)

# bootstrap data without stratification
bootstrap.data<-create.bootstrap.data(dis.data.re)
# boostrap data with stratification (here by survey which is composed of
# season and impact)
dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
bootstrap.data.str<-create.bootstrap.data(dis.data.re, stratum = "survey.id")



cleanEx()
nameEx("create.count.data")
### * create.count.data

flush(stderr()); flush(stdout())

### Name: create.count.data
### Title: Aggregate distance data into count data
### Aliases: create.count.data

### ** Examples

data(dis.data.re)
dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
           meta.data=list(width=250))
dis.data.re<-create.NHAT(dis.data.re,result)
count.data<-create.count.data(dis.data.re)



cleanEx()
nameEx("do.bootstrap.cress")
### * do.bootstrap.cress

flush(stderr()); flush(stdout())

### Name: do.bootstrap.cress
### Title: Bootstrapping function without model selection using CReSS/SALSA
###   for fitting the second stage count model
### Aliases: do.bootstrap.cress

### ** Examples

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# offshore redistribution data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data(dis.data.re)
data(predict.data.re)
data(knotgrid.off)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# distance sampling
dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
        meta.data=list(width=250))
dis.data.re<-create.NHAT(dis.data.re,result)
count.data<-create.count.data(dis.data.re)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# spatial modelling
splineParams<-makesplineParams(data=count.data, varlist=c('depth'))
#set some input info for SALSA
count.data$response<- count.data$NHAT
# make distance matrices for datatoknots and knottoknots
distMats<-makeDists(cbind(count.data$x.pos, count.data$y.pos), na.omit(knotgrid.off))
# choose sequence of radii
r_seq<-getRadiiChoices(8,distMats$dataDist)
# set initial model without the spatial term
initialModel<- glm(response ~ as.factor(season) + as.factor(impact) + offset(log(area)),
                family='quasipoisson', data=count.data)
# make parameter set for running salsa2d
salsa2dlist<-list(fitnessMeasure = 'QICb', knotgrid = knotgrid.off, knotdim=c(26,14), startKnots=4, minKnots=4,
                 maxKnots=20, r_seq=r_seq, gap=4000, interactionTerm="as.factor(impact)")
salsa2dOutput_k6<-runSALSA2D(initialModel, salsa2dlist, d2k=distMats$dataDist,
                   k2k=distMats$knotDist, splineParams=splineParams)

splineParams<-salsa2dOutput_k6$splineParams
# specify parameters for local radial function:
radiusIndices <- splineParams[[1]]$radiusIndices
dists <- splineParams[[1]]$dist
radii <- splineParams[[1]]$radii
aR <- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]
count.data$blockid<-paste(count.data$transect.id, count.data$season, count.data$impact, sep='')
# Re-fit the chosen model as a GEE (based on SALSA knot placement) and GEE p-values
geeModel<- geeglm(formula(salsa2dOutput_k6$bestModel), data=count.data, family=poisson, id=blockid)
dists<-makeDists(cbind(predict.data.re$x.pos, predict.data.re$y.pos), na.omit(knotgrid.off),
       knotmat=FALSE)$dataDist

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# bootstrapping
do.bootstrap.cress(dis.data.re, predict.data.re, ddf.obj=result, geeModel, splineParams,
              g2k=dists, resample='transect.id', rename='segment.id', stratum='survey.id',
              B=4, name="cress", save.data=FALSE, nhats=NULL, nCores=1)
load("cresspredictionboot.RData") # loading the bootstrap predictions into the workspace
# look at the first 6 lines of the bootstrap predictions (on the scale of the response)
head(bootPreds)

## Not run: 
##D # In parallel (Note: windows machines only)
##D require(parallel)
##D do.bootstrap.cress(dis.data.re, predict.data.re, ddf.obj=result, geeModel, splineParams,
##D                 g2k=dists, resample='transect.id', rename='segment.id', stratum='survey.id',
##D                 B=4, name="cress", save.data=FALSE, nhats=NULL, nCores=4)
##D load("cresspredictionboot.RData") # loading the bootstrap predictions into the workspace
##D # look at the first 6 lines of the bootstrap predictions (on the scale of the response)
##D head(bootPreds)
## End(Not run)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# nearshore redistribution data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Not run: 
##D do.bootstrap.cress(ns.data.re, ns.predict.data.re, ddf.obj=NULL, geeModel, splineParams,
##D              g2k=dists, resample='transect.id', rename='segment.id', stratum=NULL,
##D              B=2, name="cress", save.data=FALSE, nhats=NULL)
##D load("cresspredictionboot.RData") # loading the predictions into the workspace
##D # look at the first 6 lines of the bootstrap predictions (on the scale of the response)
##D head(bootPreds)
## End(Not run)



cleanEx()
nameEx("do.bootstrap.gam")
### * do.bootstrap.gam

flush(stderr()); flush(stdout())

### Name: do.bootstrap.gam
### Title: Bootstrapping function without model selection using 'gam' as
###   the second stage count model
### Aliases: do.bootstrap.gam

### ** Examples

# offshore redistribution data
data(dis.data.re)
data(predict.data.re)
dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
            meta.data=list(width=250))
dis.data.re<-create.NHAT(dis.data.re,result)
count.data<-create.count.data(dis.data.re)
require(mgcv)
gam.2<-gam(NHAT~as.factor(impact)+s(x.pos,y.pos,by=as.factor(impact))+offset(log(area)),
           data=count.data,family=quasipoisson)
do.bootstrap.gam(dis.data.re,predict.data.re,ddf.obj=result,model.obj=gam.2,resample="transect.id",
               rename="segment.id",stratum='survey.id',1,name='gam',save.data=FALSE,nhats=NULL)
load("gampredictionboot.RData") # loading the predictions into the workspace
# look at the first 6 lines of the predictions on the response scale
head(bootPreds)


## Not run: 
##D # nearshore redistribution data
##D data(ns.data.re)
##D data(ns.predict.data.re)
##D require(mgcv)
##D gam.ns2=gam(birds~as.factor(impact)+s(x.pos,y.pos,by=as.factor(impact))+offset(log(area)),
##D          data=ns.data.re,family=quasipoisson)
##D do.bootstrap.gam(ns.data.re,ns.predict.data.re,ddf.obj=NULL,model.obj=gam.ns2,resample=NULL,
##D                 rename=NULL,stratum=NULL,1,name='ns.gam',save.data=FALSE,nhats=NULL)
##D # load the replicate predictions into the workspace
##D load("ns.gampredictionboot.RData")
##D # look at the first 6 lines of the predictions on the response scale
##D head(bootPreds)
## End(Not run)



cleanEx()
nameEx("getCV_CReSS")
### * getCV_CReSS

flush(stderr()); flush(stdout())

### Name: getCV_CReSS
### Title: Calculate cross-validation score for a CReSS type model
### Aliases: getCV_CReSS

### ** Examples

# load data
data(ns.data.re)
# load prediction data
data(ns.predict.data.re)

splineParams<-makesplineParams(data=ns.data.re, varlist=c('observationhour', 'DayOfMonth'))
# set some input info for SALSA
ns.data.re$response<- ns.data.re$birds
salsa1dlist<-list(fitnessMeasure = 'QICb', minKnots_1d=c(2,2), maxKnots_1d = c(20, 20),
               startKnots_1d = c(2,2), degree=c(2,2), maxIterations = 10, gaps=c(1,1))

# set initial model without the spline terms in there
# (so all other non-spline terms)
initialModel<- glm(response ~ as.factor(floodebb) + as.factor(impact) + offset(log(area)),
                  family='quasipoisson',data=ns.data.re)

# run SALSA
salsa1dOutput<-runSALSA1D(initialModel, salsa1dlist, varlist=c('observationhour','DayOfMonth'),
            factorlist=c('floodebb', 'impact'), ns.predict.data.re, splineParams=splineParams)

# make blocking structure and fold structure
ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear,
                    ns.data.re$DayOfMonth, sep='')
ns.data.re$blockid<-as.factor(ns.data.re$blockid)
ns.data.re$foldid<-getCVids(ns.data.re, folds=5, block='blockid')

# calculate CV
cv1<-getCV_CReSS(ns.data.re, salsa1dOutput$bestModel, salsa1dOutput$splineParams)



cleanEx()
nameEx("getCVids")
### * getCVids

flush(stderr()); flush(stdout())

### Name: getCVids
### Title: IDs for running cross validation
### Aliases: getCVids

### ** Examples

# load data
data(ns.data.re)

CVids<-getCVids(ns.data.re, 5)



cleanEx()
nameEx("getDifferences")
### * getDifferences

flush(stderr()); flush(stdout())

### Name: getDifferences
### Title: Identify any significant differences between predicted data
###   before an impact event and predicted data after an impact event
### Aliases: getDifferences

### ** Examples

## Not run: 
##D getDifferences(beforePreds, afterPreds)
## End(Not run)



cleanEx()
nameEx("getPlotdimensions")
### * getPlotdimensions

flush(stderr()); flush(stdout())

### Name: getPlotdimensions
### Title: find the plotting dimensions for quilt.plot when using a regular
###   grid
### Aliases: getPlotdimensions

### ** Examples

#' # load data
data(ns.data.re)

getPlotdimensions(ns.data.re$x.pos, ns.data.re$y.pos, segmentWidth=500, segmentLength=500)



cleanEx()
nameEx("getPvalues")
### * getPvalues

flush(stderr()); flush(stdout())

### Name: getPvalues
### Title: Calculate marginal p-values from a 'model'.
### Aliases: getPvalues

### ** Examples

# load data
data(ns.data.re)

# make blocking structure
ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear,
                    ns.data.re$DayOfMonth, sep='')
ns.data.re$blockid<-as.factor(ns.data.re$blockid)

initialModel<- geeglm(birds ~ as.factor(floodebb) + as.factor(impact) + observationhour + x.pos +
              y.pos + offset(log(area)), family='poisson',data=ns.data.re, id=blockid)

getPvalues(initialModel, varlist=c('observationhour', 'x.pos', 'y.pos'),
            factorlist=c('floodebb', 'impact'))

getPvalues(initialModel)



cleanEx()
nameEx("getRadiiChoices")
### * getRadiiChoices

flush(stderr()); flush(stdout())

### Name: getRadiiChoices
### Title: Function for obtaining a sequence of range parameters for the
###   CReSS smoother
### Aliases: getRadiiChoices

### ** Examples

# load data
data(ns.data.re)
# load knot grid data
data(knotgrid.ns)

# make distance matrices for datatoknots and knottoknots
distMats<-makeDists(cbind(ns.data.re$x.pos, ns.data.re$y.pos), na.omit(knotgrid.ns))

# choose sequence of radii
r_seq<-getRadiiChoices(8, distMats$dataDist)



cleanEx()
nameEx("makeBootCIs")
### * makeBootCIs

flush(stderr()); flush(stdout())

### Name: makeBootCIs
### Title: Calculate percentile confidence intervals from a matrix of
###   bootstrapped predictions
### Aliases: makeBootCIs

### ** Examples

## Not run: 
##D makeBootCIs(bootPreds)
## End(Not run)



cleanEx()
nameEx("makeDists")
### * makeDists

flush(stderr()); flush(stdout())

### Name: makeDists
### Title: Make Euclidean distance matrices for use in CReSS and SALSA
###   model frameworks
### Aliases: makeDists

### ** Examples

# load data
data(ns.data.re)
# load knot grid data
data(knotgrid.ns)

# make distance matrices for datatoknots and knottoknots
distMats<-makeDists(cbind(ns.data.re$x.pos, ns.data.re$y.pos), na.omit(knotgrid.ns))



cleanEx()
nameEx("makesplineParams")
### * makesplineParams

flush(stderr()); flush(stdout())

### Name: makesplineParams
### Title: Constructing an object of spline parameters
### Aliases: makesplineParams

### ** Examples

# load data
data(ns.data.re)
# load prediction data
data(ns.predict.data.re)

splineParams<- makesplineParams(ns.data.re, varlist=c('observationhour', 'DayOfMonth'),
                predictionData=ns.predict.data.re)



cleanEx()
nameEx("plotCumRes")
### * plotCumRes

flush(stderr()); flush(stdout())

### Name: plotCumRes
### Title: Calculate cumulative residuals and plot.
### Aliases: plotCumRes

### ** Examples

# load data
data(ns.data.re)

model<-glm(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),
           family='quasipoisson', data=ns.data.re)

plotCumRes(model, varlist=c('observationhour'))



cleanEx()
nameEx("plotRunsProfile")
### * plotRunsProfile

flush(stderr()); flush(stdout())

### Name: plotRunsProfile
### Title: Calculate runs test and plot profile plot.  The output is a plot
###   of runs profiles (with p-value to indicate level of correlation)
### Aliases: plotRunsProfile

### ** Examples

# load data
data(ns.data.re)

model<-glm(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),
            family='quasipoisson', data=ns.data.re)

plotRunsProfile(model, varlist=c('observationhour'))



cleanEx()
nameEx("predict.cress")
### * predict.cress

flush(stderr()); flush(stdout())

### Name: predict.cress
### Title: Function for making predictions for a model containing a CReSS
###   basis (two dimensional local smooth).
### Aliases: predict.cress

### ** Examples

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# offshore redistribution data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data(dis.data.re)
data(predict.data.re)
data(knotgrid.off)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# distance sampling
dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
        meta.data=list(width=250))
dis.data.re<-create.NHAT(dis.data.re,result)
count.data<-create.count.data(dis.data.re)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# spatial modelling
splineParams<-makesplineParams(data=count.data, varlist=c('depth'))
#set some input info for SALSA
count.data$response<- count.data$NHAT
# make distance matrices for datatoknots and knottoknots
distMats<-makeDists(cbind(count.data$x.pos, count.data$y.pos), na.omit(knotgrid.off))
# choose sequence of radii
r_seq<-getRadiiChoices(8,distMats$dataDist)
# set initial model without the spatial term
initialModel<- glm(response ~ as.factor(season) + as.factor(impact) + offset(log(area)),
                family='quasipoisson', data=count.data)
# make parameter set for running salsa2d
salsa2dlist<-list(fitnessMeasure = 'QICb', knotgrid = knotgrid.off, knotdim=c(26,14), startKnots=4, minKnots=4,
                 maxKnots=20, r_seq=r_seq, gap=4000, interactionTerm="as.factor(impact)")
salsa2dOutput_k6<-runSALSA2D(initialModel, salsa2dlist, d2k=distMats$dataDist,
                   k2k=distMats$knotDist, splineParams=splineParams)

splineParams<-salsa2dOutput_k6$splineParams
# specify parameters for local radial function:
radiusIndices <- splineParams[[1]]$radiusIndices
dists <- splineParams[[1]]$dist
radii <- splineParams[[1]]$radii
aR <- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]
count.data$blockid<-paste(count.data$transect.id, count.data$season, count.data$impact, sep='')
# Re-fit the chosen model as a GEE (based on SALSA knot placement) and GEE p-values
geeModel<- geeglm(formula(salsa2dOutput_k6$bestModel), data=count.data, family=poisson, id=blockid)
dists<-makeDists(cbind(predict.data.re$x.pos, predict.data.re$y.pos), na.omit(knotgrid.off),
       knotmat=FALSE)$dataDist

# make predictions on response scale
preds<-predict.cress(predict.data.re, splineParams, dists, geeModel)



cleanEx()
nameEx("runACF")
### * runACF

flush(stderr()); flush(stdout())

### Name: runACF
### Title: run functions to create acf matrix and plot the results
### Aliases: runACF

### ** Examples

# load data
data(ns.data.re)

model<-glm(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),
           family='quasipoisson', data=ns.data.re)

ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear,
                    ns.data.re$DayOfMonth, sep='')
ns.data.re$blockid<-as.factor(ns.data.re$blockid)

runACF(ns.data.re$blockid, model)



cleanEx()
nameEx("runDiagnostics")
### * runDiagnostics

flush(stderr()); flush(stdout())

### Name: runDiagnostics
### Title: functions to create observed vs fitted and fitted vs scaled
###   pearsons residual plots
### Aliases: runDiagnostics

### ** Examples

# load data
data(ns.data.re)

model<-glm(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),
           family='quasipoisson', data=ns.data.re)

runDiagnostics(model)



cleanEx()
nameEx("runInfluence")
### * runInfluence

flush(stderr()); flush(stdout())

### Name: runInfluence
### Title: Assessing the influece of each correlated block on both the
###   precision of the parameter estimates (COVRATIO statistics) and the
###   sensitivity of model predictions (PRESS statistics).
### Aliases: runInfluence

### ** Examples

# load data
data(ns.data.re)

ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear,
                    ns.data.re$DayOfMonth, sep='')
ns.data.re$blockid<-as.factor(ns.data.re$blockid)

model<-geeglm(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),
              family='poisson', data=ns.data.re, id=blockid)

timeInfluenceCheck(model, ns.data.re$blockid)

## Not run: 
##D # **WARNING** this example takes a long time
##D influences<-runInfluence(model, ns.data.re$blockid)
## End(Not run)



cleanEx()
nameEx("runPartialPlots")
### * runPartialPlots

flush(stderr()); flush(stdout())

### Name: runPartialPlots
### Title: Plot partial plots for each of the variables listed in
###   'factorlist' or 'varlist'.
### Aliases: runPartialPlots

### ** Examples

#' # load data
data(ns.data.re)

model<-glm(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),
           family='quasipoisson', data=ns.data.re)

runPartialPlots(model, ns.data.re, factorlist=c('floodebb', 'impact'),
                varlist=c('observationhour'))



cleanEx()
nameEx("runSALSA1D")
### * runSALSA1D

flush(stderr()); flush(stdout())

### Name: runSALSA1D
### Title: Running SALSA for continuous one-dimensional covariates.
### Aliases: runSALSA1D

### ** Examples

# load data
data(ns.data.re)
# load prediction data
data(ns.predict.data.re)

splineParams<-makesplineParams(data=ns.data.re, varlist=c('observationhour', 'DayOfMonth'))
#set some input info for SALSA
ns.data.re$response<- ns.data.re$birds

#' # set initial model without the spline terms in there
# (so all other non-spline terms)
initialModel<- glm(response ~ as.factor(floodebb) + as.factor(impact) + offset(log(area)),
                    family='quasipoisson',data=ns.data.re)

salsa1dlist<-list(fitnessMeasure = 'QICb', minKnots_1d=c(2,2), maxKnots_1d = c(20, 20),
                  startKnots_1d = c(2,2), degree=c(2,2), maxIterations = 10, gaps=c(1,1))
# run SALSA
salsa1dOutput<-runSALSA1D(initialModel, salsa1dlist, varlist=c('observationhour', 'DayOfMonth'),
                 factorlist=c('floodebb', 'impact'), ns.predict.data.re, splineParams=splineParams)



cleanEx()
nameEx("runSALSA1D_withremoval")
### * runSALSA1D_withremoval

flush(stderr()); flush(stdout())

### Name: runSALSA1D_withremoval
### Title: Running SALSA for continuous one-dimensional covariates.
### Aliases: runSALSA1D_withremoval

### ** Examples

# load data
data(ns.data.re)
# load prediction data
data(ns.predict.data.re)

splineParams<-makesplineParams(data=ns.data.re, varlist=c('observationhour', 'DayOfMonth'))

# make column with foldid for cross validation calculation
ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear, ns.data.re$DayOfMonth, sep='')
ns.data.re$blockid<-as.factor(ns.data.re$blockid)
ns.data.re$foldid<-getCVids(ns.data.re, folds=5, block='blockid')

#' # set initial model without the spline terms in there
# (so all other non-spline terms)
ns.data.re$response<- ns.data.re$birds
initialModel<- glm(response ~ as.factor(floodebb) + as.factor(impact) + offset(log(area)),
                    family='quasipoisson',data=ns.data.re)

#set some input info for SALSA
salsa1dlist<-list(fitnessMeasure = 'QICb', minKnots_1d=c(2,2), maxKnots_1d = c(5, 5),
                  startKnots_1d = c(2,2), degree=c(2,2), maxIterations = 10, gaps=c(1,1))

# run SALSA
salsa1dOutput<-runSALSA1D_withremoval(initialModel, salsa1dlist, varlist=c('observationhour', 'DayOfMonth'),
                 factorlist=c('floodebb', 'impact'), ns.predict.data.re, splineParams=splineParams)



cleanEx()
nameEx("runSALSA2D")
### * runSALSA2D

flush(stderr()); flush(stdout())

### Name: runSALSA2D
### Title: Running SALSA for a spatial smooth with a CReSS basis
### Aliases: runSALSA2D

### ** Examples

# load data
data(ns.data.re)
# load prediction data
data(ns.predict.data.re)
# load knot grid data
data(knotgrid.ns)

splineParams<-makesplineParams(data=ns.data.re, varlist=c('observationhour'))

#set some input info for SALSA
ns.data.re$response<- ns.data.re$birds

# make distance matrices for datatoknots and knottoknots
distMats<-makeDists(cbind(ns.data.re$x.pos, ns.data.re$y.pos), na.omit(knotgrid.ns))

# choose sequence of radii
r_seq<-getRadiiChoices(8, distMats$dataDist)

# set initial model without the spatial term
# (so all other non-spline terms)
initialModel<- glm(response ~ as.factor(floodebb) + as.factor(impact) + offset(log(area)),
                   family='quasipoisson', data=ns.data.re)

# make parameter set for running salsa2d
salsa2dlist<-list(fitnessMeasure = 'QICb', knotgrid = knotgrid.ns, knotdim = c(7, 9),
                  startKnots=6, minKnots=4, maxKnots=20, r_seq=r_seq, gap=1,
                  interactionTerm="as.factor(impact)")

salsa2dOutput_k6<-runSALSA2D(initialModel, salsa2dlist, d2k=distMats$dataDist,
                            k2k=distMats$knotDist, splineParams=splineParams)



cleanEx()
nameEx("timeInfluenceCheck")
### * timeInfluenceCheck

flush(stderr()); flush(stdout())

### Name: timeInfluenceCheck
### Title: Timing check to see how long it will take to run 'runInfluence'.
### Aliases: timeInfluenceCheck

### ** Examples

# load data
data(ns.data.re)

ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear,
                     ns.data.re$DayOfMonth, sep='')
ns.data.re$blockid<-as.factor(ns.data.re$blockid)
model<-geeglm(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),
               family='poisson', data=ns.data.re, id=blockid)

timeInfluenceCheck(model, ns.data.re$blockid)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
