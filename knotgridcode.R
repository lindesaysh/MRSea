# knotgrid choice
setwd('//BOTO/lindesay/WWT/WWT_partII')

require(splancs)
require(data.table)
# ~~~~~~~~~~~~~~~~~~~~~~

# read in data
data<-fread('Results//fullAnalysis//BH/BH-NHATS.csv')
head(data)


wwtpoly<-c(3:7, 17, 29, 31, 40:48, 50, 71:121)

NumberOfPolygons<-length(wwtpoly)
counter<-1
for(i in wwtpoly){
  # read in polys
  eval(parse(text=paste("p",counter, "<-read.table('Data/polygons/UTM30N/PhaseIIIStudyAreaProjected10kmIslandsSimplified_Polygon1_Part", i, ".txt')",sep='' ))) 
  counter<-counter+1
}

# ~~~~~~~~~~~~~~~~~~~~~~
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
bnd<-locator()
polymap(bnd, add=T)

marker<-rep(0, nrow(Grid))
id<-which(inout(Grid, bnd)==FALSE)
knotgrid<-Grid
marker[id]<-1

points(knotgrid[-id,], pch=20, col='green')

for(i in 1:NumberOfPolygons){
  id<-NULL
  eval(parse(text=paste("id<-which(inout(Grid, p",i,")==TRUE)", sep='')))
  if(is.null(id)==F){marker[id]<-1}
}
#id<-which(inout(grid, outerbnd)==FALSE)
#marker[id]<-1

#knotgrid<-Grid
knotgrid[marker==1,]<-c(NA,NA)

points(knotgrid, pch=20, col='blue')

# ~~~~~~~~~~~~~~~~~~~~~~
#find the grid point closest to each candidate knot on the grid
d<- c()
for(i in 1:nrow(data)){
  d[i]<- which.min(sqrt((data$x.pos[i]-knotgrid[,1])**2+(data$y.pos[i]-knotgrid[,2])**2))
}

RowsToKeep<- rep(0,nrow(knotgrid))
RowsToKeep[d]<-1

GridPosIncludingNAs<- cbind(ifelse(RowsToKeep!=1, NA, RowsToKeep*knotgrid[,1]), ifelse(RowsToKeep!=1, NA, RowsToKeep*knotgrid[,2]))
knotgrid <- as.data.frame(GridPosIncludingNAs)

points(knotgrid, pch=20)

dim(knotgrid)
dim(na.omit(knotgrid))

# ~~~~~~~~~~~~~~~~~~~~~~
# still too many knots so space fill those that are left to leave about 300 locations

require(fields)

naid<- which(is.na(knotgrid))
spaceid<-cover.design(R = na.omit(knotgrid), nd = 400, nruns = 5)$best.id

realid<-(1:nrow(knotgrid))[-naid]

points(knotgrid[realid[spaceid],], pch=20, col='purple')

knotgrid[realid[-spaceid],]<- c(NA, NA)
  
plot(data$x.pos, data$y.pos, pch=20, cex=0.2)
points(knotgrid, pch=20, cex=0.5, col='red')

write.csv(knotgrid, file='Data/knotgrid_fullanalysis.csv', row.names=F)


