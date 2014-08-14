require(MRSea)

data(dis.data.re)
data(predict.data.re)

x<-seq(min(dis.data.re$x.pos), max(dis.data.re$x.pos), by=2000)
y<-seq(min(dis.data.re$y.pos), max(dis.data.re$y.pos), by=2000)
grid<-expand.grid(x.pos=x, y.pos=y)

plot(dis.data.re$x.pos, dis.data.re$y.pos, pch=20, cex=0.5, asp=1,xlab='x.pos', ylab='y.pos')
points(grid, pch=20, col='red')

# draw round boundary as no boundary polygon available
bnd<-locator() # i have drawn the coastline to the top for the boundary
polymap(bnd, add=T)

marker<-rep(0, nrow(grid))
id<-which(inout(grid, bnd)==TRUE)
# place a 1 in the marker vector for knots we do not want
marker[id]<-1
points(grid[marker==1,], pch=20, col='red')

knotgrid<-grid
knotgrid[marker==1,]<-c(NA,NA)

png('knot1.png', height=600, width=900)
plot(dis.data.re$x.pos, dis.data.re$y.pos, pch=20, cex=0.5, asp=1, xlab='x.pos', ylab='y.pos')
polymap(bnd, add=T)
points(grid, pch=20, col='red')
points(knotgrid, pch=20, col='grey')
dev.off()
# there are still knots in the open water area to the bottom of the survey
# we can use a Eulidean distance matrix of distances between knots and data to remove these.
distcheck<-makeDists(datacoords=cbind(dis.data.re$x.pos, dis.data.re$y.pos), knotcoords=na.omit(knotgrid))

naid<-which(is.na(knotgrid[,1]))

png('knot2.png', height=600, width=900)
plot(dis.data.re$x.pos, dis.data.re$y.pos, pch=20, cex=0.5, asp=1, xlab='x.pos', ylab='y.pos')
polymap(bnd, add=T)
points(grid, pch=20, col='red')
points(knotgrid, pch=20, col='lightblue')
points(knotgrid[-naid,][which(apply(distcheck$dataDist, 2, min)>1000),], pch=20, col='blue')
dev.off()

marker[-naid][which(apply(distcheck$dataDist, 2, min)>1000)]<-1
knotgrid[marker==1,]<-c(NA,NA)

png('knot3.png', height=600, width=900)
plot(dis.data.re$x.pos, dis.data.re$y.pos, pch=20, cex=0.5, asp=1, xlab='x.pos', ylab='y.pos')
polymap(bnd, add=T)
points(grid, pch=20, col='grey')
points(knotgrid, pch=20, col='purple')
dev.off()