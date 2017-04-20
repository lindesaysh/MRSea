#' Function to calculate geodesic distances
#' 
#' @param xygrid Regular grid of x and y coordinates, must have an extent greater than that of \code{datalocations} if specified.
#' @param polys Either a single polygon (defined by a dataframe with x and y) or a list object containing multiple polygons that define exclusion zones.
#' @param datalocations data frame of x and y coordinates of locations to calculate distances between. Default = NULL, in which case distances are returned between locations in \code{xygrid}
#' 
#' @details Using the raster, spancs and gdistance packages to calculate the geodesic distance between pairs of points. The polygons define the areas that an animal/object cannot be. E.g. Land if the animal is a marine mammal. 
#' 
#' The object returned is a matrix of distances between pairs of points and the x,y coordinates of each pair of points.
#' 
#' @examples
#' 
#' Nx=100
#' Ny=100
#' 
#' xygrid<-expand.grid(x=1:Nx, y=1:Ny)
#' 
#' # make exclusion polygon,
#' bnd<-data.frame(x=c(78  ,79, 52, 53 , 82 , 82, 100, 103 ,104 , 70), 
#' y= c( 99  ,65 , 66 , 43 , 42 , 16, 19, 16 ,102 ,105))
#' bnd2<-data.frame(x= c(37 ,54 ,55, 37), y=c(35 ,35, 17, 18))
#' 
#' geodistsoutput<-getGeoDist(xygrid=xygrid, polys=list(bnd, bnd2))
#' 
#' # show on plot
#' i=2150
#' a.dist.mesh = data.frame(geodistsoutput$xydata, value=geodistsoutput$distance[i,])
#' iLCdsitance = rasterFromXYZ(a.dist.mesh)
#' plot(iLCdsitance, col=topo.colors(100))
#' points(geodistsoutput$xydata[i,1],geodistsoutput$xydata[i,2])
#' polymap(bnd, add=T)
#' polymap(bnd2, add=T)
#' 
#' @export
#' 

getGeoDist<-function(xygrid, polys, datalocations=NULL){
  
  require(raster)
  require(gdistance)
  
  if(class(polys)=='list'){
    inoutid<-rep(2, length=nrow(xygrid))
    for(b in 1:length(polys)){
      bnd<-polys[[b]]
      inoutid[which(splancs::inout(xygrid, bnd)==T)]<-0
    }
  }else{
    inoutid<-ifelse(splancs::inout(xygrid, polys)==T, 0, 2)
  }
  
  rastermesh = rasterFromXYZ(data.frame(xygrid, inoutid))
  #plot(rastermesh)
  rastermesh<-ratify(rastermesh)
  
  trans<-transition(rastermesh, "areas", 16)
  # nlayers(trans)
  # plot(raster(trans[[2]]))
  
  trans = geoCorrection(trans[[2]])
  
  if(is.null(datalocations)){
    matrixmesh = as.matrix(xygrid[inoutid==2,]) # create matrix of points to and from which we want least cost distances
  }else{
    matrixmesh = as.matrix(datalocations) # create matrix of points to and from which we want least cost distances
  }
  
  leastcostdistance = costDistance(trans, fromCoords=matrixmesh, toCoords=matrixmesh)
  
  return(list(xydata=matrixmesh, distance=leastcostdistance))
}
