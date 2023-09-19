
#
#
# Start SCR script right here
#
#

# need a lot of libraries to do the spatial analysis -- some of these are going obsolete I think.
library(rgdal)
library(rgeos)
library(sp)
library(raster)
library(chron)
library(spatialEco)
 
# Load the spatial data from the search tracks. Buffering the search tracks has to be done as an initial step (Step 1).
# Has buffered.tracks, 'gr', 'allpoints' and nfiles = vector of how many GPS tracks per search day (also probably not needed in the clean script)
load("spatial_data.RData")

# This function makes a grid cell into a polygon given the centroid and spacing
gridMaker <- function(anchor,grdSize,id, proj=NULL){
  anchor <- anchor-grdSize/2
  x <- c(anchor[1],anchor[1],anchor[1] + grdSize, anchor[1] + grdSize, anchor[1])
  y <- c(anchor[2],anchor[2] + grdSize, anchor[2] + grdSize, anchor[2], anchor[2])
  pp <- matrix(c(x,y),5,2,byrow=F)
  pp <- SpatialPolygons(list(
    Polygons(list(
      Polygon(pp)),id)))
  if(!is.null(proj)) projection(pp) <- proj
  return(pp)
  
}


# Make a plot of all search paths represented by points
plot(allpoints, pch=".", xlab = "Easting",ylab = "Northing")

# We define a rectangle that defines the state-space boundary for the SCR model
grid.boundary.utm <- list( x= c( -76.82817, -76.82823, -76.80073, -76.80067),
                           y = c( 39.04724, 39.06254, 39.06251, 39.04724) ) 
# Have to convert to UTM
utms5=project(as.matrix(cbind(grid.boundary.utm$x,grid.boundary.utm$y)), "+proj=utm +zone=18 ellps=WGS84")
grid.boundary.utm$x<- utms5[,1]
grid.boundary.utm$y<- utms5[,2]
grid.boundary<- grid.boundary.utm
xlim<- range(grid.boundary.utm$x)
ylim<- range(grid.boundary.utm$y)

# Make centroids of the grid which will define 'traps' for the SCR model
delta<- 10 # This is used a few times and represents 1/3 of the grid size. Should re-code things and have delta = 30
xgr<- seq(xlim[1],xlim[2],delta*3)
ygr<- seq(ylim[1],ylim[2],delta*3)
gr<-expand.grid(xgr,ygr)

plot(gr,cex=0.7,col="blue",pch=".")
points(allpoints,pch=".",col="red")

# Convert grid centers to a polygon converage
pg<- list()
for(i in 1:nrow(gr)){
  pg[[i]]<- gridMaker(as.matrix(gr[i,]), grdSize = 3*delta, id = as.character(i))
  plot(pg[[i]],add=TRUE)
 }

 

#
# This block of code intersects each search track with the trap polygons and computes how much area overlap there is
# uses this to create "effortmatrix" which I have saved and you can just load, so don't have to run this code
# Takes a long time!
# 

Xarea0<-array(0,c(nrow(gr), length(buffered.tracks), max(nfiles)))
for(days in 1:60){
cat("day: ", days, fill=TRUE)
    for(j in 1:nfiles[days]){
      cat("file: ", j, fill=TRUE)
      # for days = 1, j = 2, we have some polygon validity problems which I don't understand.... starting at i = 1015 I think
      for(i in 1:length(pg)){
       aa<- buffered.tracks[[days]][[j]]
       bb<- pg[[i]]
       zz<- gIntersection(aa,bb)
      if(is.null(zz) |  (class(zz) != "SpatialPolygons")  ){  # I think 2nd situation happens if there is no area overlap but one or more points are the same
        Xarea0[i,days,j] <- NA
      }else{
        Xarea0[i,days,j] <-  zz@polygons[[1]]@Polygons[[1]]@area    ### REPLACE with gArea()
      }
    }  # end j
  } # end i
 } # end days

Xarea<- Xarea0/( (3*delta)^2)   # expressed as a proportion of the grid cell
total<- apply(Xarea,c(1),sum,na.rm=TRUE)
effortmatrix<- apply(Xarea, c(1,2), sum,na.rm=TRUE)

### BE SURE TO SAVE, uncomment this line
#save(effortmatrix, Xarea, total, file="effortmatrix.RData")

# Make a plot to make sure things check out
 total[total==0]<- NA  # do this for plotting
# Basic illustration of the effort covariate
par(mar=c(3,3,3,8))
load("utils.RData")
spatial.plot(gr,  (total), col=terrain.colors(20) )







