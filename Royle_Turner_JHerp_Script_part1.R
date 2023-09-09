 

delta<- 10 ## for UTM calcs -- buffer width for GPS tracks

# Need a few libraries. I think some of these are becoming obsolete. 

library(rgdal)
library(rgeos)
library(sp)
library(raster)
library(chron)
library(spatialEco)
 


# It is assumed that you have the gpx files stored in directories BY DATE and that they are named
# the same as I have here. Obviously you can change this. "lod" = list of directories

 
# Grab the daily search directory names
 lod<- list.files("GPStracks")   
 lod<-lod[grep("2020",lod)]
 lod0<- lod
 # Split into year month and day components
 lod<- strsplit(lod,"[.]")
 lodmat<- matrix(NA,nrow=length(lod),ncol=3)
 for(i in 1:nrow(lodmat)){
   lodmat[i,1:3]<- as.numeric(lod[[i]])
 }
 lod0<- lod0[order(lodmat[,1],lodmat[,2])]
 lodmat<-lodmat[order(lodmat[,1],lodmat[,2]),]

 buffered.tracks<- line.tracks<- list()
 nfiles<- rep(0,length(lod0))

# We loop over search days and read in each track (here only up to 4 tracks). Then each track needs 
# converted to a uniform number of points PER UNIT TIME
# the objective here is to standardize each search track in terms of number of points per time searched
# and then use a buffered line coverage, intersect that with a grid and tally up how much 
# buffered line is in each grid cell , use that as a measure of search effort in the SCR model. 

# Number of search days (i.e., directories)
 ndays<- length(lod0)
 allpoints<- NULL

 for(i in 1:length(lod0)){
 
  fnames<- NULL
  lof<- list.files(paste("GPStracks/",lod0[i],sep=""))
  fnames<- lof[ grep("gpx",lof) ]
  nfiles[i]<- length(fnames)
  coords1<- coords2<- coords3<- coords4<- NULL
  regpoints1<- regpoints2<- regpoints3<- regpoints4<- NULL
  
  if(nfiles[i] >=1){
    f1 <- readOGR(paste("GPStracks/",lod0[i],"/",fnames[1],sep=""),"track_points")
    t1<-f1@data$time
    
    t1<- unlist(strsplit(t1,split=" "))[seq(2,2*length(t1),2)]
    t1<- substring(t1,1,8)
    thetimes<- chron(times=t1, format = c('h:m:s'))
    tot1<- 60*24*as.numeric(times(thetimes))
    tot1<- diff(range(tot1))
    
    coords1<- f1@coords
    line1 <- f1
    points1 <- SpatialPoints(line1)
    sLine1 <- Line(points1)
    
    regpoints1<- spsample(sLine1, 4*round(tot1,0), type="regular")  # 4 points per minute of sampling

    #   UTM conversion
    utms1=project(as.matrix(regpoints1@coords), "+proj=utm +zone=18 ellps=WGS84")
    regpoints1@coords<- utms1
    points1<- as(sp::SpatialPoints(regpoints1) ,"SpatialLines")

  }
  if(nfiles[i] >=2){
    f2 <- readOGR(paste("GPStracks/",lod0[i],"/",fnames[2],sep=""),"track_points")
    t2<-f2@data$time
    t2<- unlist(strsplit(t2,split=" "))[seq(2,2*length(t2),2)]
    t2<- substring(t2,1,8)
    thetimes<- chron(times=t2, format = c('h:m:s'))
    tot2<- 60*24*as.numeric(times(thetimes))
    tot2<- diff(range(tot2))
     # Sometimes a GPS track has some weird stuff that needs clipped out. Usually not important
     ##   if(i==33) f2@coords<- f2@coords[7:nrow(f2@coords),]

    coords2<- f2@coords
    line2 <- f2
    points2 <- SpatialPoints(line2)

    sLine2 <- Line(points2)
    
    plot(f2)

    regpoints2<- spsample(sLine2, 4*round(tot2,0), type="regular")  # 4 points per minute of sampling
    points(regpoints2)
    utms2=project(as.matrix(regpoints2@coords), "+proj=utm +zone=18 ellps=WGS84")
    regpoints2@coords<- utms2

    # Sausage making here, I had to clip a track for some reason
    # This is not important usually so I just commented it out here   
    if(i == 54){
      # regpoints2@coords <- regpoints2@coords[1:290,]
    }

    points2<- as(sp::SpatialPoints(regpoints2) ,"SpatialLines")

  }
 
  if(nfiles[i] >=3){
    f3 <- readOGR(paste("GPStracks/",lod0[i],"/",fnames[3],sep=""),"track_points")
    t3<-f3@data$time
    t3<- unlist(strsplit(t3,split=" "))[seq(2,2*length(t3),2)]
    t3<- substring(t3,1,8)
    thetimes<- chron(times=t3, format = c('h:m:s'))
    tot3<- 60*24*as.numeric(times(thetimes))
    tot3<- diff(range(tot3))
    
    coords3<- f3@coords
    line3 <- f3
 
    points3 <- SpatialPoints(line3)
    sLine3 <- Line(points3)
    
    regpoints3<- spsample(sLine3, 4*round(tot3,0), type="regular")  # 4 points per minute of sampling
    utms3=project(as.matrix(regpoints3@coords), "+proj=utm +zone=18 ellps=WGS84")
    regpoints3@coords<- utms3
    points3<- as(sp::SpatialPoints(regpoints3) ,"SpatialLines")
    
  }
    if(nfiles[i] >=4){
    f4 <- readOGR(paste("GPStracks/",lod0[i],"/",fnames[4],sep=""),"track_points")
    t4<-f4@data$time
    t4<- unlist(strsplit(t4,split=" "))[seq(2,2*length(t4),2)]
    t4<- substring(t4,1,8)
    thetimes<- chron(times=t4, format = c('h:m:s'))
    tot4<- 60*24*as.numeric(times(thetimes))
    tot4<- diff(range(tot4))
    
    coords4<- f4@coords
    line4 <- f4
    points4 <- SpatialPoints(line4)
    sLine4 <- Line(points4)
    
    regpoints4<- spsample(sLine4, 4*round(tot4,0), type="regular")  # 4 points per minute of sampling
    utms4=project(as.matrix(regpoints4@coords), "+proj=utm +zone=18 ellps=WGS84")
    regpoints4@coords<- utms4
    points4<- as(sp::SpatialPoints(regpoints4) ,"SpatialLines")
    
  }
 
  
  ## Keep the regularized tracks
  line.tracks[[i]]<- list(regpoints1=regpoints1, regpoints2=regpoints2,regpoints3=regpoints3, regpoints4=regpoints4)
  
  ## for each GPS file on a given day, do the buffering
  tmp<- NULL
  ##
    tmp[[1]]<- gBuffer(points1,width=delta)   
    allpoints<- rbind(allpoints, regpoints1@coords)  
    plot(tmp[[1]])
    if(nfiles[i]>=2){
      tmp[[2]]<- gBuffer(points2,width=delta)   
      allpoints<- rbind(allpoints, regpoints2@coords)
    }
    if(nfiles[i]>=3){
      tmp[[3]]<- gBuffer(points3,width=delta)  
      allpoints<- rbind(allpoints, regpoints3@coords)
    }
  if(nfiles[i]>=4){
      tmp[[4]]<- gBuffer(points4,width=delta)   
      allpoints<- rbind(allpoints, regpoints4@coords)
  }
 
  # Save the buffered tracks. Here buffered.tracks[[i]] are all the tracks for a given day. 
  # tmp is a LIST -- having all the tracks for each day
  buffered.tracks[[i]]<- tmp
  
}  # Ends loop over all directories and gpx files

save(buffered.tracks, nfiles, line.tracks, allpoints, file = "spatial_data.RData")

 

