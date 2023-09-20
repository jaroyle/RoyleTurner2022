
#
# Script 3: process the data into encounter data file (EDF) and trap deployment file (TDF), fit some SCR models using oSCR
#  
#
#
  

# Contains the effort covariate and other things output from Script 2
load("effortmatrix.RData")
# effortmatrix = ngridcell x ndays of survey effort.

# Total is the total across the whole season
total[total==0]<- NA  # do this for plotting
# Basic illustration of the effort covariate
par(mar=c(3,3,3,8))
load("utils.RData")
spatial.plot(gr,  (total), col=terrain.colors(20) )


# Figure for paper -- showing one day of effort (total area coverage by all searchers in each grid cell)
par(mar=c(3,3,3,8))
v<- effortmatrix[,6]; v[v==0]<- NA
r<- rasterFromXYZ(cbind(gr, v) )
image(r, col=terrain.colors(20) ) 
image.scale( seq(min(v,na.rm=TRUE),max(v,na.rm=TRUE),,10), col=terrain.colors(10) )
title("survey effort for occasion 6, 3 observers")

# You can look at each day's total effort
for(day in 1:60){
 r<- rasterFromXYZ(cbind(gr, effortmatrix[,day]) )
 image(r, col=terrain.colors(20) ) 
 image.scale( seq(-1.7e-07, 30,,20), col=terrain.colors(20) )
 title(paste("effort for day = ",day, sep=" ") )
  #browser()
}
 

# Load the 2020 encounter data 
load("capturedata.RData")

# Grab long and lat coords of each capture, convert to UTM
lat <- capdata$Latitude
long <- capdata$Longitude
capture <- data.frame(lat, long)
### points(capture$long, capture$lat, col = "black", pch =20)   # Long/lat coordinates here
xy <- (capdata[,c("Longitude","Latitude")])
utms6=project(as.matrix(xy), "+proj=utm +zone=18 ellps=WGS84")
points(utms6,pch=20,col="black")
xy<- data.frame(utms6)

# Figure out which grid cell each encounter is in (look at closest grid cell center) (grid cell == trap)
intrap<- rep(NA, nrow(xy))
for(i in 1:nrow(xy)){
  d<- sqrt( (xy[i,1] - gr[,1])^2 +  (xy[i,2] - gr[,2])^2  )
  intrap[i]<- (1:length(d))[d==min(d)]
}

#
# Figure 2 in the paper
#
#
#png("Fig2_new_v2.png", width=720,height=0.75*720, pointsize=14)
par(mar=c(5,5,5,8))
r<- rasterFromXYZ(cbind(gr, total) )
v<-values(r)
image(r, col=terrain.colors(20) ,xlab="UTM Easting",ylab="UTM Northing") 
image.scale( seq(min(v,na.rm=TRUE), max(v,na.rm=TRUE),,10), col=terrain.colors(10) )
##title("total grid cell effort (effective passes through cell)")
plot(buffered.tracks[[58]][[1]],add=TRUE)
plot(buffered.tracks[[5]][[1]],add=TRUE)
points(xy[,1:2],pch=20)
#dev.off()

 

# Now we're setting things up to run SCR models

##############
# make the trap deployment file TDF
##############
 
tdf<- data.frame(cbind(1:nrow(gr), gr))
names(tdf) <- c("De", "X", "Y")

# This needs to have 60 "trap operation" columns added to it. 0 if no effort for that grid cell... 1 if grid cell contains at least some effort
#   To do this we look at the entries of the effortmatrix:

trapopp<- effortmatrix
trapopp[trapopp>0]<- 1
tdf<- cbind(tdf,trapopp)


##############
#make the encounter data file EDF
##############

library(dplyr)  
 #make single column of turtle ID
 ID <- subset(capdata, select = c(Turtle.ID))
 day<- subset(capdata, select = c(day))

edf <- cbind(rep(1,length(ID)), intrap, ID, day)
colnames(edf)<- c("session","trap","turtleID","occasion")
 
edf

# Histogram of captures per day
 plot(table(edf[,4]),xlab="Day of season",ylab="Frequency of captures")
 

# Gotta install oSCR this direct from github it's not on CRAN:   https://github.com/jaroyle/oSCR
library("oSCR")

edf<- data.frame(edf)
edf[,4]<- as.numeric(edf[,4])
colnames(tdf)<- c("ID","X","Y",1:60)

#
# Make a covariate for each trap called "day" that runs 1:ndays and then standardize it. Also create a quadratic day covariate, and EFFORT (per grid cell)
# Then if you feed this tdf into data2oscr it will process the covariate appropritately. Extra columns of the TDF are recognized as covariates
#

day<- col(tdf[,4:ncol(tdf)])
day<- (day-30)/10
day2<- day^2
colnames(day)<- paste("day",1:ncol(day),sep=".")
colnames(day2)<- paste("day2", 1:ncol(day2), sep=".")
effort<- effortmatrix
colnames(effort)<- paste("effort",1:60,sep=".") 

# Add these to the TDF file
tdf<- cbind(tdf, sep = "/", day,day2, effort)
 


# Use data2oscr to process the edf and tdf.
# NOTE THAT NTRAPS SEEMS NOT TO BE RELEVANT -- I think this is obsolete in data2oscr will update data2oscr later
sf <-  data2oscr(edf, sess.col = 1, id.col = 3, occ.col = 4, trap.col = 2, sex.col = NULL, tdf = list(tdf), 
                  K = c(60), ntraps = 5076, remove.zeros = TRUE, 
                  trapcov.names=c("day","day2","effort"),remove.extracaps = TRUE, sex.nacode = NULL, tdf.sep = "/")

# Complex data structure
str(sf$scrFrame)

# Grab the scrFrame
scrFrame<- sf$scrFrame

# Shows the pattern of spatial recaptures. Uses a custom function.
plot(tdf[,2:3],  pch=".",axes=FALSE,xlab=" ",ylab=" ")  # You can see 1 unit has about 3 grid cells == 100 m

sp.out <- spiderplot(sf$scrFrame,session=1, add=TRUE)

 
# Create the state-space data frame
#  Do this automatically.  This creates a statespace object that is too fine (takes a long time to run)
#  This is too fine, should use a coarser grid
ssDF <- make.ssDF(sf$scrFrame, buffer=200, res = 30)

# Or just thin out the trap locations... here I take every other grid point and define this subset to be the statespace
# This produces 2390 state-space points
ssDF<- list(tdf[seq(1,nrow(tdf),2),2:3] )

#fit the NULL model: First run with getStarts = TRUE to see the order of the parameters, this helps you give the right starting values
#  (getStarts = TRUE) 
m0 <- oSCR.fit( list(D~1,p0~1,sig~1), 
      start.vals=c(0,0,0), scrFrame=sf$scrFrame, ssDF=ssDF , getStarts=TRUE)
print(m0)

# Then run it for real . Hear using starts that are close to the MLEs
m0 <- oSCR.fit( list(D~1,p0~1,sig~1), 
      start.vals=c(-5.7,4.02,-1.49), scrFrame=sf$scrFrame, ssDF=ssDF, encmod = "M", trimS=500)
# I save each model as I go
save(m0, file="m0.RData")

d0<- m0$coef.mle[3,2]  # estimated density, log scale
# pop size:  
 exp(d0)*nrow(ssDF[[1]])

# sigma: on the log scale , this is in meters b/c coordinate system is UTM
sig0<- m0$coef.mle[2,2]
  exp(sig0) 

# Model with linear date effect on detection probability. Again run with getStarts = True to see order of parameters
m1 <- oSCR.fit( list(D~1,p0~day,sig~1), 
      start.vals=c(-6,-0.54,-2.5), scrFrame=sf$scrFrame, ssDF=ssDF,trimS=10,getStarts=TRUE)

# output looks like this: (tells you order of start.vals)

$parameters
[1] "p0.(Intercept)"  "sig.(Intercept)" "t.beta.day"      "d0.(Intercept)" 

# Then run it for real
m1 <- oSCR.fit( list(D~1,p0~day,sig~1), 
      start.vals=c(-5.7, 4.0, 0, -1.5), scrFrame=sf$scrFrame, ssDF=ssDF, encmod = "M", trimS=500)

save(m1, file="m1.RData")
 
# Run model 2: effort, order of starts: p0, sig, t.beta.effort, d0
m2 <- oSCR.fit( list(D~1,p0~effort,sig~1), 
      start.vals=c(-5.76,-0.54, 0, -2.85), scrFrame=sf$scrFrame, ssDF=ssDF,trimS=10,  getStarts=TRUE)

# Run model 2: effort, order of starts: p0, sig, t.beta.effort, d0
m2 <- oSCR.fit( list(D~1,p0~effort,sig~1), 
      start.vals=c(-5.76,4.02, 0, -2.6), scrFrame=sf$scrFrame, ssDF=ssDF, encmod = "M" , trimS=500)

save(m2,file="m2.RData")
 
# Run model 3: effort, order of starts: p0, sig, t.beta.effort, d0
m3 <- oSCR.fit( list(D~1,p0~effort + day,sig~1), 
      start.vals=c(-5.76,-0.54, 0, -2.85), scrFrame=sf$scrFrame, ssDF=ssDF,trimS=10, getStarts=TRUE)

# Run model 2: effort, order of starts: p0, sig, t.beta.effort, d0
m3 <- oSCR.fit( list(D~1,p0~effort + day,sig~1), 
      start.vals=c(-5.76,4.42, 1.80, 0,-3.81), scrFrame=sf$scrFrame, ssDF=ssDF, encmod = "M", trimS=500)

save(m3,file="m3.RData")

# Run model 4: effort, order of starts: p0, sig, t.beta.effort, d0
m4 <- oSCR.fit( list(D~1,p0~effort + day + I(day^2),sig~1), 
      start.vals=c(-5.76,-0.54, 0, -2.85), scrFrame=sf$scrFrame, ssDF=ssDF,trimS=10, getStarts=TRUE)

# Run model 2: effort, order of starts: p0, sig, t.beta.effort, d0
m4 <- oSCR.fit( list(D~1,p0~effort + day + day2,sig~1), 
      start.vals=c(-5.76,4.42, 1.8, -.2,0.2, -3.81), scrFrame=sf$scrFrame, ssDF=ssDF, encmod = "M", trimS=500)

save(m4,file="m4.RData")

# Run model 2: effort, order of starts: p0, sig, t.beta.effort, d0
m4b <- oSCR.fit( list(D~1,p0~ day + day2,sig~1), 
      start.vals=c(-5.76,4.42, -.2,0.2, -3.81), scrFrame=sf$scrFrame, ssDF=ssDF, encmod = "M", trimS=500)
save(m4b, file="m4b.RData")

load("m0.RData"); load("m1.RData"); load("m2.RData"); load("m3.RData");
load("m4.RData"); load("m4b.RData")
  
fl <- fitList.oSCR(list(m0,m1,m2,m3, m4, m4b),rename=TRUE)
ms <- modSel.oSCR(fl)

# encounter probablitly at d(x,s)=0
# define the values we want to make prediction for
pred.df.det <- data.frame(effort = c(0.5,1,2,3,4,5,6,7,8,9,10), day = 0, day2=0)
#make predictions on the real scale
(pred.det <-get.real(m4,type ="det",newdata =pred.df.det))
(red.sig <- get.real(m4, type="sig",newdata=pred.df.det))
(pred.D <- get.real(m4, type="dens", newdata=pred.df.det, d.factor = 5.555556))  #1/0.18 = 5.55  , I think this is how many grid cells per ha

# Now per study area:
(pred.D <- get.real(m4, type="dens", newdata=pred.df.det, d.factor = (5.555556/2.2)*175  ))  #1/0.18 = 5.55.  175 -- ha?  not sure where this came from 


# Conversion to natural scale 
> (pred.D <- get.real(m4, type="dens", newdata=pred.df.det, d.factor = 5.555556))  #1/0.18 = 5.55 
  estimate        se      lwr      upr effort day day2
1 1.213916 0.1925106 0.889605 1.656457      1   0    0
 (pred.sig <- get.real(m4, type="sig",newdata=pred.df.det))
  estimate      se      lwr      upr effort day day2
1 63.31043 5.74344 52.99733 75.63043      1   0    0
  (pred.det <-get.real(m4,type ="det",newdata =pred.df.det))
  effort day day2    estimate           se         lwr         upr
1      1   0    0 0.003736079 0.0009287887 0.002294236 0.006078547
  pred.df.det <- data.frame(effort = c(0.5,1,2,3), day = 0, day2=0)
  (pred.det <-get.real(m4,type ="det",newdata =pred.df.det))

  effort day day2    estimate           se          lwr         upr
1    0.5   0    0 0.001480837 0.0003770978 0.0008988246 0.002438798
2    1.0   0    0 0.003736079 0.0009287887 0.0022942360 0.006078547
3    2.0   0    0 0.023417057 0.0066800697 0.0133457523 0.040774533
4    3.0   0    0 0.132939518 0.0441227355 0.0675164640 0.245094467
 

 

