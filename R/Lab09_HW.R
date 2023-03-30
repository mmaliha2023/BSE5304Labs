Sys.getenv('USER')
LabNo="/Lab09"
#
# What needs to be loaded
#
if (!require("pacman")) install.packages("pacman")
myhomedir=Sys.getenv("HOME")
datadir=paste0(myhomedir,"/data",LabNo)
dir.create(datadir,recursive = T)
srcdir=paste0(myhomedir,"/src")
dir.create(srcdir,recursive = T)
# Setting the directory for where the GitHub project exists. 
# This depends on where you set up your git, and what you called it locally, 
# but when you start a new git project, it will be the first directory you 
# are placed in... or if later in the project:
# 
mygitdir=rstudioapi::getActiveProject()
mypdfdir=paste0(mygitdir,"/pdfs",LabNo)
dir.create(mypdfdir)
# 
setwd(mygitdir)
system("git config --global user.email 'mushtari@vt.edu' ") 
system("git config --global user.name 'Mushtari Maliha' ")
system("git config pull.rebase false")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,elevatr,raster,rgdal,
               data.table,foreign,maptools,dataRetrieval,gdistance)
setwd(datadir)



##############################################
# 0205551460 LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA
##############################################
make_usgs_gage_list=function(siteNo = "0205551460",
                             parameterCd = c("00060","00065"),
                             start.date = "2017-05-01",
                             end.date = "2017-11-01")#grad students need to change the dates 
{
  USGS=list()   # Organize the data in a nice list as in previous labs
  USGS[["flowdata"]]<- readNWISuv(siteNumbers = siteNo,parameterCd = parameterCd,startDate = start.date,endDate = end.date)
  #head(USGS$flowdata)  # Note that we have 00060 and 00065...
  
  # And of course we want to work in SI units so:
  USGS$flowdata$depth_m=USGS$flowdata$X_00065_00000*0.3048
  # m/ft depth
  USGS$flowdata$cms=USGS$flowdata$X_00060_00000*.02832
  # m3/ft3 flow
  #
  # Let's add in the USGS gage site information to the list and inspect
  USGS[["site"]]=readNWISsite(siteNo)
  #head(USGS$site)
  #class(USGS$site$dec_lat_va)
  #
  # Set the Manning Coefficient in the USGS Gage's Site Table
  #
  USGS$site$man_n=.035/1.49
  #
  # Create a SpatialPointsDataFrame out of the site dataframe in the USGS list
  coordinates(USGS$site)=~dec_long_va+dec_lat_va
  
  return(USGS)
}

testlist=make_usgs_gage_list()
# If done correctly, your function should be able to populate
# the lists for the remaining gages!
#
USGS02056000=make_usgs_gage_list(siteNo = "02056000")
USGS0205551460=make_usgs_gage_list(siteNo ="0205551460" )
USGS02055100=make_usgs_gage_list(siteNo ="02055100" )
USGS02055000=make_usgs_gage_list(siteNo ="02055000" )
USGS02054530=make_usgs_gage_list(siteNo ="02054530" )

##########HW_1##########################
plot(USGS0205551460$flowdata$cms, USGS0205551460$flowdata$depth_m)
View(USGS0205551460$flowdata)
#Fixing the rating curve
original_USGS0205551460 <- USGS0205551460
View(original_USGS0205551460$flowdata)



USGS0205551460[["rating"]]=readNWISrating(USGS0205551460$site$site_no)
plot(USGS0205551460$rating$DEP,USGS0205551460$rating$INDEP,xlab="DEP",ylab="INDEP")
View(USGS0205551460$rating)
USGS0205551460$flowdata$X_00065_00000=approx(USGS0205551460$rating$DEP,
                                             USGS0205551460$rating$INDEP, xout = USGS0205551460$flowdata$X_00060_00000, ties = min)$y
points(USGS0205551460$flowdata$X_00060_00000,USGS0205551460$flowdata$X_00065_00000,
         col="red")

USGS0205551460$flowdata$depth_m=USGS0205551460$flowdata$X_00065_00000*0.3048
View(USGS0205551460$flowdata)

###################################HW_2##################################################

ab_ll=rbind(USGS02056000$site,
            USGS0205551460$site,
            USGS02055100$site,
            USGS02055000$site,
            USGS02054530$site)
class(ab_ll)
ab_ll@proj4string
proj4_utm = paste0("+proj=utm +zone=",
                   trunc((180+coordinates(USGS02056000$site)[1])/6+1), #changed USGS02055000 to USGS02056000
                   " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)
# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
proj4string(ab_ll)=proj4_ll
ab_utm=spTransform(ab_ll,crs_utm)
ab_utm@coords
mydem=get_aws_terrain(locations=ab_utm@coords, 
                      z = 12, prj = proj4_utm,expand=1)#grad students need to change z

plot(mydem)
plot(ab_utm,add=T)
text(ab_utm, labels=ab_utm@data$site_no, cex=0.6, font=2,pos=1)
#USGS H
url="https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/Shape/NHD_H_03010101_HU8_Shape.zip"
curl_download(url,"NHD_H_03010101_HU8_Shape.zip")
unzip("NHD_H_03010101_HU8_Shape.zip",exdir="03010101")
streams=readOGR("03010101/Shape/NHDFlowline.dbf")
streams_utm=spTransform(streams,crs_utm)
plot(streams_utm,col="blue",add=T)


A=SpatialPoints(USGS02055100$site)# Up gradient site at the top gage (changed, was 0205551460)
B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),1100)) #grad students be careful about setting the extent. make sure the dem contains the whole stream
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
plot(cropmydem)
plot(ab_utm,add=T)
text(ab_utm, labels=ab_utm@data$site_no, cex=0.6, font=2,pos=1)

#getting same data for USGS02054530
C=SpatialPoints(USGS02054530$site)# roanoke river at glenvar

proj4string(C)=proj4_ll

C_utm=spTransform(C,crs_utm)

#getting same data for USGS02055000
D=SpatialPoints(USGS02055000$site)# roanoke river at glenvar

proj4string(D)=proj4_ll

D_utm=spTransform(D,crs_utm)

##########################getting dimensions##################################################################
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)

# Find and plot the flow path and slope for USGS02055100
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
USGS02055100$site$L=SpatialLinesLengths(AtoB) # km to m
USGS02055100$site$L # reach length in m

USGS02055100$site$slope=(extract(mydem,A_utm)-extract(mydem,B_utm))/USGS02055100$site$L
USGS02055100$site$slope

# Find and plot the flow path and slope for USGS02054530
CtoB <- shortestPath(Conductance, C_utm, B_utm, output="SpatialLines")
plot(CtoB,add=T)
SpatialLinesLengths(CtoB)
USGS02054530$site$L=SpatialLinesLengths(CtoB) # km to m
USGS02054530$site$L # reach length in m

USGS02054530$site$slope=(extract(mydem,C_utm)-extract(mydem,B_utm))/USGS02054530$site$L
USGS02054530$site$slope

# Find and plot the flow path and slope for USGS02055000
DtoB <- shortestPath(Conductance, D_utm, B_utm, output="SpatialLines")
plot(DtoB,add=T)

SpatialLinesLengths(DtoB)
USGS02055000$site$L=SpatialLinesLengths(DtoB) # km to m
USGS02055000$site$L # reach length in m

USGS02055000$site$slope=(extract(mydem,C_utm)-extract(mydem,B_utm))/USGS02055000$site$L
USGS02055000$site$slope


# ck for USGS02055100
USGS02055100$flowdata$ck = 5/3*((sqrt(USGS02055100$site$slope))/USGS02055100$site$man_n)*USGS02055100$flowdata$depth_m^(2/3)
mean(USGS02055100$flowdata$ck,na.rm=T)

# ck for USGS02054530
USGS02054530$flowdata$ck = 5/3*((sqrt(USGS02054530$site$slope))/USGS02054530$site$man_n)*USGS02054530$flowdata$depth_m^(2/3)
mean(USGS02054530$flowdata$ck,na.rm=T)


# ck for USGS02055000
USGS02055000$flowdata$ck = 5/3*((sqrt(USGS02055000$site$slope))/USGS02055000$site$man_n)*USGS02055000$flowdata$depth_m^(2/3)
mean(USGS02055000$flowdata$ck,na.rm=T)

#dt for USGS02055100
USGS02055100$flowdata$dt = USGS02055100$site$L/USGS02055100$flowdata$ck
mean(USGS02055100$flowdata$dt,na.rm=T)
plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$dt)
USGS02055100$flowdata$outTime=USGS02055100$flowdata$dateTime+
  USGS02055100$flowdata$dt


#dt for USGS02054530
USGS02054530$flowdata$dt = USGS02054530$site$L/USGS02054530$flowdata$ck
mean(USGS02054530$flowdata$dt,na.rm=T)
plot(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$dt)
USGS02054530$flowdata$outTime=USGS02054530$flowdata$dateTime+
  USGS02054530$flowdata$dt

#dt for USGS02055000
USGS02055000$flowdata$dt = USGS02055000$site$L/USGS02055000$flowdata$ck
mean(USGS02055000$flowdata$dt,na.rm=T)
plot(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$dt)
USGS02055000$flowdata$outTime=USGS02055000$flowdata$dateTime+
  USGS02055000$flowdata$dt
##Top gauge
# Find the beginning of  Waves assuming a new wave starts at 110% of prior 
# flow. This might need to change for your homework
WaveStartDecPercent=2.5
USGS02055100$flowdata$newwave=
  USGS02055100$flowdata$cms *WaveStartDecPercent <
  data.table::shift(USGS02055100$flowdata$cms)
summary(USGS02055100$flowdata$newwave)
# Add plot of the point found
len=length(USGS02055100$flowdata$newwave)
USGS02055100$flowdata$newwave[is.na(USGS02055100$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(USGS02055100$flowdata$newwave[i]==T &
     USGS02055100$flowdata$newwave[i-1]==T){
    USGS02055100$flowdata$newwave[i]=F
  }
}
plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$cms,type="l")
points(USGS02055100$flowdata$dateTime[USGS02055100$flowdata$newwave],
       USGS02055100$flowdata$cms[USGS02055100$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(USGS02055100$flowdata$newwave == TRUE)
plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$cms,
     type="l",xlim=c(USGS02055100$flowdata$dateTime[1109],
                     USGS02055100$flowdata$dateTime[1109+200]))
lines(USGS02055100$flowdata$outTime,USGS02055100$flowdata$cms,col=2)
title(main='Time series plot for the top gauge')
View(USGS02055100$flowdata$outTime)
#bottom gauge
# Find the beginning of  Waves assuming a new wave starts at 110% of prior 
# flow. This might need to change for your homework
WaveStartDecPercent=2.5
USGS02055000$flowdata$newwave=
  USGS02055000$flowdata$cms *WaveStartDecPercent <
  data.table::shift(USGS02055000$flowdata$cms)
summary(USGS02055000$flowdata$newwave)
# Add plot of the point found
len=length(USGS02055000$flowdata$newwave)
USGS02055000$flowdata$newwave[is.na(USGS02055000$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(USGS02055000$flowdata$newwave[i]==T &
     USGS02055000$flowdata$newwave[i-1]==T){
    USGS02055000$flowdata$newwave[i]=F
  }
}
plot(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$cms,type="l")
points(USGS02055000$flowdata$dateTime[USGS02055000$flowdata$newwave],
       USGS02055000$flowdata$cms[USGS02055000$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(USGS02055000$flowdata$newwave == TRUE)
plot(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$cms,
     type="l",xlim=c(USGS02055000$flowdata$dateTime[1109],
                     USGS02055000$flowdata$dateTime[1109+200]))
lines(USGS02055000$flowdata$outTime,USGS02055000$flowdata$cms,col=2)
title(main='Time series plot for the bottom gauge')
View(USGS02055000$flowdata$outTime)


####################Grad student homework#########################################

upstream_USGS05247500=make_usgs_gage_list(siteNo ="05247500", start.date = "2020-05-01", #crow wing river near pillager, mn
                                          end.date = "2020-11-01" )
middle_USGS05267000=make_usgs_gage_list(siteNo ="05267000", start.date = "2020-05-01", #mississippi near royalton, mn
                                 end.date = "2020-11-01" )
downstream_USGS05270700=make_usgs_gage_list(siteNo ="05270700", start.date = "2020-05-01", #mississippi at st cloud, mn
                                        end.date = "2020-11-01" )
#correction for missing data at the upstream station
original_upstream_USGS05247500 <- upstream_USGS05247500
#View(original_upstream_USGS05247500$flowdata)
plot(original_upstream_USGS05247500$flowdata$X_00060_00000,original_upstream_USGS05247500$flowdata$X_00065_00000)
plot(upstream_USGS05247500$flowdata$X_00060_00000, upstream_USGS05247500$flowdata$X_00065_00000)

upstream_USGS05247500[["rating"]]=readNWISrating(upstream_USGS05247500$site$site_no)
plot(upstream_USGS05247500$rating$DEP,upstream_USGS05247500$rating$INDEP,xlab="DEP",ylab="INDEP")
#View(upstream_USGS05247500$rating)
upstream_USGS05247500$flowdata$X_00065_00000=approx(upstream_USGS05247500$rating$DEP,
                                             upstream_USGS05247500$rating$INDEP, xout = upstream_USGS05247500$flowdata$X_00060_00000, ties = min)$y
points(upstream_USGS05247500$flowdata$X_00060_00000,upstream_USGS05247500$flowdata$X_00065_00000,
       col="red")

upstream_USGS05247500$flowdata$depth_m=upstream_USGS05247500$flowdata$X_00065_00000*0.3048

#correction for missing data at the middle station
original_middle_USGS05267000 <- middle_USGS05267000
##View(original_middle_USGS05267000$flowdata)



middle_USGS05267000[["rating"]]=readNWISrating(middle_USGS05267000$site$site_no)
plot(middle_USGS05267000$rating$DEP,middle_USGS05267000$rating$INDEP,xlab="DEP",ylab="INDEP")
#View(middle_USGS05267000$rating)
middle_USGS05267000$flowdata$X_00065_00000=approx(middle_USGS05267000$rating$DEP,
                                             middle_USGS05267000$rating$INDEP, xout = middle_USGS05267000$flowdata$X_00060_00000, ties = min)$y
points(middle_USGS05267000$flowdata$X_00060_00000,middle_USGS05267000$flowdata$X_00065_00000,
       col="red")

middle_USGS05267000$flowdata$depth_m=middle_USGS05267000$flowdata$X_00065_00000*0.3048

#correction for missing data at the downstream station
original_downstream_USGS05270700 <- downstream_USGS05270700
#View(original_downstream_USGS05270700$flowdata)



downstream_USGS05270700[["rating"]]=readNWISrating(downstream_USGS05270700$site$site_no)
plot(downstream_USGS05270700$rating$DEP,downstream_USGS05270700$rating$INDEP,xlab="DEP",ylab="INDEP")
#View(downstream_USGS05270700$rating)
downstream_USGS05270700$flowdata$X_00065_00000=approx(downstream_USGS05270700$rating$DEP,
                                             downstream_USGS05270700$rating$INDEP, xout = downstream_USGS05270700$flowdata$X_00060_00000, ties = min)$y
points(downstream_USGS05270700$flowdata$X_00060_00000,downstream_USGS05270700$flowdata$X_00065_00000,
       col="red")

downstream_USGS05270700$flowdata$depth_m=downstream_USGS05270700$flowdata$X_00065_00000*0.3048





ab_ll=rbind(upstream_USGS05247500$site,
            middle_USGS05267000$site,
            downstream_USGS05270700$site)
class(ab_ll)
ab_ll@proj4string
proj4_utm = paste0("+proj=utm +zone=",
                   trunc((180+coordinates(downstream_USGS05270700$site)[1])/6+1), #changed 
                   " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)
# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
proj4string(ab_ll)=proj4_ll
ab_utm=spTransform(ab_ll,crs_utm)
ab_utm@coords
mydem=get_aws_terrain(locations=ab_utm@coords, 
                      z = 9, prj = proj4_utm,expand=1)#grad students need to change z

plot(mydem)
plot(ab_utm,add=T)
text(ab_utm, labels=ab_utm@data$site_no, cex=0.6, font=2,pos=1)
#USGS H
#url="https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU10/Shape/NHD_H_0701010409_HU10_Shape.zip"
#curl_download(url,"NHD_H_0701010409_HU10_Shape.zip")
#unzip("NHD_H_0701010409_HU10_Shape.zip",exdir="0701010409")
#streams=readOGR("0701010409/Shape/NHDFlowline.dbf")
#streams_utm=spTransform(streams,crs_utm)
#plot(streams_utm,col="blue",add=T)


A=SpatialPoints(upstream_USGS05247500$site)# Up gradient site at the top gage
B=SpatialPoints(downstream_USGS05270700$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),7000)) #grad students be careful about setting the extent. make sure the dem contains the whole stream
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
plot(cropmydem)
plot(ab_utm,add=T)
text(ab_utm, labels=ab_utm@data$site_no, cex=0.6, font=2,pos=1)

#getting same data for middle one
C=SpatialPoints(middle_USGS05267000$site)

proj4string(C)=proj4_ll

C_utm=spTransform(C,crs_utm)


##########################getting dimensions##################################################################
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)

# Find and plot the flow path and slope for upstream to downstream
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
upstream_USGS05247500$site$L=SpatialLinesLengths(AtoB) # km to m
upstream_USGS05247500$site$L # reach length in m

upstream_USGS05247500$site$slope=(extract(mydem,A_utm)-extract(mydem,B_utm))/upstream_USGS05247500$site$L
upstream_USGS05247500$site$slope

# Find and plot the flow path and slope for middle one to downstream
CtoB <- shortestPath(Conductance, C_utm, B_utm, output="SpatialLines")
plot(CtoB,add=T)
SpatialLinesLengths(CtoB)
middle_USGS05267000$site$L=SpatialLinesLengths(CtoB) # km to m
middle_USGS05267000$site$L # reach length in m

middle_USGS05267000$site$slope=(extract(mydem,C_utm)-extract(mydem,B_utm))/middle_USGS05267000$site$L
middle_USGS05267000$site$slope


# ck for upstream
upstream_USGS05247500$flowdata$ck = 5/3*((sqrt(upstream_USGS05247500$site$slope))/upstream_USGS05247500$site$man_n)*upstream_USGS05247500$flowdata$depth_m^(2/3)
mean(upstream_USGS05247500$flowdata$ck,na.rm=T)

# ck for middle one
middle_USGS05267000$flowdata$ck = 5/3*((sqrt(middle_USGS05267000$site$slope))/middle_USGS05267000$site$man_n)*middle_USGS05267000$flowdata$depth_m^(2/3)
mean(middle_USGS05267000$flowdata$ck,na.rm=T)


#dt for upstream
upstream_USGS05247500$flowdata$dt = upstream_USGS05247500$site$L/upstream_USGS05247500$flowdata$ck
mean(upstream_USGS05247500$flowdata$dt,na.rm=T)
plot(upstream_USGS05247500$flowdata$dateTime,upstream_USGS05247500$flowdata$dt)
upstream_USGS05247500$flowdata$outTime=upstream_USGS05247500$flowdata$dateTime+
  upstream_USGS05247500$flowdata$dt


#dt for middle_USGS05267000
middle_USGS05267000$flowdata$dt = middle_USGS05267000$site$L/middle_USGS05267000$flowdata$ck
mean(middle_USGS05267000$flowdata$dt,na.rm=T)
plot(middle_USGS05267000$flowdata$dateTime,middle_USGS05267000$flowdata$dt)
middle_USGS05267000$flowdata$outTime=middle_USGS05267000$flowdata$dateTime+
  middle_USGS05267000$flowdata$dt


##Top gauge
# Find the beginning of  Waves assuming a new wave starts at 110% of prior 
# flow. This might need to change for your homework
WaveStartDecPercent=1.1
upstream_USGS05247500$flowdata$newwave=
  upstream_USGS05247500$flowdata$cms *WaveStartDecPercent <
  data.table::shift(upstream_USGS05247500$flowdata$cms)
summary(upstream_USGS05247500$flowdata$newwave)
# Add plot of the point found
len=length(upstream_USGS05247500$flowdata$newwave)
upstream_USGS05247500$flowdata$newwave[is.na(upstream_USGS05247500$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(upstream_USGS05247500$flowdata$newwave[i]==T &
     upstream_USGS05247500$flowdata$newwave[i-1]==T){
    upstream_USGS05247500$flowdata$newwave[i]=F
  }
}
plot(upstream_USGS05247500$flowdata$dateTime,upstream_USGS05247500$flowdata$cms,type="l")
points(upstream_USGS05247500$flowdata$dateTime[upstream_USGS05247500$flowdata$newwave],
       upstream_USGS05247500$flowdata$cms[upstream_USGS05247500$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(upstream_USGS05247500$flowdata$newwave == TRUE)
plot(upstream_USGS05247500$flowdata$dateTime,upstream_USGS05247500$flowdata$cms,
     type="l",xlim=c(upstream_USGS05247500$flowdata$dateTime[1109],
                     upstream_USGS05247500$flowdata$dateTime[1109+200]))
lines(upstream_USGS05247500$flowdata$outTime,upstream_USGS05247500$flowdata$cms,col=2)
title(main='Time series plot for the top gauge')
View(upstream_USGS05247500$flowdata$outTime)
#middle gauge
# Find the beginning of  Waves assuming a new wave starts at 110% of prior 
# flow. This might need to change for your homework
WaveStartDecPercent=1.1
middle_USGS05267000$flowdata$newwave=
  middle_USGS05267000$flowdata$cms *WaveStartDecPercent <
  data.table::shift(middle_USGS05267000$flowdata$cms)
summary(middle_USGS05267000$flowdata$newwave)
# Add plot of the point found
len=length(middle_USGS05267000$flowdata$newwave)
middle_USGS05267000$flowdata$newwave[is.na(middle_USGS05267000$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(middle_USGS05267000$flowdata$newwave[i]==T &
     middle_USGS05267000$flowdata$newwave[i-1]==T){
    middle_USGS05267000$flowdata$newwave[i]=F
  }
}
plot(middle_USGS05267000$flowdata$dateTime,middle_USGS05267000$flowdata$cms,type="l")
points(middle_USGS05267000$flowdata$dateTime[middle_USGS05267000$flowdata$newwave],
       middle_USGS05267000$flowdata$cms[middle_USGS05267000$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(middle_USGS05267000$flowdata$newwave == TRUE)
plot(middle_USGS05267000$flowdata$dateTime,middle_USGS05267000$flowdata$cms,
     type="l",xlim=c(middle_USGS05267000$flowdata$dateTime[1109],
                     middle_USGS05267000$flowdata$dateTime[1109+200]))
lines(middle_USGS05267000$flowdata$outTime,middle_USGS05267000$flowdata$cms,col=2)
title(main='Time series plot for the bottom gauge')
View(middle_USGS05267000$flowdata$outTime)



