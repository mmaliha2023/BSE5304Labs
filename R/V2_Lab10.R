if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,dplyr,patchwork,rnoaa)
pacman::p_load(operators,topmodel,DEoptim,soilDB,sp,curl,httr,
               rnoaa,raster,shapefiles,rgdal,elevatr,terra,progress,lubridate)
LabNo="/Lab10"
#
# Getting our organization on for where we want to put
# Data, external programs, and our project files.
# Things are going to get messy if we don't start issolating
# our data files by Lab
#
myhomedir=Sys.getenv("HOME")
datadir=paste0(myhomedir,"/data",LabNo)
dir.create(datadir,recursive = T)
srcdir=paste0(myhomedir,"/src")
dir.create(srcdir,recursive = T)

# Setting the directory for where the GitHub project exists. 
# This depends on where you set up your git, and what you called it locally, 
# but when you start a new git project, it will be the first directory you 
# are placed in... or if later in the project:
# WOOO HOOO... took me a few hours to find this function!
# 
mygitdir=rstudioapi::getActiveProject()
mypdfdir=paste0(mygitdir,"/pdfs",LabNo)
dir.create(mypdfdir)
# 
setwd(mygitdir)
system("git config --global user.email 'mushtari@vt.edu' ") 
system("git config --global user.name 'Mushtari Maliha' ")
system("git config pull.rebase false")
#
# This was already done before, and doesn't need to be repeated unless there
# is an update to R or the EcoHydRology Package... but 
#
setwd(srcdir)
system("svn checkout svn://scm.r-forge.r-project.org/svnroot/ecohydrology/"); 
install.packages(c("ecohydrology/pkg/EcoHydRology/"),repos = NULL)
pacman::p_load(EcoHydRology)

setwd(datadir)

myflowgage_id="0158397967" #Mine Bank Run at MD
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",
                         end_date = "2019-01-01")
#
# This is where some folks had issues... they forgot to check their 
# watershed areas per the homework... though there were ways to fix
# it later with lower resolution DEM pull
#
print(myflowgage$area)
# For most watershed modelling purposes we normalize Q in mm/day for basins
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

# In the Lab02, we introduced you to a way to quickly get your WX Data 
# for any location in the world way easier than traditional download and
# parsing methods most old people use.
#
WXData=FillMissWX(declat=myflowgage$declat, declon=myflowgage$declon,
                  StnRadius=30,minstns=10,date_min="2010-01-01",
                  date_max="2023-02-01",targElev=1,
                  method = "IDEW",alfa=2)

BasinData=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
#
# Setting the projection information for the specific location
#
proj4_utm = paste0("+proj=utm +zone=", trunc((180+myflowgage$declon)/6+1), " +datum=WGS84 +units=m +no_defs")

# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"

# Now we will build our proj4strings which define our “Coordinate 
# Reference Systems” or CRS in future geographic manipulations. 
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
#
# Double chec

latlon <- cbind(myflowgage$declon,myflowgage$declat)
myflowgage$gagepoint_ll <- SpatialPoints(latlon)
proj4string(myflowgage$gagepoint_ll)=proj4_ll
myflowgage$gagepoint_utm=spTransform(myflowgage$gagepoint_ll,crs_utm)
# Open up maps.google.com to guesstimate area/lengths
url=paste0("https://www.google.com/maps/@",
           myflowgage$declat,",",myflowgage$declon,",18z")
browseURL(url)
# We are going to over estimate our area
# For our search we are going to multiply the area by 6 and
# to get the distance
searchlength=sqrt(myflowgage$area*8)*1000 
pourpoint=SpatialPoints(myflowgage$gagepoint_utm@coords,proj4string = crs_utm)
bboxpts=myflowgage$gagepoint_utm@coords
bboxpts=rbind(bboxpts,bboxpts+searchlength)
bboxpts=rbind(bboxpts,bboxpts-searchlength)
bboxpts
bboxpts=rbind(bboxpts,c(min(bboxpts[,1]),max(bboxpts[,2])))
bboxpts=rbind(bboxpts,c(max(bboxpts[,1]),min(bboxpts[,2])))
bboxpts
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# From Lab04, get your DEM
mydem2=get_aws_terrain(locations=bboxpts@coords, 
                      z = 12, prj = proj4_utm,src ="aws",expand=1)
res(mydem2)
plot(mydem2)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")

# Write our raster to a geotiff file that can be used with
# OS level hydrological models 
writeRaster(mydem2,filename = "mydem2.tif",overwrite=T)
# Our quick intro to terminal where the cloud offerings are usually Linux
# ls; cd ~; pwd;  # Linux/Mac 
# dir; cd ; # Windows

#
# I am going to set two different zoom levels so I can inspect 
# the TauDEM Processing below.
#

zoomext=myflowgage$gagepoint_utm@coords
zoomext=rbind(zoomext,zoomext+res(mydem2)*100)
zoomext=rbind(zoomext,zoomext-res(mydem2)*100)
zoomext=SpatialPoints(zoomext,proj4string = crs_utm)  
zoomext2=myflowgage$gagepoint_utm@coords
zoomext2=rbind(zoomext2,zoomext2+res(mydem2)*10)
zoomext2=rbind(zoomext2,zoomext2-res(mydem2)*10)
zoomext2=SpatialPoints(zoomext2,proj4string = crs_utm)  
zoom(mydem2,ext=zoomext)
plot(pourpoint,add=T,col="red")

# If you already installed this in your ~/src directory and it 
# worked... you really 
# cd ~/src/      # Set your directory to your home directory
# git clone https://github.com/dtarb/TauDEM.git
# mkdir ~/src/TauDEM/bin
# cd ~/src/TauDEM/src
# sed -i -e 's/MPI_Type_struct/MPI_Type_create_struct/g' linklib.h
## yes, this next line is very small font, but it is one line so...
# sed -i -e 's/MPI_Type_extent(MPI_LONG, \&extent)/MPI_Aint lb\;MPI_Type_get_extent(MPI_LONG, \&lb, \&extent)/g' linklib.h
## Now let's try make again!
# make

rm("old_path")
old_path <- Sys.getenv("PATH")
old_path

if( ! grepl("~/src/TauDEM/bin",old_path)){
  Sys.setenv(PATH = paste(old_path,
                          paste0(Sys.getenv("HOME"),"/src/TauDEM/bin"), 
                          sep = ":"))
}

system("mpirun aread8")

setwd(datadir)
z=raster("mydem2.tif")
plot(z)

# Pitremove
system("mpiexec -n 2 pitremove -z mydem2.tif -fel mydem2fel.tif")
fel=raster("mydem2fel.tif")
plot(fel-z)


# D8 flow directions
system("mpiexec -n 2 d8flowdir -p mydem2p.tif -sd8 mydem2sd8.tif -fel mydem2fel.tif",show.output.on.console=F,invisible=F)
p=raster("mydem2p.tif")
plot(p)
sd8=raster("mydem2sd8.tif")
plot(sd8)

# Contributing area
system("mpiexec -n 2 aread8 -p mydem2p.tif -ad8 mydem2ad8.tif")
ad8=raster("mydem2ad8.tif")
plot(log(ad8))
zoom(log(ad8),ext=zoomext2)

plot(pourpoint,add=T)

# Grid Network 
system("mpiexec -n 2 gridnet -p mydem2p.tif -gord mydem2gord.tif -plen mydem2plen.tif -tlen mydem2tlen.tif")
gord=raster("mydem2gord.tif")
plot(gord)
zoom(gord,ext=zoomext2)

# DInf flow directions
system("mpiexec -n 2 dinfflowdir -ang mydem2ang.tif -slp mydem2slp.tif -fel mydem2fel.tif",show.output.on.console=F,invisible=F)
ang=raster("mydem2ang.tif")
plot(ang)
slp=raster("mydem2slp.tif")
plot(slp)

# Dinf contributing area
system("mpiexec -n 2 areadinf -ang mydem2ang.tif -sca mydem2sca.tif")
sca=raster("mydem2sca.tif")
plot(log(sca))
zoom(log(sca),ext=zoomext2)

targetbasins = 3
res(mydem2)
myflowgage$area
myflowgage$area * 10^3 * 10^3 / res(mydem2)[1]^2 / targetbasins
subthreshold=myflowgage$area * 10^3 * 10^3 / res(mydem2)[1]^2 / targetbasins
subthreshold=as.integer(subthreshold)

# Threshold
syscmd=paste0("mpiexec -n 2 threshold -ssa mydem2ad8.tif -src mydem2src.tif -thresh ", subthreshold)
system(syscmd)
src=raster("mydem2src.tif")
plot(src)
zoom(src,ext=zoomext2)
plot(pourpoint, add=T)
outlet=SpatialPointsDataFrame(myflowgage$gagepoint_utm,
                              data.frame(Id=c(1),outlet=paste("outlet",1,sep="")))
writeOGR(outlet,dsn=".",layer="approxoutlets",
         driver="ESRI Shapefile", overwrite_layer=TRUE)
#

# Move Outlets
system("mpiexec -n 2 moveoutletstostrm -p mydem2p.tif -src mydem2src.tif -o approxoutlets.shp -om outlet.shp")

approxpt=readOGR("approxoutlets.shp")
plot(approxpt,add=T, col="blue")
outpt=readOGR("outlet.shp")
plot(outpt,add=T, col="red")

# Contributing area upstream of outlet
# Now that we know the location of an outlet, we can isolate our basin 
#
system("mpiexec -n 2 aread8 -p mydem2p.tif -o outlet.shp -ad8 mydem2ssa.tif")
ssa=raster("mydem2ssa.tif")
plot(ssa) 

# Threshold
syscmd=paste0("mpiexec -n 2 threshold -ssa mydem2ssa.tif -src mydem2src1.tif -thresh ", subthreshold)
system(syscmd)
src1=raster("mydem2src1.tif")
plot(src1)
zoom(src1,ext=zoomext)

# Stream Reach and Watershed
system("mpiexec -n 2 streamnet -fel mydem2fel.tif -p mydem2p.tif -ad8 mydem2ad8.tif -src mydem2src1.tif -o outlet.shp -ord mydem2ord.tif -tree mydem2tree.txt -coord mydem2coord.txt -net mydem2net.shp -w mydem2w.tif")
plot(raster("mydem2ord.tif"))
zoom(raster("mydem2ord.tif"),ext=zoomext2)
mydem2w=raster("mydem2w.tif")
zoom(mydem2w,ext=zoomext)
summary(mydem2w)
plot(mydem2w)

# Trimming, Cropping, and Masking to make life prettier and easier
mydem2w=raster("mydem2w.tif")
mybasinmask=trim(mydem2w,padding=2)
mydem2=raster("mydem2.tif")
mybasindem=crop(mydem2,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
plot(mybasindem)

# Make a poly with raster library (slow)
# or from thee command line gdal (fast)
# gdal_polygonize.py -8 mydem2w.tif mydem2w_poly_gdal.shp
mydem2w=rast("mydem2w.tif")
mydem2w_poly=as.polygons(mydem2w,na.rm=T)
plot(mydem2w_poly,add=T,col=rainbow(6))

writeVector(mydem2w_poly,dsn=".",layer="mydem2w",driver="ESRI Shapefile", overwrite_layer=TRUE)
writeVector(mydem2w_poly, filename="mydem2w.shp", filetype="ESRI Shapefile", layer="mydem2w", insert=FALSE,
            overwrite=TRUE)

ssurgo.geom <- SDA_spatialQuery(
  mydem2w_poly,
  what = 'mupolygon',
  db = 'SSURGO',
  geomIntersection = TRUE
)

ssurgo.geom_utm=project(ssurgo.geom,crs_utm)
plot(ssurgo.geom_utm,col=rainbow(length(ssurgo.geom_utm)))
plot(mydem2w_poly,add=T)
ssurgo.geom_utm_crop=crop(ssurgo.geom_utm,mydem2w_poly)
#zip("mydem2w.zip",list.files(pattern="mydem2w[:.:]"))
plot(ssurgo.geom_utm_crop,col=rainbow(length(ssurgo.geom_utm)))
plot(mydem2w_poly,add=T)

# Get a list of the Map Units to use in accessing the soil data
#
unique(ssurgo.geom_utm_crop$mukey)
mukey_statement = format_SQL_in_statement(unique(ssurgo.geom_utm_crop$mukey))
print(mukey_statement)
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
print(q_mu2co)
mu2co = SDA_query(q_mu2co)
head(mu2co)
summary(mu2co)

# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
print(q_co2ch)
co2ch = SDA_query(q_co2ch)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
View(mu2ch)
summary(mu2ch)
mu2chmax=aggregate(mu2ch,list(mu2ch$mukey),max)
summary(mu2chmax)   	# What should we do with NAs?
# What do we have here vs our model for AWC?
soilmap=merge(ssurgo.geom_utm_crop,mu2chmax,all=T)

par(mfrow=c(2,2))
plot(soilmap,y="ksat_r")
plot(soilmap,y="awc_r")
plot(soilmap,y="hzdepb_r")
plot(mydem2w_poly,main="Subbasins",col=rainbow(length(mydem2w_poly)))

soilmap$ksat_r[is.na(soilmap$ksat_r)]=mean(soilmap$ksat_r,na.rm=T)
soilmap$awc_r[is.na(soilmap$awc_r)]=mean(soilmap$awc_r,na.rm=T)
soilmap$hzdepb_r[is.na(soilmap$hzdepb_r)]=mean(soilmap$hzdepb_r,na.rm=T)

par(mfrow=c(1,1))
# Lab 04 TI 
slp=raster("mydem2slp.tif")
plot(slp,ext=zoomext)
sca=raster("mydem2sca.tif")
plot(log(sca),e=zoomext)
TI = log( (sca+1)/(slp+0.00001) )
plot(TI)
zoom(log(TI),e=zoomext)

# To make things easier, while adding some confusion
# we will convert raster objects to "terra" supported SpatRasters
# and back again to use features that are unique to each
#
TI_terra=rast(TI)
TI_terra=crop(mask(TI_terra,mydem2w_poly),mydem2w_poly)
#

# Why would we want to mask the TI to the watershed
# boundaries before we build TI Classes?
# 
TI=raster(TI_terra)
pacman::p_load(classInt)

nTIclass=7 #number of TI classes, currently equal area, can adjust method various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI)
v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
#this brks.qt needs to be changed based on hydrology and objective.
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)
plot(TIC)
TIC_terra=rast(TIC)
plot(TIC_terra)

# Lab 04 Calibration
# Building more complex functions
# 

TMWB=BasinData
#
# Our model will
# 1) Calculate PET for the basin via Function
# 2) Calculate the Snow Accumulation and Melt via Function
# 3) Run TMWB via Function 
#
#
# First functions from last week we already have, Wetting, Drying, 
# and Wetting above capacity 

source("https://raw.githubusercontent.com/vtdrfuka/BSE5304Labs/main/R/TMWBFuncs.R")
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304Labs/a69125a7f3f133675928a5434024daef6c581e21/R/TISnow.R")

# Lets make one out of our Temperature Index Snow Model

SNO_df=TISnow(TMWB, SFTmp=2, bmlt6=1, bmlt12=0, Tmlt=3, Tlag=1)
TMWB$SNO=SNO_df$SNO
TMWB$SNOmlt=SNO_df$SNOmlt
TMWB$SNOfall=SNO_df$SNOfall
TMWB$Tsno=SNO_df$Tsno
#detach(TMWB)
#
# Our PET Model we will borrow from EcoHydrology
#

attach(TMWB)
TMWB$PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),Tmax_C = MaxTemp,Tmin_C = MinTemp,
                      lat_radians = myflowgage$declat*pi/180) * 1000
plot(date,TMWB$PET)
###############################################################FUNCTION################################
detach(TMWB)
source("https://raw.githubusercontent.com/mmaliha2023/BSE5304Labs/main/R/MyTMWBModel.R")
TMWBnew=TMWBmodel(TMWB, fcres=.2,FldCap=.20,WiltPt=.12,Z=200)
#Make a plot that has Qmm, P,and Qpred over time
plot(TMWBnew$date,TMWBnew$P,col="black")
plot(TMWBnew$date,TMWBnew$Qmm,type = "l",col="red", xlab= 'Year', ylab="Qmm", main="TMWB model obtained Qmm" )

plot(TMWBnew$date,TMWBnew$Qpred,type='l', col="blue",xlab= 'Year', ylab="Qpred", main="TMWB model obtained Qpred")
plot(TMWBnew$Qmm, TMWBnew$Qpred, xlab='Qmm obtained from TMWB model', ylab='Qpred obtained from TMWB model', main='x-y plot of Predicted vs Observed from TMWB')

NSE=function(Qobs,Qsim){
  return(1-sum((Qobs-Qsim)^2,na.rm=TRUE)/sum((Qobs-mean(Qobs, na.rm=TRUE))^2, na.rm=TRUE))
}
NSE(TMWBnew$Qmm, TMWBnew$Qpred)


#
# Functionalizing big big big time
# Here is a great place to make this into a function!
# return(TMWB)

BasinTMWB_JO=TMWBnew[(month(TMWBnew$date) > 5 
                      & month(TMWBnew$date) < 11),]
attach(BasinTMWB_JO)
plot(dP,Qmm)
detach(BasinTMWB_JO)

(1000/85-10)*25.4   # our CN estimate in bold
#[1] Storage=44.82353
(1000/50-10)*25.4   # our CN estimate in bold
#[1] Storage=254
#
# So we are going to visually "guestimate" that S should be somewhere between 
# 45mm and 260mm… repeat plotting until your solution covers the 
# largest Qmm vs dP event (upper right hand corner of plot). 
# 

# Assuming that (P-Ia) ~ dP, we can visually compare 
attach(BasinTMWB_JO)
plot(dP,Qmm)
points(dP,dP^2/(dP+45),col="red")  # S guestimates in bold
points(dP,dP^2/(dP+260),col="blue")# S guestimates in bold

# Vary S to maximize NSE using Eq. 4 of Lyon 2004 as our predictor of Q
#   Qpred=dP^2/(dP+S)
#
NSE(Qmm,dP^2/(dP+260))

NSE(Qmm,dP^2/(dP+45))

# Keep iterating until NSE is as high as you can get for your 

# best estimate to S (Sest)
#
f <- function (x) {
  Sest=x
  NSE(Qmm,dP^2/(dP+Sest))
}
y <- optimize(f, c(50,500), tol = 0.0001,maximum = TRUE)$maximum
Sest=f(y)
plot(dP,Qmm)
points(dP,dP^2/(dP+Sest),col="red") 
########
detach(BasinTMWB_JO)

mysoil_utm <- as(ssurgo.geom_utm_crop, "Spatial")
# Rasterizing for categorical analysis! 
rmysoil_utm=rasterize(mysoil_utm,TIC,field=as.numeric(mysoil_utm$mukey))
unique(rmysoil_utm)
pacman::p_load(circlize)
plot(rmysoil_utm,col=rand_color(length(unique(values(rmysoil_utm)))))
unique(rmysoil_utm)

ratify( TIC*10^9)

mybasinslp=mask(crop(slp,rmysoil_utm),rmysoil_utm)
unique(ratify(round(mybasinslp*10+1)))
#[1] 1 2 3 4 5 6 7
unique(ratify((rmysoil_utm*10^3)))
#[1] 520231000 520232000 520233000 520239000 520244000 520245000 520251000 520255000

#[25] 520386000 520414000 520432000
#
# Now build an HRU table with the combination of the 1) raster Soils, 2) TIC,
# and 3) slope layers. 
#
hru=ratify(TIC*10^9 + (rmysoil_utm*10^3) + round(mybasinslp*10+1))
unique(values(hru))
length(unique(values(hru)))
#[1] 485
plot(hru,col=rand_color(length(unique(values(hru)))))
# Think of how you will color this plot based on the sediment runoff you will
# calculate later.
#
# Build an HRU attribute table
hru_table = levels(hru)[[1]]
origID = hru_table$ID # preserve data order for later reordering
# metadata parameters from a string... this will make more sense
# after the next "head()" command
hru_table$TIclass = as.numeric(substr(sprintf("%10.0f", hru_table$ID), 1,1))
hru_table$mukey = as.numeric(substr(sprintf("%10.0f", hru_table$ID), 2,7))
hru_table$slp = (as.numeric(substr(sprintf("%10.0f", 
                                             hru_table$ID), 8,10))-1)/10
#
# Calculate the area for each unique soil (mukey) X TIClass combination
# using res(raster) for x and y resolution in units of m
# Note that table() function returns the count of the occurrences of
# unique values in the hru raster cells.
hru_table$areaSQKM = as.vector(round(res(hru)[1]*res(hru)[2]*
                                         table(values(hru))/10^6, 3))
#
# To better understand what happened, look at the new hru_table
head(hru_table)
summary(hru_table)
################### Not absolutely necessary #################
rm("mu2co")
for (mymukey in unique(hru_table$mukey)){
  print(mymukey)
  mukey_statement = format_SQL_in_statement(mymukey)
  q_mu2co = paste("SELECT mukey,cokey FROM component 
           WHERE mukey IN ", mukey_statement, sep="")
  if(!exists("mu2co")){
    mu2co=SDA_query(q_mu2co)} 
  else{
    mu2co=rbind(mu2co,SDA_query(q_mu2co))
  } 
}
View(mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
# cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
# q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r,frag3to10_r  
# FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
# co2ch = SDA_query(q_co2ch)
rm("co2ch")
for (mycokey in unique(mu2co$cokey)){
  print(mycokey)
  cokey_statement = format_SQL_in_statement(mycokey)
  q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r,frag3to10_r FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
  print(q_co2ch)
  if(!exists("co2ch")){
    co2ch=SDA_query(q_co2ch)
  } else{
    try((co2ch=rbind(co2ch,SDA_query(q_co2ch))))
  } 
}
View(co2ch)
rm("co2co")
for (mycokey in unique(mu2co$cokey)){
  print(mycokey)
  cokey_statement = format_SQL_in_statement(mycokey)
  q_co2co = paste("SELECT cokey,slopelenusle_r FROM component WHERE cokey IN ", cokey_statement, sep="")
  print(q_co2co)
  if(!exists("co2co")){
    co2co=SDA_query(q_co2co)} 
  else{
    try((co2co=rbind(co2co,SDA_query(q_co2co))))} 
}
View(co2co)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
mu2ch=merge(mu2ch,co2co)
View(mu2ch)

#
# Complete the table for running the MUSLE (Modified
# Universal Soil Loss Equation) to determine daily sediment
# loss, SWAT Theory eq. 4:1.1.1 . Assume:
# Kusle=.28
# Cusle=.2
# Pusle=.5
#
# Merge and then aggregate spatially derived hru_table with
# SSURGO summary table
 MUSLE_mrg=merge(hru_table,mu2ch)   
 MUSLE_mrg$ksat_r=as.numeric(MUSLE_mrg$ksat_r)
 MUSLE_mrg$awc_r=as.numeric(MUSLE_mrg$awc_r)
 MUSLE_mrg$hzdepb_r=as.numeric(MUSLE_mrg$hzdepb_r)
 MUSLE_mrg$slopelenusle_r=as.numeric(MUSLE_mrg$slopelenusle_r)
 MUSLE_mrg$frag3to10_r=as.numeric(MUSLE_mrg$frag3to10_r)
 MUSLE=aggregate(MUSLE_mrg,list(MUSLE_mrg$TIclass),mean,na.rm=T)
#
# Easiest first! Eq. 4:1.1.15 Course Fragment Factor
 MUSLE$CFRG=exp(-0.053*MUSLE$frag3to10_r)
 MUSLE
#
# LSusle is calculated using eq. 4.1.12
 MUSLE$alpha=atan(MUSLE$slp/100)
 MUSLE$LSm=.6*(1-exp(-35.835*MUSLE$slp/100))
 MUSLE$LS=(MUSLE$slopelenusle_r/22.1)^MUSLE$LSm * (65.41*sin(MUSLE$alpha)^2+4.56*sin(MUSLE$alpha)+0.065)
#
# Pusle
 MUSLE$Pusle=.50
#
# Cusle
 MUSLE$Cusle=.20
#
# Kusle
 MUSLE$Kusle=0.28
#
# Build a constant for those we are not changing day to day
 attach(MUSLE)
 MUSLE$KCPLSCFRG118=11.8*Kusle*Cusle*Pusle*LS*CFRG
 detach(MUSLE)
 MUSLE # Make sure values look correct, Pusle, Cusle, Kusle
#
# Now we need to use each of the TIClass Q solutions from Lab06 to calculate
# peak flows (qpeak) and complete the MUSLE Sediment Loss for each class.
# Run Model
#
# Now we need to use Q solutions from Lab06 to calculate
# peak flows (qpeak) and complete the MUSLE Sediment Loss for each class.
# Run Model

 #source CNmodel function
  source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/CNmodel")
  pacman::p_load(data.table)
 # We will split into 5 VSA areas represented by 5 TI Classes
  nTIclass=7
  VSAsol=data.table(TIClass=seq(from=nTIclass,to=3),
                     As=seq(1:nTIclass)*(1/nTIclass),Wetfrac=(1/nTIclass))
  VSAsol[,sSratio:=2*(sqrt(1-shift(As))-sqrt(1-As))/Wetfrac-1]
 #
  VSAsol$sSratio[1]=2*(sqrt(1-0)-sqrt(1-VSAsol$As[1]))/VSAsol$Wetfrac[1]-1
 # Calculate TI Class localized sigma and Curve Number
  VSAsol[,sigma:=Sest*sSratio]
  VSAsol[,CN:=25400/(sigma+254)]
  VSAsol
  VSAsol1 <- VSAsol[1:5,]
 
  VSAParams=merge(VSAsol1,MUSLE,by.x="TIClass",by.y="TIclass")
  View(VSAParams)
  ################################### Plots with original Pusle=0.5 ############################################# 
  TIC01=TMWB
  TIC02=TMWB
  TIC03=TMWB
  TIC04=TMWB
  TIC05=TMWB
 # For TIC01 CNavg=VSAParams$CN[1] but confirm
  TIC01 = CNmodel(CNmodeldf = TIC01, CNavg=VSAParams$CN[1], 
                   declat=myflowgage$declat,declon=myflowgage$declon)
  TIC01$qpeak=TIC01$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec
 
 TIC01$sed=(TIC01$Qpred*TIC01$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[1]    # Eq. 4:1.1.1 SWAT Theory

View(TIC01)
 glimpse(TIC01)
 
plot(TIC01$date, TIC01$sed, main = "", xlab="Date", ylab="Sediment Loss",
     col.axis="blue") 
title(main = "Time series of sediment loss for TI Class 1", sub = "Original model with Pusle = 0.5",
      cex.main = 2,   font.main= 2, col.main= "black",
      cex.sub = 1, font.sub = 3, col.sub = "blue",
      col.lab ="darkblue"
) 


# For TIC02 CNavg=VSAParams$CN[2] but confirm
TIC02 = CNmodel(CNmodeldf = TIC02, CNavg=VSAParams$CN[2], 
                declat=myflowgage$declat,declon=myflowgage$declon)
TIC02$qpeak=TIC02$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec

TIC02$sed=(TIC02$Qpred*TIC02$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[2]    # Eq. 4:1.1.1 SWAT Theory

glimpse(TIC02)

plot(TIC02$date, TIC02$sed, main = "", xlab="Date", ylab="Sediment Loss",
     col.axis="blue") 
title(main = "Time series of sediment loss for TI Class 2", sub = "Original model with Pusle = 0.5",
      cex.main = 2,   font.main= 2, col.main= "black",
      cex.sub = 1, font.sub = 3, col.sub = "blue",
      col.lab ="darkblue"
) 

# For TIC03 CNavg=VSAParams$CN[3] but confirm
TIC03 = CNmodel(CNmodeldf = TIC03, CNavg=VSAParams$CN[3], 
                declat=myflowgage$declat,declon=myflowgage$declon)
TIC03$qpeak=TIC03$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec

TIC03$sed=(TIC03$Qpred*TIC03$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[3]    # Eq. 4:1.1.1 SWAT Theory

glimpse(TIC03)

plot(TIC03$date, TIC03$sed, main = "", xlab="Date", ylab="Sediment Loss",
     col.axis="blue") 
title(main = "Time series of sediment loss for TI Class 3", sub = "Original model with Pusle = 0.5",
      cex.main = 2,   font.main= 2, col.main= "black",
      cex.sub = 1, font.sub = 3, col.sub = "blue",
      col.lab ="darkblue"
) 
# For TIC04 CNavg=VSAParams$CN[4] but confirm
TIC04 = CNmodel(CNmodeldf = TIC04, CNavg=VSAParams$CN[4], 
                declat=myflowgage$declat,declon=myflowgage$declon)
TIC04$qpeak=TIC04$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec

TIC04$sed=(TIC04$Qpred*TIC04$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[4]    # Eq. 4:1.1.1 SWAT Theory

glimpse(TIC04)
dev.off()
plot(TIC04$date, TIC04$sed, main = "", xlab="Date", ylab="Sediment Loss",
     col.axis="blue") 
title(main = "Time series of sediment loss for TI Class 4", sub = "Original model with Pusle = 0.5",
      cex.main = 2,   font.main= 2, col.main= "black",
      cex.sub = 1, font.sub = 3, col.sub = "blue",
      col.lab ="darkblue"
) 

# For TIC05 CNavg=VSAParams$CN[5] but confirm
TIC05 = CNmodel(CNmodeldf = TIC05, CNavg=VSAParams$CN[5], 
                declat=myflowgage$declat,declon=myflowgage$declon)
TIC05$qpeak=TIC05$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec

TIC05$sed=(TIC05$Qpred*TIC05$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[5]    # Eq. 4:1.1.1 SWAT Theory

glimpse(TIC05)

plot(TIC05$date, TIC05$sed, main = "", xlab="Date", ylab="Sediment Loss",
     col.axis="blue") 
title(main = "Time series of sediment loss for TI Class 5", sub = "Original model with Pusle = 0.5",
      cex.main = 2,   font.main= 2, col.main= "black",
      cex.sub = 1, font.sub = 3, col.sub = "blue",
      col.lab ="darkblue"
) 

TotAgg_origModel <- TIC01$sed + TIC02$sed + TIC03$sed + TIC04$sed + TIC05$sed
glimpse(TotAgg_origModel)
plot(TIC05$date, TotAgg_origModel, main = "", xlab="Date", ylab="Total Aggregate",
     col.axis="blue") 
title(main = "Time series of Total Aggregate", sub = "Original model with Pusle=0.50",
      cex.main = 2,   font.main= 2, col.main= "black",
      cex.sub = 1, font.sub = 3, col.sub = "blue",
      col.lab ="darkblue"
) 
################################ Creating spatial map for the original model #############################

TIC01$Year <- format(TIC01$date, format="%Y")

meanannualsedloss_orig_TIC01 <- TIC01 %>% 
  group_by(Year) %>% 
  summarize(mean(sed)) %>% 
  rename("annualsed_TIC01" = "mean(sed)")

TIC02$Year <- format(TIC02$date, format="%Y")

meanannualsedloss_orig_TIC02 <- TIC02 %>% 
  group_by(Year) %>% 
  summarize(mean(sed)) %>% 
  rename("annualsed_TIC02" = "mean(sed)")

TIC03$Year <- format(TIC03$date, format="%Y")

meanannualsedloss_orig_TIC03 <- TIC03 %>% 
  group_by(Year) %>% 
  summarize(mean(sed)) %>% 
  rename("annualsed_TIC03" = "mean(sed)")

TIC04$Year <- format(TIC04$date, format="%Y")

meanannualsedloss_orig_TIC04 <- TIC04 %>% 
  group_by(Year) %>% 
  summarize(mean(sed)) %>% 
  rename("annualsed_TIC04" = "mean(sed)")

TIC05$Year <- format(TIC05$date, format="%Y")

meanannualsedloss_orig_TIC05 <- TIC05 %>% 
  group_by(Year) %>% 
  summarize(mean(sed)) %>% 
  rename("annualsed_TIC05" = "mean(sed)")

meanannualsedloss_orig <- merge(meanannualsedloss_orig_TIC01,meanannualsedloss_orig_TIC02, by="Year")

meanannualsedloss_orig <- merge(meanannualsedloss_orig,meanannualsedloss_orig_TIC03, by="Year")

meanannualsedloss_orig <- merge(meanannualsedloss_orig,meanannualsedloss_orig_TIC04, by="Year")

meanannualsedloss_orig <- merge(meanannualsedloss_orig,meanannualsedloss_orig_TIC05, by="Year")






################################### Plots with changed Pusle ############################################# 

#Changing Pusle

#Row 1: for slp = 0.265, Pusle = 0.9 (Table 4:1-2)
#Row 2: for slp = 0.23, Pusle = 0.9
#Row 3: for slp = 0.17, Pusle = 0.8
#Row 4: for slp = 0.14, Pusle = 0.7
#Row 5: for slp = 0.10, Pusle = 0.6


#VSAParams1$Pusle <- c(0.9,0.9,0.8,0.7,0.6)

#View (MUSLE)
MUSLE1 <- MUSLE
glimpse(MUSLE1)
MUSLE1$Pusle <- c(0.9,0.9,0.8,0.7,0.6)

attach(MUSLE1)
MUSLE1$KCPLSCFRG118=11.8*Kusle*Cusle*Pusle*LS*CFRG
detach(MUSLE1)

VSAParams1=merge(VSAsol1,MUSLE1,by.x="TIClass",by.y="TIclass")
glimpse(VSAParams1)
#now plotting with changed Pusle

TIC01=TMWB
TIC02=TMWB
TIC03=TMWB
TIC04=TMWB
TIC05=TMWB
# For TIC01 CNavg=VSAParams1$CN[1] but confirm
TIC01 = CNmodel(CNmodeldf = TIC01, CNavg=VSAParams1$CN[1], 
                declat=myflowgage$declat,declon=myflowgage$declon)
TIC01$qpeak=TIC01$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec

TIC01$sed=(TIC01$Qpred*TIC01$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE1$KCPLSCFRG118[1]    # Eq. 4:1.1.1 SWAT Theory
VSAParams1$Pusle[1]
# View(TIC01)
glimpse(TIC01)

plot(TIC01$date, TIC01$sed, main = "", xlab="Date", ylab="Sediment Loss",
     col.axis="blue") 
title(main = "Time series of sediment loss for TI Class 1", sub = "Changed model with Pusle = 0.9",
      cex.main = 2,   font.main= 2, col.main= "black",
      cex.sub = 1, font.sub = 3, col.sub = "blue",
      col.lab ="darkblue"
) 

# For TIC02 CNavg=VSAParams1$CN[2] but confirm
TIC02 = CNmodel(CNmodeldf = TIC02, CNavg=VSAParams1$CN[2], 
                declat=myflowgage$declat,declon=myflowgage$declon)
TIC02$qpeak=TIC02$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec

TIC02$sed=(TIC02$Qpred*TIC02$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE1$KCPLSCFRG118[2]    # Eq. 4:1.1.1 SWAT Theory
VSAParams1$Pusle[2]
glimpse(TIC02)

plot(TIC02$date, TIC02$sed, main = "", xlab="Date", ylab="Sediment Loss",
     col.axis="blue") 
title(main = "Time series of sediment loss for TI Class 2", sub = "Changed model with Pusle = 0.9",
      cex.main = 2,   font.main= 2, col.main= "black",
      cex.sub = 1, font.sub = 3, col.sub = "blue",
      col.lab ="darkblue"
) 

# For TIC03 CNavg=VSAParams1$CN[3] but confirm
TIC03 = CNmodel(CNmodeldf = TIC03, CNavg=VSAParams1$CN[3], 
                declat=myflowgage$declat,declon=myflowgage$declon)
TIC03$qpeak=TIC03$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec

TIC03$sed=(TIC03$Qpred*TIC03$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE1$KCPLSCFRG118[3]    # Eq. 4:1.1.1 SWAT Theory
VSAParams1$Pusle[3]
glimpse(TIC03)

plot(TIC03$date, TIC03$sed, main = "", xlab="Date", ylab="Sediment Loss",
     col.axis="blue") 
title(main = "Time series of sediment loss for TI Class 3", sub = "Changed model with Pusle = 0.8",
      cex.main = 2,   font.main= 2, col.main= "black",
      cex.sub = 1, font.sub = 3, col.sub = "blue",
      col.lab ="darkblue"
) 
# For TIC04 CNavg=VSAParams1$CN[4] but confirm
TIC04 = CNmodel(CNmodeldf = TIC04, CNavg=VSAParams1$CN[4], 
                declat=myflowgage$declat,declon=myflowgage$declon)
TIC04$qpeak=TIC04$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec

TIC04$sed=(TIC04$Qpred*TIC04$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE1$KCPLSCFRG118[4]    # Eq. 4:1.1.1 SWAT Theory
VSAParams1$Pusle[4]
glimpse(TIC04)

plot(TIC04$date, TIC04$sed, main = "", xlab="Date", ylab="Sediment Loss",
     col.axis="blue") 
title(main = "Time series of sediment loss for TI Class 4", sub = "Changed model with Pusle = 0.7",
      cex.main = 2,   font.main= 2, col.main= "black",
      cex.sub = 1, font.sub = 3, col.sub = "blue",
      col.lab ="darkblue"
) 

# For TIC05 CNavg=VSAParams1$CN[5] but confirm
TIC05 = CNmodel(CNmodeldf = TIC05, CNavg=VSAParams1$CN[5], 
                declat=myflowgage$declat,declon=myflowgage$declon)
TIC05$qpeak=TIC05$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec

TIC05$sed=(TIC05$Qpred*TIC05$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE1$KCPLSCFRG118[5]    # Eq. 4:1.1.1 SWAT Theory
VSAParams1$Pusle[5]
glimpse(TIC05)

plot(TIC05$date, TIC05$sed, main = "", xlab="Date", ylab="Sediment Loss",
     col.axis="blue") 
title(main = "Time series of sediment loss for TI Class 5", sub = "Changed model with Pusle = 0.6",
      cex.main = 2,   font.main= 2, col.main= "black",
      cex.sub = 1, font.sub = 3, col.sub = "blue",
      col.lab ="darkblue"
) 

TotAgg_changedModel <- TIC01$sed + TIC02$sed + TIC03$sed + TIC04$sed + TIC05$sed
glimpse(TotAgg_changedModel)
plot(TIC05$date, TotAgg_changedModel, main = "", xlab="Date", ylab="Total Aggregate",
     col.axis="blue") 
title(main = "Time series of Total Aggregate", sub = "Changed model",
      cex.main = 2,   font.main= 2, col.main= "black",
      cex.sub = 1, font.sub = 3, col.sub = "blue",
      col.lab ="darkblue"
) 

