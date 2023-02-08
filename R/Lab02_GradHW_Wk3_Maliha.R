Sys.getenv()
Sys.getenv("HOME")
myhomedir=Sys.getenv("HOME")
# Settin up git directory
getwd()
mygitdir=getwd()

# Location for homework PDFs
mypdfdir=paste0(mygitdir,"/pdfs")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,dplyr,patchwork,rnoaa)
pacman::p_load(operators,topmodel,DEoptim,soilDB,sp,curl,httr,
               rnoaa,raster,shapefiles,rgdal,elevatr,terra,progress,lubridate)
system("git config --global user.email 'mushtari@vt.edu' ") 
system("git config --global user.name 'Mushtari Maliha' ")

datadir=paste0(myhomedir,"/data")
dir.create(datadir,recursive = T)
srcdir=paste0(myhomedir,"/src")
dir.create(srcdir,recursive = T)

setwd(srcdir)
system("svn checkout svn://scm.r-forge.r-project.org/svnroot/ecohydrology/"); 
install.packages(c("ecohydrology/pkg/EcoHydRology/"),repos = NULL)
pacman::p_load(EcoHydRology)

setwd(datadir)

declon <- c(-80.4472973010981)
declat <- c(37.208226443212645)
area <- c(15.02)

myflowgage <- data.frame(declat, declon, area)

glimpse(myflowgage)

trunc((180+myflowgage$declon)/6+1)
proj4_utm = paste0("+proj=utm +zone=", trunc((180+myflowgage$declon)/6+1), " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)

# Obtaining Lat/Lon
proj4_ll = "+proj=longlat"

# Defining “Coordinate  Reference Systems” or CRS in future geographic manipulations. 
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
print(crs_ll)
print(crs_utm)

myflowgage$area

latlon <- cbind(myflowgage$declon,myflowgage$declat)
myflowgage$gagepoint_ll <- SpatialPoints(latlon)
proj4string(myflowgage$gagepoint_ll)=proj4_ll
myflowgage$gagepoint_utm=spTransform(myflowgage$gagepoint_ll,crs_utm)

# Opening up maps.google.com to estimate area/lengths
url=paste0("https://www.google.com/maps/@",
           myflowgage$declat,",",myflowgage$declon,",18z")
browseURL(url)

# We are going to over estimate our area
sqrt(myflowgage$area)   # assuming square watershed
# Multiplying the area to search for DEM
sqrt(myflowgage$area*8)
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

streamlabDEM=get_aws_terrain(locations=bboxpts@coords, 
                      z = 12, prj = proj4_utm,src ="aws",expand=1)
res(streamlabDEM)
plot(streamlabDEM)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")


# Writing obtained raster to a geotiff file that can be used with OS level hydrological models 
writeRaster(streamlabDEM,filename = "streamlabDEM.tif",overwrite=T)


zoomext=myflowgage$gagepoint_utm@coords
zoomext=rbind(zoomext,zoomext+res(streamlabDEM)*100)
zoomext=rbind(zoomext,zoomext-res(streamlabDEM)*100)
zoomext=SpatialPoints(zoomext,proj4string = crs_utm)  
zoom(streamlabDEM,ext=zoomext)

zoomext2=myflowgage$gagepoint_utm@coords
zoomext2=rbind(zoomext2,zoomext2+res(streamlabDEM)*10)
zoomext2=rbind(zoomext2,zoomext2-res(streamlabDEM)*10)
zoomext2=SpatialPoints(zoomext2,proj4string = crs_utm)  


plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")

#---In Terminal---
# cd ~/src/      # Set your directory to your home directory
#  git clone https://github.com/dtarb/TauDEM.git
# mkdir ~/src/TauDEM/bin
# cd ~/src/TauDEM/src
# make
# sed -i -e 's/MPI_Type_struct/MPI_Type_create_struct/g' linklib.h
# yes, this next line is very small font, but it is one line so...
# sed -i -e 's/MPI_Type_extent(MPI_LONG, \&extent)/MPI_Aint lb\;MPI_Type_get_extent(MPI_LONG, \&lb, \&extent)/g' linklib.h
# Now let's try make again!
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

# Setting working directory to your location
setwd(datadir)

z=raster("streamlabDEM.tif")
plot(z)



# Pitremove: gives the filled DEM
system("mpiexec -n 2 pitremove -z streamlabDEM.tif -fel mydemfel.tif")
fel=raster("mydemfel.tif")
plot(fel, main='Filled DEM of StREAM Lab: Full Extent')
zoom(fel,ext=zoomext)
title(main='Filled DEM of StREAM Lab: Zoomed to 100 cells')
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")

#Difference between DEM and filled DEM

difference= fel-z
plot(difference)
title(main='Difference between DEM and Filled DEM of StREAM Lab: Full extent')
zoom(difference,ext=zoomext)
title(main='Difference between DEM and Filled DEM of StREAM Lab: Zoomed to 100 cells')
