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

myflowgage_id="01589330" #Dead Run at Franklintown, MD
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",
                         end_date = "2019-01-01")

#Normalizing Q in mm/day for basins
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

#Getting WX Data
WXData=FillMissWX(declat=myflowgage$declat, declon=myflowgage$declon,
                  StnRadius=30,minstns=10,date_min="2010-01-01",
                  date_max="2023-02-01",targElev=1,
                  method = "IDEW",alfa=2)

BasinData=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
# A few constants
MinTempCol <- "#0000ff"
  MaxTempCol <- "#ff0000"
    PCol <- "#000000"
      QCol <- PCol
      
      coeff=1
      p1= ggplot(BasinData, aes(x=date)) +
        geom_line( aes(y=MaxTemp), linewidth=1, color=MaxTempCol) + 
        geom_line( aes(y=MinTemp), linewidth=1, color=MinTempCol) + 
        geom_line( aes(y=Qmm), linewidth=1, color=QCol) +
        scale_y_continuous(
          # Features of the first axis
          name = "Temp(C)",
          
          # Adding a second axis and specifying its features
          sec.axis = sec_axis(~.*coeff, name="Depth(mm)")
        ) + 
        theme(
          axis.title.y = element_text(color = "black", size=13),
          axis.title.y.right = element_text(color = QCol, size=13)
        ) +
        ggtitle(myflowgage$gagename)
      
      p1

basestr=format(Sys.time(),"/%Y%m%d%H%M")
filename=paste0(mypdfdir,basestr,"graphHW2_Wk3.pdf")
pdf(filename) 
plot(p1)
dev.off()
print("file size")
print(file.size(filename))
print("I finished!")

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

mydemHW2_2=get_aws_terrain(locations=bboxpts@coords, 
                      z = 12, prj = proj4_utm,src ="aws",expand=1)
res(mydemHW2_2)
plot(mydemHW2_2)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")


# Writing obtained raster to a geotiff file that can be used with OS level hydrological models 
writeRaster(mydemHW2_2,filename = "mydemHW2_2.tif",overwrite=T)


zoomext=myflowgage$gagepoint_utm@coords
zoomext=rbind(zoomext,zoomext+res(mydemHW2_2)*100)
zoomext=rbind(zoomext,zoomext-res(mydemHW2_2)*100)
zoomext=SpatialPoints(zoomext,proj4string = crs_utm)  
zoom(mydemHW2_2,ext=zoomext)

zoomext2=myflowgage$gagepoint_utm@coords
zoomext2=rbind(zoomext2,zoomext2+res(mydemHW2_2)*10)
zoomext2=rbind(zoomext2,zoomext2-res(mydemHW2_2)*10)
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

z=raster("mydemHW2_2.tif")
plot(z)



# Pitremove: gives the filled DEM
system("mpiexec -n 2 pitremove -z mydemHW2_2.tif -fel mydemfel.tif")
fel=raster("mydemfel.tif")
plot(fel, main='Filled DEM of Dead Run: Full Extent')
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")
zoom(fel,ext=zoomext2)
title(main='Filled DEM of Dead Run: Zoomed to 10 cells')

#Difference between DEM and filled DEM

difference= fel-z
plot(difference)
title(main='Difference between DEM and Filled DEM of Dead Run: Full Extent')
zoom(difference,ext=zoomext2)
title(main='Difference between DEM and Filled DEM of Dead Run: Zoomed to 10 cells')

