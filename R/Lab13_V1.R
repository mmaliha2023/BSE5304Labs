if (!require("pacman")) install.packages("pacman")
pacman::p_load(meteoForecast)
?meteoForecast::getPointDays

############################## Precipitation #####################################################
gfsvars = grepVar('precip', service='gfs', complete=TRUE)


today <- Sys.Date()
myflowgage_id="0158397967"  # Mine Bank Run at MD

myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",
                         end_date = today)


library(lattice)

## Multiple variables
Precipvars <- getPoint(c(myflowgage$declon,myflowgage$declat),
        vars = "Total_precipitation_surface_Mixed_intervals_Accumulation",
        service ="gfs",
        day = today-1)

GFSprecip <- aggregate(Precipvars,as.Date(time(Precipvars)),sum)
plot(GFSprecip, xlab='Date',ylab='Total Precipitation (cm)', main="Precipitation Forecast")
# xyplot(vars)


#################################### Temperature ###############################################
gfsvarstemp = grepVar('temp', service='gfs', complete=TRUE)

Tempvars <- getPoint(c(myflowgage$declon,myflowgage$declat),
                 vars = "Temperature_surface",
                 service ="gfs",
                 day = today)

GFSMaxTemp <- aggregate(Tempvars,as.Date(time(Tempvars)),max)
GFSMinTemp <- aggregate(Tempvars,as.Date(time(Tempvars)),min)

plot(GFSMinTemp, ylim=c(276,315),col='blue', xlab='Date', ylab='Temperature(K)', main='Temperature Forecast')
lines(GFSMaxTemp, col='red')
legend("topleft", legend = c("Minimum", "Maximum"), col = c("blue", "red"), lty = 1)
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
source("https://raw.githubusercontent.com/Rojakaveh/FillMissWX/main/FillMissWX.R")
WXData=FillMissWX(declat=myflowgage$declat, declon=myflowgage$declon,
                  StnRadius=30,minstns=10,date_min="2015-01-01",
                  date_max=today,targElev=myflowgage$elev,
                  method = "IDW",alfa=2)

BasinData=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")

GFSprecip_df <- data.frame(Date=index(GFSprecip), Value = coredata(GFSprecip))
colnames(GFSprecip_df) <- c("date", "P")
GFSMaxTemp_df <- data.frame(Date=index(GFSMaxTemp), Value = coredata(GFSMaxTemp))
colnames(GFSMaxTemp_df) <- c("date", "MaxTemp")
GFSMaxTemp_df$MaxTemp <- GFSMaxTemp_df$MaxTemp - 273.15
GFSMinTemp_df <- data.frame(Date=index(GFSMinTemp), Value = coredata(GFSMinTemp))
colnames(GFSMinTemp_df) <- c("date", "MinTemp")
GFSMinTemp_df$MinTemp <- GFSMinTemp_df$MinTemp - 273.15
########################## Merging Precipitation ###################################################
WXData_Precip <- data.frame(WXData$date, P = WXData$P * 2.54)
colnames(WXData_Precip) <- c("date","P")
WXData_Precip <- na.omit(WXData_Precip)

Precip_forecast <- rbind(WXData_Precip,GFSprecip_df)
Precip_forecast <- Precip_forecast[-1,]
plot(Precip_forecast, type='l', xlab='Year', ylab='Total Precipitation in cm', main='Precipitation through the years')

Precip_2023 <- Precip_forecast[Precip_forecast$date >= "2023-01-01",]
plot(Precip_2023, type='l',xlab='Year', ylab='Total Precipitation in cm', main='Precipitation in 2023')

################################################## Merging Temperature ###############################
WXData_MaxTemp <- data.frame(WXData$date,WXData$MaxTemp)
WXData_MaxTemp <- na.omit(WXData_MaxTemp)
colnames(WXData_MaxTemp) <- c("date","MaxTemp")

MaxTemp_forecast <- rbind(WXData_MaxTemp,GFSMaxTemp_df)
MaxTemp_forecast <- MaxTemp_forecast[-1,]
#plot(MaxTemp_forecast, type='l', xlab='Year', ylab='Maximum Temperature in Celcius', main='Maximum temperature through the years')

MaxTemp_2023 <- MaxTemp_forecast[MaxTemp_forecast$date >= "2023-01-01",]
#plot(MaxTemp_2023, type='l',xlab='Year', ylab='Maximum Temperature in Celcius', main='Maximum temperature in 2023')

WXData_MinTemp <- data.frame(WXData$date,WXData$MinTemp)
WXData_MinTemp <- na.omit(WXData_MinTemp)
colnames(WXData_MinTemp) <- c("date","MinTemp")

MinTemp_forecast <- rbind(WXData_MinTemp,GFSMinTemp_df)
MinTemp_forecast <- MinTemp_forecast[-1,]
plot(MinTemp_forecast, type='l', xlab='Year', ylab='Temperature in Celcius', main='Temperature through the years', col='blue')
lines(MaxTemp_forecast, col='red')
legend("topleft", legend = c("Minimum", "Maximum"), col = c("blue", "red"), lty = 1)

MinTemp_2023 <- MinTemp_forecast[MinTemp_forecast$date >= "2023-01-01",]
plot(MinTemp_2023, type='l',xlab='Year', ylab='Temperature in Celcius', main='Temperature in 2023', col= 'blue')
lines(MaxTemp_2023, col='red')
legend("topleft", legend = c("Minimum", "Maximum"), col = c("blue", "red"), lty = 1)

dir.create("~/pngs")
setwd("~/pngs")
graphdir="~/pngs"
png(paste0(graphdir,"TempForecast.png"))





























