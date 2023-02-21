TMWBmodel=function(TMWB,fcres=.3,FldCap=.45,WiltPt=.15,Z=1000){
  
  # Our TMWB Model
  detach(TMWB)
  
  TMWB$ET = TMWB$PET # in mm/day
  TMWB$AWC=(0.45-0.15)*1000 #Fld Cap = .45, Wilt Pt = .15, z=1000mm
  TMWB$dP = TMWB$P-TMWB$ET -TMWB$SNO + TMWB$SNOmlt 
  
  attach(TMWB)# Remember to detach or it gets ugly
  plot(date,Qmm,type = "l",col="black")
  lines(date,P,type = "l",col="red")
  lines(date,Qmm,type = "l",col="black") # We repeat to have Qmm on top of P
  lines(date,ET,type = "l",col="blue")
  legend("topright", c("P", "Qmm", "ET"), col = c("red", "black", "blue"),
         lty = 1:2, cex = 0.8)
  detach(TMWB) # IMPORTANT TO DETACH
  
  
  TMWB$AWC=(0.45-0.15)*1000 #Fld Cap = .45, Wilt Pt = .15, z=1000mm
  
  
  TMWB$AW=NA  #Assigns all values in column with “NA” (Not available)
  TMWB$AW[1]=250
  TMWB$Excess=NA
  TMWB$Excess[1]=0
  head(TMWB)
  
  # Here we go looping through our functions….
  
  attach(TMWB)
  for (t in 2:length(date)){
    if (dP[t]< 0) {  
      values<-soildrying(AW[t-1],dP[t],AWC[t])
    } else if (AW[t-1]+dP[t]>AWC[t]) {
      values<-soil_wetting_above_capacity(AW[t-1],dP[t],AWC[t])
    } else {
      values<-soilwetting (AW[t-1],dP[t],AWC[t])
    }
    AW[t]<-values[1]
    Excess[t]<-values[2]
  }
  
  detach(TMWB)
  TMWB$AW <-AW
  TMWB$Excess<-Excess
  rm(list=c("AW","Excess"))
  
  # Calculate Watershed Storage and River Discharge: 
  TMWB$Qpred=NA
  TMWB$Qpred[1]=0
  TMWB$S=NA
  TMWB$S[1]=0
  
  attach(TMWB)
  
  for (t in 2:length(date)){
    S[t]=S[t-1]+Excess[t]     
    Qpred[t]=fcres*S[t]
    S[t]=S[t]-Qpred[t]
  }
  detach(TMWB) # IMPORTANT TO DETACH
  TMWB$S=S
  TMWB$Qpred=Qpred # UPDATE vector BEFORE DETACHING
  rm(list=c("S","Qpred"))
  View(TMWB)
  dev.off()
  plot(TMWB$date,TMWB$Qmm,col="black",ylab ="Qmm(mm)",xlab="date",type="l")
  lines(TMWB$date,TMWB$Qpred,col="blue",type="l", 
        xlab = "", ylab = "")
  legend("topright", c("Qmm(mm)", "Qpred(mm)"), col = c("black", "blue"),
         lty = 1:2, cex = 0.8)
  
  
  TMWB$AWC=(FldCap-WiltPt)*Z # 
  TMWB$dP = 0 # Initializing Net Precipitation
  TMWB$ET = 0 # Initializing ET
  TMWB$AW = 0 # Initializing AW
  TMWB$Excess = 0 # Initializing Excess
  
  
  # Loop to calculate AW and Excess
  attach(TMWB)
  for (t in 2:length(AW)){
    # This is where Net Precipitation is now calculated
    # Do you remember what Net Precip is? Refer to week 2 notes
    # Update this to reflect the ET model described above
    
    ET[t] = (AW[t-1]/AWC[t-1])*PET[t] # New Model
    dP[t] = P[t] - ET[t] + SNOmlt[t] - SNOfall[t] 
    # From here onward, everything is the same as Week2’s lab
    if (dP[t]<=0) {
      values<-soildrying(AW[t-1],dP[t],AWC[t])
    } else if((dP[t]>0) & (AW[t-1]+dP[t])<=AWC[t]) {
      values<-soilwetting(AW[t-1],dP[t],AWC[t])
    } else {
      values<-soil_wetting_above_capacity(AW[t-1],dP[t],AWC[t])
    }
    AW[t]<-values[1]
    Excess[t]<-values[2]
    print(t)
  }
  TMWB$AW=AW
  TMWB$Excess=Excess
  TMWB$dP=dP
  TMWB$ET=ET
  rm(list=c("AW","dP","ET", "Excess"))
  detach(TMWB) # IMPORTANT TO DETACH
  
  # Calculate Watershed Storage and River Discharge, S and Qpred, playing with the reservoir coefficient to try to get Qpred to best match Qmm
  
  TMWB$Qpred=NA
  TMWB$Qpred[1]=0
  TMWB$S=NA
  TMWB$S[1]=0
  attach(TMWB)
  
  for (t in 2:length(date)){
    S[t]=S[t-1]+Excess[t]     
    Qpred[t]=fcres*S[t]
    S[t]=S[t]-Qpred[t]
  }
  TMWB$S=S
  TMWB$Qpred=Qpred # UPDATE vector BEFORE DETACHING
  detach(TMWB) # IMPORTANT TO DETACH
  rm(list=c("Qpred","S"))
  return(TMWB)
}
#Make a plot that has Qmm, P,and Qpred over time
plot(TMWB$date,TMWB$P,col="black")
lines(TMWB$date,TMWB$Qmm,type = "l",col="red")
lines(TMWB$date,TMWB$Qpred,col="blue")
plot(TMWB$Qmm, TMWB$Qpred)

NSE=function(Qobs,Qsim){
  return(1-sum((Qobs-Qsim)^2,na.rm=TRUE)/sum((Qobs-mean(Qobs, na.rm=TRUE))^2, na.rm=TRUE))
}