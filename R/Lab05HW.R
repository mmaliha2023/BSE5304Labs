#?DEoptim
#citation("DEoptim")
install.packages("DEoptim")
library("DEoptim")
Rosenbrock <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  return(100 * (x2 - x1 * x1)^2 + (1 - x1)^2)
}
Rosenbrock(c(3.7,5))

## DEoptim searches for minima of the objective function between
## lower and upper bounds on each parameter to be optimized. Therefore
## in the call to DEoptim we specify vectors that comprise the
## lower and upper bounds; these vectors are the same length as the
## parameter vector.
lower <- c(-10,-10)
upper <- -lower

Rosenbrock(lower)
Rosenbrock(upper)
## run DEoptim and set a seed first for replicability
set.seed(1234)
outDEoptim=DEoptim(Rosenbrock, lower, upper)

## increase the population size
DEoptim(Rosenbrock, lower, upper, DEoptim.control(NP = 100))

## change other settings and store the output
outDEoptim <- DEoptim(Rosenbrock, lower, upper, DEoptim.control(NP = 80,
                                                                itermax = 400, F = 1.2, CR = 0.7))

plot(outDEoptim$member$bestmemit, col="black", xlim=c(-10,10), ylim=c(-10,10))
points(outDEoptim$member$bestmemit, col='red')

## plot the output
plot(outDEoptim)


# Things we repeat can be saved and accessed without messing up our working area
# https://raw.githubusercontent.com/mmaliha2023/BSE5304Labs/main/R/HW_Lab04_Wk5.R
# becomes: 
url="https://raw.githubusercontent.com/mmaliha2023/BSE5304Labs/main/R/HW_Lab04_Wk5.R"
# This will grab the solution for last weeks Lab03 Homework
download.file(url,"HW_Lab04_Wk5.R")
file.edit("HW_Lab04_Wk5.R")

# Grab out models for Snow and TMWB
# https://raw.githubusercontent.com/vtdrfuka/BSE5304Labs/main/R/TMWBFuncs.R
# becomes: 
url="https://raw.githubusercontent.com/vtdrfuka/BSE5304Labs/main/R/TMWBFuncs.R"
# This will grab the solution for last weeks Lab03 Homework
download.file(url,"TMWBFuncs.R")
file.edit("TMWBFuncs.R")
# I actually am starting to trust my snow model
url="https://raw.githubusercontent.com/vtdrfuka/BSE5304Labs/main/R/TISnow.R"
# This will grab the solution for last weeks Lab03 Homework
source(url)
download.file(url,"TISnow.R")
file.edit("TISnow.R")
#outTMWB = TMWBmodel(TMWBdf=TMWB)
#1-NSE(Yobs=outTMWB$Qmm, Ysim= outTMWB$Qpred )

TMWB=BasinData

detach(TMWBdf)
detach(WBData)

TMWBoptFunc <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  x4 <- x[4]
  x5 <- x[5]
  x6 <- x[6]
  x7 <- x[7]
  outTMWB=TMWBmodel(TMWBdf = TMWB,fcres = x1,Z = x2,SFTmp = x3,bmlt6=x4, bmlt12=x5,Tmlt=x6,Tlag=x7)
  return(1-NSE(Yobs = outTMWB$Qmm,Ysim = outTMWB$Qpred))
}
lower <- c(.01,300,1,.1,.1,1,.1)
upper <- c(.95,3000,6,5,5,6,1)
outDEoptim <- DEoptim(TMWBoptFunc, lower, upper, 
                      DEoptim.control(NP = 80,
                                      itermax = 100, F = 1.2, CR = 0.7))
par(mar = c(1, 1, 1, 1))
##Optimizing curve number model

CNmodeldf=BasinData
source("https://raw.githubusercontent.com/mmaliha2023/BSE5304Labs/main/CNmodel2.R")
download.file(url,"CNmodel2.R")
file.edit("CNmodel2.R")
#CNmodelnew=CNModel(CNmodeldf=BasinData,CNavg = 75, func_z=1000,fnc_fcres=.2) 
#NSE(CNmodelnew$Qpred,CNmodelnew$Qmm)
#CNmodeloptFunc <- function(x){
 # x1 <- x[1]
 # x2 <- x[2]
#  x3 <- x[3]
 # outCN=CNModel(BasinData,CNavg = x3, func_z=x2,fnc_fcres=x1)
#  return(1-NSE(Yobs=outTMWB$Qmm, Ysim=outTMWB$Qpred))
#}
#lower <- c(0.01, 300, 70)
#upper <- c(0.95, 3000, 90)
#outDEoptim <- DEoptim(CNmodeloptFunc,lower,upper, DEoptim.control(NP=80, itermax=70, F=1.2, CR=0.7))

CNmodeloptFunc <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  x4 <- x[4]
  x5 <- x[5]
  outCN=CNModel(BasinData,CNavg = x1,IaFrac = x2,func_DAWC=x3,func_z=x4,fnc_fcres=x5)
  return(1-NSE(Yobs=outTMWB$Qmm, Ysim=outTMWB$Qpred))
}
lower <- c(70,0.01, 0.01,100,0.01)
upper <- c(90,0.1, 0.05,1000,0.95)
outDEoptim <- DEoptim(CNmodeloptFunc,lower,upper, DEoptim.control(NP=80, itermax=100, F=1.2, CR=0.7))
detach(CNmodeldf)
detach(TMWBdf)


