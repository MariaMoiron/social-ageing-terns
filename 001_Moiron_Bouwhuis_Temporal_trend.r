# Code for "Age-dependent shaping of the social environment in a longlived seabird – a quantitative genetic approach"
# Unpublished manuscript, doi: tba
# Moiron M, Bouwhuis S

######################################################
# DATA ANALYSIS OF TEMPORAL TREND IN NUMBER OF NB
######################################################

# Loading packages
library(tidyr)
library(dplyr)
library(ggplot2)

# Loading phenotypic data
data <- read.table("social_ageing_data.txt", header=TRUE)

#Response variable
data$N.NB<-as.numeric(data$n.neigbors.4w.ahead.2m) #social data, number of NB in a 2m radious

#Covariate
data$year = as.numeric(data$year)  

#Get statististics
Yrs<-min(Data$year):max(Data$year)
trait.means<-tapply(Data$N.NB,Data$year,mean,na.rm=T)
standard_error <- function(x) sd(x) / sqrt(length(x))
trait.se<-tapply(Data$N.NB,Data$year,standard_error)

# Plotting linear temporal trend (Figure 2)
ye <- rgb(0.128,0.128,0.128,0.5) #set colour

par(mar=c(5,6,2,2))
  plot(Yrs, trait.means,xlab= "Year",ylab= expression("Average number of neighbours"),
       pch=21,bg=1,cex.lab=3,col=ye,cex.axis=2,
       ylim=c(1,12),
       cex=2)
  segments(Yrs,trait.means-trait.se,Yrs,trait.means+trait.se, col = ye)

mbd <- lm(size ~ 1 + year,data = Data, na.action = "na.omit")
m0<-lm(trait.means~Yrs, na.action = "na.omit")
y<-predict(mbd, data.frame(year=Yrs),interval="confidence")
lines( Yrs, y[,1] , col=ye, lwd=5,bg=2)
lines( Yrs, y[,2],lty=2 , col=ye, lwd=4,bg=2)
lines( Yrs, y[,3],lty=2 , col=ye, lwd=4,bg=2)

# Prediction of change per year
coef(mbd)[2]
#0.2198639

# Prediction of the overall change across years:
#Estimate:
(max(Yrs)-min(Yrs))*coef(mbd)[2]
#6.15619 

#Confidence interval:
CImbd <- confint(mbd)[2,]
(max(Yrs)-min(Yrs))*
  (CImbd[1])
## 2.5 %
## 5.524219 

(max(Yrs)-min(Yrs))*
  (CImbd[2])
## 97.5 %
## 6.788161

