# Code for "Age-dependent shaping of the social environment in a longlived seabird – a quantitative genetic approach"
# Unpublished manuscript, doi: 10.1098/rstb.2022.0465
# Moiron M, Bouwhuis S

######################################################
#DATA ANALYSIS OF INDIVIDUAL VARIATION IN THE CHANGE OF NUMBER OF NB ALONG AN AGE GRADIENT
######################################################

# Loading packages
library(tidyr)
library(dplyr)
library (broom)
library(nadiv)
library(MCMCglmm)

# Loading phenotypic data
Data <- read.table("social_ageing_data.txt", header=TRUE)

#Response variable
Data$N.NB<-as.numeric(Data$n.neigbors.4w.ahead.2m) #social data, number of NB in a 2m radious

#Random effects
Data$ID=as.factor(Data$ID)  
Data$animal = as.factor(Data$ID)
Data$YEAR=as.factor(Data$year)
Data$ISLE=as.factor(Data$island)

#Fixed effects
Data$LDZ=as.numeric(scale(Data$Laying.date))
Data$AGE=as.numeric(scale(Data$age))
Data$yearZ=as.numeric(scale(Data$year))
Data$colony.sizeZ=as.numeric(scale(Data$Npairs))

#Estimate mean and delta age values
#Extract mean individual values of age
meansAge <- as.data.frame(tapply(Data$age, Data$ID, mean, rm.na=TRUE))
ID <- as.data.frame(rownames(meansAge))
df1 <- cbind(ID, meansAge)
names(df1) <- c("ID","IDmeanAge")
df1$IDmeanAge=as.numeric(df1$IDmeanAge)

#Merge both dataframes (individual mean dataframe and original dataframe)
df2 <- as.data.frame(merge(Data, df1,by="ID"))

#Extract deviations from individual mean values of Temperature
df2$AgeSD <- as.numeric(df2$age-df2$IDmeanAge)
Data=df2

#Specifying heterogeneous residuals
nblocks5 <-5	# number of 'residual blocks'
Data$envclass5 <- as.numeric(arules::discretize(Data$year,breaks= nblocks5, method='frequency'))
Data$envclass5=as.factor(Data$envclass5)

# Setting number of samples and iterations
nsamp <- 1000 
BURN <- 10000; THIN <- 10000
(NITT <- BURN + THIN*nsamp)

#Specifying prior
prior <- list(R = list(V = diag(5), nu = 1),
              G = list(G1 = list(V = diag(2)*0.1, nu= 3.002,alpha.mu=c(0,0), alpha.V=diag(2)*625),
                       G2 = list(V = diag(1)*0.1, nu = 1.002, alpha.mu=0, alpha.V=1000),
                       G3 = list(V = diag(1)*0.1, nu = 1.002, alpha.mu=0, alpha.V=1000)))

# Running model
mod<-MCMCglmm(N.NB~IDmeanAge*AgeSD+yearZ+colony.sizeZ,
               random=~us(1+AGE):ID+YEAR+ISLE,
               rcov= ~idh(envclass5):units,    
               data=Data,family="poisson",
               nitt = NITT, thin = THIN, burnin = BURN,
               prior=prior,verbose=TRUE, pr=FALSE)

#save(mod, file = "Individual.RR.analysis_N.NB.rda")
#load("Individual.RR.analysis_N.NB.rda")

summary(mod)
mod$DIC 

# Assesing convergence convergence and auto-correlation of posterior distributions
plot(mod$VCV)
effectiveSize(mod$VCV)
heidel.diag(mod$VCV)
autocorr.diag(mod$VCV)

#Fixed effects
posterior.mode(mod$Sol)
HPDinterval(mod$Sol)
plot(mod$Sol)

#Table of Fixed effects
allfixef <- list(mod$Sol[,"(Intercept)"],
                 mod$Sol[,"IDmeanAge"],
                 mod$Sol[,"AgeSD"],
                 mod$Sol[,"IDmeanAge:AgeSD"],
                 mod$Sol[,"yearZ"],
                 mod$Sol[,"colony.sizeZ"]
)

tabfixef <- data.frame(param = c("Intercept","Mean Age", "Delta Age",
                                 "Mean age * Delta age",
                                 "Year","colony size"),
                       mode=round(unlist(lapply(allfixef, posterior.mode)),3),
                       median=unlist(lapply(allfixef, function(x){
                         paste0("(",round(median(x),3),")")})),
                       CI=unlist(lapply(allfixef, function(x){
                         paste0("[",round(HPDinterval(x)[1],3), ", ",
                                round(HPDinterval(x)[2], 3), "]")})))
tabfixef

#Random effects
round(posterior.mode(mod$VCV),3)
round(HPDinterval(mod$VCV), 3)    
plot(mod$VCV)     

#Random effects
allranef <- list(
  mod$VCV[,"(Intercept):(Intercept).ID"],
  mod$VCV[,"AGE:AGE.ID"],
  mod$VCV[,"AGE:(Intercept).ID"],
  mod$VCV[, 2]/sqrt(modG$VCV[, 1] * modG$VCV[, 4]),
  mod$VCV[,"YEAR"],
  mod$VCV[,"ISLE"],
  mod$VCV[,"envclass51.units"],
  mod$VCV[,"envclass52.units"],
  mod$VCV[,"envclass53.units"],
  mod$VCV[,"envclass54.units"],
  mod$VCV[,"envclass55.units"])

predm0_bdate <- predict(mod, marginal=mod$Random$formula,posterior = "all")

tabranef <- data.frame(mode=round(unlist(lapply(allranef, posterior.mode)),3),
                       median=unlist(lapply(allranef, function(x){paste0("(",round(median(x),3),")")})),
                       CI=unlist(lapply(allranef, function(x){
                         paste0("[",round(HPDinterval(x)[1], 3), ", ",
                                round(HPDinterval(x)[2], 3), "]")})),
                       row.names = c("Va Intercepts","Va Slope", "Va Covariance", "Va Correlation",
                                     "Year","Island ID",
                                     "Residuals1","Residuals2","Residuals3","Residuals4","Residuals5"))


tabranef 
