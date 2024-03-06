# Code for "Age-dependent shaping of the social environment in a long-lived seabird – a quantitative genetic approach"
# Unpublished manuscript, doi: tba
# Moiron M, Bouwhuis S

######################################################
# DATA ANALYSIS OF SOURCES OF VARIATION IN THE NUMBER OF NB
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

#Plot distribution of data (Figure 1A)
hist(Data$N.NB, xlab="Number of neighbours", main="", breaks=50)

#Random effects
Data$ID=as.factor(Data$ID)  
Data$animal = as.factor(Data$ID)
Data$YEAR=as.factor(Data$year)
Data$ISLE=as.factor(Data$island)

#Fixed effects
Data$LDZ=as.numeric(scale(Data$Laying.date))
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

# Loading pedigree
pedigree <- read.table("pedigree.txt",header=TRUE)
ped=prunePed(pedigree, Data$animal, make.base=TRUE)
my_inverse <- inverseA(ped)$Ainv

# Setting number of samples and iterations
nsamp <- 1000 
BURN <- 10000; THIN <- 10000
(NITT <- BURN + THIN*nsamp)

#Specifying prior
prior <-list(R = list(V = 1*0.1, nu = 1.002),
             G = list(G1 = list(V = 1*0.1, nu = 1.002, alpha.mu = 0, alpha.V = 1000),
                      G2 = list(V = 1*0.1, nu = 1.002, alpha.mu = 0, alpha.V = 1000),
                      G3 = list(V = 1*0.1, nu = 1.002, alpha.mu = 0, alpha.V = 1000),
                      G4 = list(V = 1*0.1, nu = 1.002, alpha.mu = 0, alpha.V = 1000)))

# Running model
mod<-MCMCglmm(DD~IDmeanAge*AgeSD+yearZ+colony.sizeZ,
               random=~animal+ID+YEAR+ISLE,
               family="poisson",
               data=Data, 
               ginverse=list(animal = my_inverse),
               nitt = NITT, thin = THIN, burnin = BURN,
               prior=prior, 
               verbose=TRUE,pr=FALSE)

#save(mod, file = "analysis_N.NB.rda")
#load("analysis_N.NB.rda")

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
                 mod$Sol[,"colony.sizeZ"],
                 mod$Sol[,"IDmeanAge"] - mod$Sol[,"AgeSD"]
)

tabfixef <- data.frame(param = c("Intercept","Mean Age", "Delta Age",
                                 "Mean age * Delta age",
                                 "Year","colony size", "average-delta diff"),
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

#Table of Random effects
allranef <- list(modG$VCV[,"animal"],
                 modG$VCV[,"ID"],
                 modG$VCV[,"YEAR"],
                 modG$VCV[,"ISLE"],
                 modG$VCV[,"units"],
                 (modG$VCV[,"animal"]+modG$VCV[,"ID"])/rowSums(modG$VCV),
                 modG$VCV[,"animal"]/rowSums(modG$VCV),
                 100*sqrt(modG$VCV [ , "animal"]) / (mean(Data$DD)))

predm0_bdate <- predict(modG, marginal=modG$Random$formula,posterior = "all")

tabranef <- data.frame(mode=round(unlist(lapply(allranef, posterior.mode)),3),
                       median=unlist(lapply(allranef, function(x){paste0("(",round(median(x),3),")")})),
                       CI=unlist(lapply(allranef, function(x){
                         paste0("[",round(HPDinterval(x)[1], 3), ", ",
                                round(HPDinterval(x)[2], 3), "]")})),
                       row.names = c("Additive genetic", "Permanent environment","Year", "Island", "Residual",
                                     "Repeatability", "Heritability", "Evolvability"))

tabranef

