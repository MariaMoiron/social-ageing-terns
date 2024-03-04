# Code for "Age-dependent shaping of the social environment in a longlived seabird – a quantitative genetic approach"
# Unpublished manuscript, doi: tba
# Moiron M, Bouwhuis S

# The code provided here is sufficient to replicate the results presented in the above paper


######################################################
# DATA ANALYSIS OF GENETIC VARIATION IN THE CHANGE OF NUMBER OF NB ALONG AN AGE GRADIENT
######################################################

# Loading packages
library(tidyr)
library(dplyr)
library (broom)
library(nadiv)
library(MCMCglmm)

# Loading phenotypic data
data <- read.table("social_ageing_data.txt", header=TRUE)

#Response variable
data$N.NB<-as.numeric(data$Data$n.neigbors.4w.ahead.2m) #social data, number of NB in a 2m radious

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

# Loading pedigree
pedigree <- read.table("pedigree.txt",header=TRUE)
ped=prunePed(pedigree, Data$animal, make.base=TRUE)  #to prune the pedigree, i.e., only keep those records on the pedigree that have phenotypic information for our analyses
my_inverse <- inverseA(ped)$Ainv

#Specifying heterogeneous residuals
nblocks5 <-5	# number of 'residual blocks'

Data$envclass5 <- as.numeric(arules::discretize(Data$year,breaks= nblocks5, method='frequency'))
table(Data$envclass5)
Data$envclass5=as.factor(Data$envclass5)

# Setting number of samples and iterations
nsamp <- 1000 
BURN <- 10000; THIN <- 10000
(NITT <- BURN + THIN*nsamp)

#Specifying prior
prior <- list(R = list(V = diag(5), nu = 1),
              G = list(G1 = list(V = diag(2)*0.01, nu= 3.002,alpha.mu=c(0,0), alpha.V=diag(2)*625),
                       G2 = list(V = diag(2)*0.01, nu= 3.002,alpha.mu=c(0,0), alpha.V=diag(2)*625),
                       G3 = list(V = diag(1)*0.01, nu = 1, alpha.mu=0, alpha.V=625),
                       G4 = list(V = diag(1)*0.01, nu = 1, alpha.mu=0, alpha.V=625)))

# Running model
mod<-MCMCglmm(N.NB~IDmeanAge*AgeSD+yearZ+colony.sizeZ,
              random=~us(1+AGE):animal+us(1+AGE):ID+YEAR+ISLE,
              rcov= ~idh(envclass5):units, 
              data=Data,family="poisson",
              ginverse=list(animal = my_inverse),
               nitt = NITT, thin = THIN, burnin = BURN,
               prior=prior,verbose=TRUE, pr=TRUE)

#save(mod, file = "Genetic.RR.analysis_N.NB.rda")
#load("Genetic.RR.analysis_N.NB.rda")

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
allranef <- list(mod$VCV[,"(Intercept):(Intercept).animal"],
                 mod$VCV[,"AGE:AGE.animal"],
                 mod$VCV[,"AGE:(Intercept).animal"],
                 mod$VCV[, 2]/sqrt(mod$VCV[, 1] * mod$VCV[, 4]),
                 
                 mod$VCV[,"(Intercept):(Intercept).ID"],
                 mod$VCV[,"AGE:AGE.ID"],
                 mod$VCV[,"AGE:(Intercept).ID"],
                 mod$VCV[, 6]/sqrt(mod$VCV[, 5] * mod$VCV[, 8]),
                 
                 mod$VCV[,"YEAR"],
                 mod$VCV[,"ISLE"],
                 mod$VCV[,"envclass51.units"],
                 mod$VCV[,"envclass52.units"],
                 mod$VCV[,"envclass53.units"],
                 mod$VCV[,"envclass54.units"],
                 mod$VCV[,"envclass55.units"])


predm0_bdate <- predict(mod, marginal=mod$Random$formula,
                        posterior = "all")

tabranef <- data.frame(mode=round(unlist(lapply(allranef, posterior.mode)),3),
                       median=unlist(lapply(allranef, function(x){
                         paste0("(",round(median(x),3),")")})),
                       CI=unlist(lapply(allranef, function(x){
                         paste0("[",round(HPDinterval(x)[1], 3), ", ",
                                round(HPDinterval(x)[2], 3), "]")})),
                       
                       row.names = c("Va Intercepts","Va Slope", "Va Covariance", "Va Correlation",
                                     "Vpe Intercepts","Vpe Slope", "Vpe Covariance", "Vpe Correlation",
                                     "Year","Island ID",#"residuals"))
                                     "Residuals1","Residuals2","Residuals3","Residuals4","Residuals5"))


tabranef 

# Figure 3
df_rr_ind <- cbind(Data,
                   fit = predict(modG, marginal = NULL)) %>%
  group_by(animal,AGE) %>%
  summarise(fit = mean(fit))

# Plot predictions and overlay original data points
ye <- rgb(0.128,0.128,0.128,0.5)

plot=ggplot(df_rr_ind, aes(x =AGE, y = fit)) +
  geom_jitter(data = Data,width = 0.01,
              aes(y = DD),size=2,
              alpha = 0.4, col="grey") +
  scale_x_continuous(breaks = c(-1, 0,1,2,3, 4))+
  scale_y_continuous(breaks = c(0,10,20,30))+
  geom_smooth(aes(y = fit,group=animal,color=animal),
              method=lm,formula = y ~x,   # Add linear regression lines
              se=FALSE,    # Don't add shaded confidence region
              fullrange=FALSE, alpha=0.1)+
  scale_color_viridis_d(option = "magma")

plot+theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size=25, colour = "black"),
        axis.title=element_text(size=25),
        plot.margin = unit(c(0,0.5,0,0), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position= "none",
        legend.justification=c(0.9,0.8),
        legend.title = element_blank(),
        legend.background =element_rect(fill="gray95"))+
  labs(y = "Number of neighbours", x = "Age (standardized)")

