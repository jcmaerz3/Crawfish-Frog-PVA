#------------------------------------------------------------------#
###############Crawfish frog stage-based matrix PVA model#############
#################Coded by B.A. Crawford and J.C. Maerz##############
#for Terrell et al. Population dynamics of threatened crawfish frogs informs management decisions#
#------------------------------------------------------------------#
setwd("C:/Users/Brian/Dropbox/Collaborative research/Kinney et al crawfish frog paper/Crawfish frog PVA in R") 
rm(list=ls())

library("plyr")
library("dplyr")
library("reshape")
library("ggplot2")
library("gridExtra")
library("grid")
library("RColorBrewer")
library("stringr")
library("fitdistrplus")

#--------------------------------------------------#
###Import data
#--------------------------------------------------#
vitals.orig <- read.csv("craw_vital_rates.csv")   # Vital rates

#Because tadpole survival was modeled separately for Nates v. Cattail Pond, we coded the model to first run Nate's Pond using the density dependent tadpole survival and then the model runs for Cattail Pond before merging the two outputs into a single file

#-----------------------#
### BEGIN NATES MODEL ###
#-----------------------#

vitals <- vitals.orig[which(vitals.orig$Site=="Nates"),]

# Extract mean vital rates from input data table
phi.emb <- vitals$Mean[vitals$Rate=="phi.emb"]
phi.tad <- vitals$Mean[vitals$Rate=="phi.tad"]
phi.juv <- vitals$Mean[vitals$Rate=="phi.juv"]
phi.adult.1 <- vitals$Mean[vitals$Rate=="phi.adult.1"]
phi.adult.2 <- vitals$Mean[vitals$Rate=="phi.adult.2"]
p <- vitals$Mean[vitals$Rate=="repro.p"]
q <- vitals$Mean[vitals$Rate=="repro.q"]
repro.age.2.f <- vitals$Mean[vitals$Rate=="repro.age.2.f"]
repro.age.3.f <- vitals$Mean[vitals$Rate=="repro.age.3.f"]
repro.age.4.f <- vitals$Mean[vitals$Rate=="repro.age.4.f"]
clutch <- vitals$Mean[vitals$Rate=="clutch"]
prop.clutch.f <- vitals$Mean[vitals$Rate=="prop.clutch.f"]
gam.BNB <- vitals$Mean[vitals$Rate=="gam.BNB"]
gam.NBNB <- vitals$Mean[vitals$Rate=="gam.NBNB"]

# Extract SD vital rates from input data table
phi.tad.SD <- vitals$SD[vitals$Rate=="phi.tad"]
phi.juv.SD <- vitals$SD[vitals$Rate=="phi.juv"]
phi.adult.1.SD <- vitals$SD[vitals$Rate=="phi.adult.1"]
phi.adult.2.SD <- vitals$SD[vitals$Rate=="phi.adult.2"]
p.SD <- vitals$SD[vitals$Rate=="repro.p"]
q.SD <- vitals$SD[vitals$Rate=="repro.q"]
clutch.SD <- vitals$SD[vitals$Rate=="clutch"]
gam.BNB.SD <- vitals$SD[vitals$Rate=="gam.BNB"]
gam.NBNB.SD <- vitals$SD[vitals$Rate=="gam.NBNB"]

#Calculate stable age distribution for initial population sizes
fert <- clutch*prop.clutch.f
gam <- 1-gam.BNB
s1 <- phi.emb*phi.tad*(phi.juv^0.75)
s2 <- phi.juv*(1-p)
s3 <- phi.juv*(1-gam)*p
s4 <- phi.juv*gam*p
s5 <- phi.juv*(1-q)
s6 <- phi.juv*(1-gam)*q
s7 <- phi.juv*gam*q
s8 <- phi.adult.1*(1-gam)
s9 <- phi.adult.1*gam
s10 <- phi.juv*(1-gam)
s11 <- phi.juv*gam
s12 <- phi.adult.1*(1-gam)
s13 <- phi.adult.1*gam
s14 <- phi.adult.2*(1-gam)
s15 <- phi.adult.2*gam

fert2 <- fert*s4
fert3 <- fert*s7
fert4 <- fert*s9
fert5 <- fert*s11
fert6 <- fert*s13
fert7 <- fert*s15

# Matrix columns:
# 1) year 1: embryo, metamorph, juveniles
# 2) year 2: s2 = juveniles, not mature, s3 = juveniles, mature, not breed, s4 = juveniles, mature, breed
# 3) year 3: s5 = juveniles, not mature, s6 = juveniles, mature, not breed, s7 = juveniles, mature, breed, s8 = adults, not breed, s9 = adults, breed
# 4) year 4: s10 = juveniles, mature, not breed, s11 = juveniles, mature, breed, s12 = 1st adults, not breed, s13 = 1st adults, breed, s14 = 2nd adults, not breed, s15 = 2nd adults, breed

Lefkovitch <- c(0,     fert2,   fert3,  fert4, fert4,  fert5,  fert6,  fert6,  fert7, fert7,  fert6,  fert6, fert7,  fert7,  fert7, fert7,
                s1,    0,       0,      0,     0,      0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0, 
                0,     s2,      0,      0,     0,      0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     s3,      0,      0,     0,      0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,        
                0,     s4,      0,      0,     0,      0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       s5,     0,     0,      0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       s6,     0,     0,      0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       s7,     0,     0,      0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       0,      s8,    s8,     0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       0,      s9,    s9,     0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       0,      0,     0,      s10,    0,      0,      0,     0,      0,      0,     0,      0,      0,     0, 
                0,     0,       0,      0,     0,      s11,    0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       0,      0,     0,      0,      s12,    s12,    0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       0,      0,     0,      0,      s13,    s13,    0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       0,      0,     0,      0,      0,      0,      s14,   s14,    s12,    s12,   s14,    s14,    s14,   s14, 
                0,     0,       0,      0,     0,      0,      0,      0,      s15,   s15,    s13,    s13,   s15,    s15,    s15,   s15)
craw <- matrix(Lefkovitch,nrow=16,byrow=T)
craw_eigen <- eigen(craw)
stable_stage <- craw_eigen$vectors[,1]/sum(craw_eigen$vectors[,1]) 
stable_stage <-as.numeric(stable_stage)

# Immigration distribution
imm.dat <- c(3,5,8,13)  # Observed number of immigrant females to Nate's pond
# fit the negative binomial distribution
fit <- fitdist(imm.dat, "nbinom")

# get the fitted densities. mu and size from fit.
fitD <- rnbinom(1000, size=fit$estimate[1], mu=fit$estimate[2])

#--------------------------------------------------------------------------------------------------#
##Population Viability Analysis - stage-based model for 50 years and 1000 iterations per scenario###
#--------------------------------------------------------------------------------------------------#

iter <- 1000 # Number of iterations per scenario/strategy
n.yrs <- 50 # Number of years to project pop for each iteration
breed.prob <- 0.8 # Probability of Breeding in a Particular Year
breed.state <- array(dim=c(n.yrs+1,iter))
breed.state[,] <- sapply(1:iter, function (i) rbinom(n.yrs+1,1,breed.prob))

# Set up main pop abundance (N) array
N<-array(dim=c(n.yrs+1,18,iter))   # Main abundance array for 18 classes (16 stage classes described in above matrix + 1 column for total adult breeders + 1 column for total abundance)
dim(N)    #Check if dimensions are right (1=years, 2=stages, 3=iterations, 4=mgmt scenarios)
lambda<-array(dim=c(n.yrs,iter))
crash <- array(dim=c(iter))
mean.lambda.iter<-array(dim=c(iter))
mean.lambda<-array(dim=1)
perc.change.iter<-array(dim=c(iter))
perc.change<-array(dim=1)
persist<-array(dim=c(iter))
persistence<-array(dim=c(iter))
persistence.prob<-array(dim=1)
N.tads<-array(dim=c(n.yrs+1,iter))
N.met<-array(dim=c(n.yrs+1,iter))
N.imms<-array(dim=c(n.yrs+1,iter))
phi.tad.DD<-array(dim=c(n.yrs+1,iter))

#### Population Viability Model
N.AB.0<-36    # Starting adult pop size
Total0 <- round(N.AB.0/sum(stable_stage[c(5,8,10,12,14,16)]),0)
stage.N <- round(stable_stage*Total0,0)
N0 <- c(stage.N,sum(stage.N[c(5,8,10,12,14,16)]),Total0)

# Begin PVA simulations
for (i in 1:iter){
  # Begin forecasting population abundance over 30 years
  N[1,,i] <- N0                                                                        #Initial abundances
  for (t in 1:n.yrs){
    # Incorporate parametric uncertainty around sensitive parameters (initial abundance, terrestrial survival rates)
    phi.S.juv <- rnorm(1,phi.juv,phi.juv.SD)
    phi.S.juv <- ifelse(phi.S.juv<phi.juv-2*phi.juv.SD,phi.juv-2*phi.juv.SD,ifelse(phi.S.juv>phi.juv+2*phi.juv.SD,phi.juv+2*phi.juv.SD,phi.S.juv))  # Bounded juvenile survival to +/- 2 SDs
    phi.S.adult.2 <- rnorm(1,phi.adult.2,phi.adult.2.SD)
    phi.S.adult.2 <- ifelse(phi.S.adult.2>0.85,0.85,ifelse(phi.S.adult.2<0.2,0.2,phi.S.adult.2))
    phi.S.adult.1 <- 0.9089*phi.S.adult.2 - 0.1924
    phi.S.adult.1 <- ifelse(phi.S.adult.1>0.7,0.7,ifelse(phi.S.adult.1<0.15,0.15,phi.S.adult.1))
    gam.BNB.S <- rnorm(1,gam.BNB,gam.BNB.SD)
    gam.NBNB.S <- rnorm(1,gam.NBNB,gam.NBNB.SD)
    
    N.tads[t+1,i] <- rbinom(1,N[t,1,i],phi.emb)      # Eggs surviving to tadpoles                               
    phi.tad.DD[t+1,i] <- 0.0346 / (1 + (exp(0.085*(N[t,17,i]-27)) + rnorm(1,0,(0.008 * exp(-(((N[t,17,i]-25)^2)/(2*9^2))))))) #Density dependence at Nates using negative sigmoid function based on # breeding females in year prior
    N.met[t+1,i] <- rbinom(1,N.tads[t+1,i],phi.tad.DD[t+1,i])          # Metamorphs
    N[t+1,2,i] <- rbinom(1,N.met[t+1,i],phi.S.juv^0.75)                # Juvs (age 1) = first year survival through embryo, tadpole, and metamorph/juvenile stages
    
    N[t+1,3,i] <- rbinom(1,N[t,2,i],phi.S.juv*(1-p))    #Juvs (age 2) = juvs that didn't mature at age 2 (p) and survived
    N[t+1,4,i] <- rbinom(1,N[t,2,i],phi.S.juv*p*(gam.BNB.S))              #Juvs to ads, nb (age 2) = juvs that matured, survived, and didn't breed
    N[t+1,5,i] <- rbinom(1,N[t,2,i],phi.S.juv*p*(1-gam.BNB.S))                  #Juvs to ads, b (age 2) = juvs that matured, survived, and bred
    
    N[t+1,6,i] <- rbinom(1,N[t,3,i],phi.S.juv*(1-q))    #Juvs (age 3) = juvs that didn't mature at age 3 (q) and survived
    N[t+1,7,i] <- rbinom(1,N[t,3,i],phi.S.juv*q*(gam.BNB.S))              #Juvs to ads, nb (age 3) = juvs that matured, survived, and didn't breed
    N[t+1,8,i] <- rbinom(1,N[t,3,i],phi.S.juv*q*(1-gam.BNB.S))                  #Juvs to ads, b (age 3) = juvs that matured, survived, and bred
    N[t+1,9,i] <- rbinom(1,N[t,4,i],phi.S.adult.1*(gam.NBNB.S)) + rbinom(1,N[t,5,i],phi.S.adult.1*(gam.BNB.S))     #1st yr Ads, nb (age 3) = ads that survived and didn't breed
    N[t+1,10,i] <- rbinom(1,N[t,4,i],phi.S.adult.1*(1-gam.NBNB.S)) + rbinom(1,N[t,5,i],phi.S.adult.1*(1-gam.BNB.S))         #1st yr Ads, b (age 3) = ads that survived and bred
    
    N[t+1,11,i] <- rbinom(1,N[t,6,i],phi.S.juv*(gam.BNB.S))              #Juvs to ads, nb (age 4) = juvs that matured, survived, and didn't breed
    N[t+1,12,i] <- rbinom(1,N[t,6,i],phi.S.juv*(1-gam.BNB.S))                  #Juvs to ads, b (age 4) = juvs that matured, survived, and bred
    imms <- rnbinom(100, size=fit$estimate[1], mu=fit$estimate[2]) # Draw random numbers of breeding adult immigrants
    N.imms[t+1,i] <- imms[min(which(imms<=max(imm.dat)))] # Draw annual number of breeding adult immigrants, restricted to being less than or equal to the max observed number of immigrants
    N[t+1,13,i] <- rbinom(1,sum(N[t,c(7,11),i]),phi.S.adult.1*(gam.NBNB.S)) + rbinom(1,sum(N[t,c(8,12),i]),phi.S.adult.1*(gam.BNB.S))            #1st yr Ads, nb (age 4) = 1st yr ads that survived and didn't breed
    N[t+1,14,i] <- rbinom(1,sum(N[t,c(7,11),i]),phi.S.adult.1*(1-gam.NBNB.S)) + rbinom(1,sum(N[t,c(8,12),i]),phi.S.adult.1*(1-gam.BNB.S))        #1st yr Ads, b (age 4) = 1st yr ads that survived and bred
    N[t+1,14,i] <- N[t+1,14,i] + N.imms[t+1,i] #1st yr Ads, b (age 4) = 1st yr ads that survived and bred + immigrants
    N[t+1,15,i] <- rbinom(1,sum(N[t,c(9,13,15),i]),phi.S.adult.2*(gam.NBNB.S)) + rbinom(1,sum(N[t,c(10,14,16),i]),phi.S.adult.2*(gam.BNB.S))                #2nd+ yr Ads, nb (age 4) = 2nd+ yr ads that survived and didn't breed
    N[t+1,16,i] <- rbinom(1,sum(N[t,c(9,13,15),i]),phi.S.adult.2*(1-gam.NBNB.S)) + rbinom(1,sum(N[t,c(10,14,16),i]),phi.S.adult.2*(1-gam.BNB.S))                   #2nd+ yr Ads, b (age 4) = 2nd+ yr ads that survived and bred
  
    N[t+1,1,i] <- round(sum(N[t,c(5,8,10,12,14,16),i]) * clutch*prop.clutch.f * breed.state[t+1,i],0)  # nobreed is indicator variable simulating different lengths/frequencies of no recruitment years (e.g., droughts)
    
    # Track total abundance, lambda, time to extinction
    N[t+1,17,i] <- sum(N[t+1,c(5,8,10,12,14,16),i])   #Annual adult breeder abundance
    N[t+1,17,i]<-(N[t,17,i]>=2)*N[t+1,17,i]        #Makes N= 0 if adult abundance falls below 2 (quasi-extinction threshold)
    N[t+1,18,i] <- sum(N[t+1,1:16,i])                 #Annual total abundance
    N[t+1,18,i]<-(N[t,17,i]>=2)*N[t+1,18,i]        #Makes N= 0 if adult abundance falls below 2 (quasi-extinction threshold)
    lambda[t,i] <- N[t+1,17,i]/N[t,17,i]           #Annual population growth rate
  } #t
  lambda[lambda[,i]==0,i]<-NA                   #Replaces any 0s with NAs
  crash[i] <- ifelse(!is.na(lambda[n.yrs,i]), n.yrs+1, min(which(is.na(lambda[,i]))))          #Function to get time of crash (lambda=0 or NA) in each simulation
  
  mean.lambda.iter[i]<- ifelse(crash[i]>15,mean(lambda[16:(crash[i]-1),i],na.rm=T),next)  #Mean population growth rate per iteration
  persist[i]<-ifelse(N[n.yrs+1,18,i]>0,1,0)                     #Indicator (0 or 1) if population persisted at year 50
} #i

mean.lambda <- mean(mean.lambda.iter,na.rm=T)       #Mean population growth rate
perc.change.iter <- sapply(1:iter, function(i) ifelse(N[n.yrs+1,17,i]==0,-1,(N[n.yrs+1,17,i]-N[1,17,i])/N[1,17,i]))      #Percentage change in abundance from beginning of study
perc.change <- mean(perc.change.iter)      #Mean percentage change in abundance from beginning of study
persist[is.na(persist)] <- 0
persistence.prob <- mean(persist,na.rm=T)    #Mean persistence probability per strategy
#### End PVA Model

# Extract adult breeding abundance over time for 10 random iterations (for supplemental figure)
N.scen.dat_nate <- data.frame(N[,17,sample(1:iter,10,replace=F)])
N.scen.dat_nate$year <- 1:nrow(N.scen.dat_nate)
N.scen.dat_nate <- melt(N.scen.dat_nate,id="year")
colnames(N.scen.dat_nate)[2:3] <- c("iter","N")
N.scen.dat_nate$iter <- factor(rep(1:10,each=max(N.scen.dat_nate$year)))
N.scen.dat_nate$pond <- "Nate's Pond"

#Generate summary output at step 45 in all simulations
craw.out.nat <- data.frame(pond="Nates",iter=1:iter,N.breeder=N[45,17,],s.tad=phi.tad.DD[45,],lambda=lambda[45,])
write.csv(craw.out.nat,"craw.out.nat.imm_2.1.22.csv")

#-----------------------#
### CREATE RESULTS TABLE - NATES ###
#-----------------------#
summary.fun = function(data) {  df = data.frame(site="Nates")
df$mean = mean(data, na.rm=T)
df$median = median(data, na.rm=T)
df$low95 = quantile(data, prob=0.025, na.rm=T)
df$up95 = quantile(data, prob=0.975, na.rm=T)
df$low90 = quantile(data, prob=0.05, na.rm=T)
df$up90 = quantile(data, prob=0.9, na.rm=T)
return(df)
}

scen.lambda <- summary.fun(mean.lambda.iter)
scen.lambda[is.na(scen.lambda)] <- 0
scen.time.to.extinct <- summary.fun(crash)
crash.test <- crash
scen.abundance <- summary.fun(N[n.yrs+1,17,])   # Adult breeder abundance at end of model time frame
scen.percchange <- summary.fun(perc.change.iter)

craw.out.summary.nat <- data.frame(Site=c("Nates"),N0.ad.b=N.AB.0,persist=persistence.prob,
                               time.to.ext.mean=scen.time.to.extinct$mean, time.to.ext.low90=scen.time.to.extinct$low90, time.to.ext.up90=scen.time.to.extinct$up90,
                               N.adult.median=scen.abundance$median,N.adult.low90=scen.abundance$low90,N.adult.up90=scen.abundance$up90,
                               lambda.median=scen.lambda$median,lambda.low90=scen.lambda$low90,lamda.up90=scen.lambda$up90,
                               percchange.med=scen.percchange$median,percchange.low90=scen.percchange$low90,percchange.up90=scen.percchange$up90)
#-----------------------#
#### END NATES MODEL ####
#-----------------------#

#-----------------------#
## BEGIN CATTAIL MODEL ##
#-----------------------#

vitals <- vitals.orig[which(vitals.orig$Site=="Cattail"),]

# Extract mean vital rates from input data table
phi.emb <- vitals$Mean[vitals$Rate=="phi.emb"]
phi.tad <- vitals$Mean[vitals$Rate=="phi.tad"]
phi.juv <- vitals$Mean[vitals$Rate=="phi.juv"]
phi.adult.1 <- vitals$Mean[vitals$Rate=="phi.adult.1"]
phi.adult.2 <- vitals$Mean[vitals$Rate=="phi.adult.2"]
p <- vitals$Mean[vitals$Rate=="repro.p"]
q <- vitals$Mean[vitals$Rate=="repro.q"]
repro.age.2.f <- vitals$Mean[vitals$Rate=="repro.age.2.f"]
repro.age.3.f <- vitals$Mean[vitals$Rate=="repro.age.3.f"]
repro.age.4.f <- vitals$Mean[vitals$Rate=="repro.age.4.f"]
clutch <- vitals$Mean[vitals$Rate=="clutch"]
prop.clutch.f <- vitals$Mean[vitals$Rate=="prop.clutch.f"]
gam.BNB <- vitals$Mean[vitals$Rate=="gam.BNB"]
gam.NBNB <- vitals$Mean[vitals$Rate=="gam.NBNB"]

# Extract SD vital rates from input data table
phi.tad.SD <- vitals$SD[vitals$Rate=="phi.tad"]
phi.juv.SD <- vitals$SD[vitals$Rate=="phi.juv"]
phi.adult.1.SD <- vitals$SD[vitals$Rate=="phi.adult.1"]
phi.adult.2.SD <- vitals$SD[vitals$Rate=="phi.adult.2"]
p.SD <- vitals$SD[vitals$Rate=="repro.p"]
q.SD <- vitals$SD[vitals$Rate=="repro.q"]
clutch.SD <- vitals$SD[vitals$Rate=="clutch"]
gam.BNB.SD <- vitals$SD[vitals$Rate=="gam.BNB"]
gam.NBNB.SD <- vitals$SD[vitals$Rate=="gam.NBNB"]

#Calculate stable age distribution for initial population sizes
fert <- clutch*prop.clutch.f
gam <- 1-gam.BNB
s1 <- phi.emb*phi.tad*(phi.juv^0.75)
s2 <- phi.juv*(1-p)
s3 <- phi.juv*(1-gam)*p
s4 <- phi.juv*gam*p
s5 <- phi.juv*(1-q)
s6 <- phi.juv*(1-gam)*q
s7 <- phi.juv*gam*q
s8 <- phi.adult.1*(1-gam)
s9 <- phi.adult.1*gam
s10 <- phi.juv*(1-gam)
s11 <- phi.juv*gam
s12 <- phi.adult.1*(1-gam)
s13 <- phi.adult.1*gam
s14 <- phi.adult.2*(1-gam)
s15 <- phi.adult.2*gam

fert2 <- fert*s4
fert3 <- fert*s7
fert4 <- fert*s9
fert5 <- fert*s11
fert6 <- fert*s13
fert7 <- fert*s15
# Matrix columns:
# 1) year 1: embryo, metamorph, juveniles
# 2) year 2: s2 = juveniles, not mature, s3 = juveniles, mature, not breed, s4 = juveniles, mature, breed
# 3) year 3: s5 = juveniles, not mature, s6 = juveniles, mature, not breed, s7 = juveniles, mature, breed, s8 = adults, not breed, s9 = adults, breed
# 4) year 4: s10 = juveniles, mature, not breed, s11 = juveniles, mature, breed, s12 = 1st adults, not breed, s13 = 1st adults, breed, s14 = 2nd adults, not breed, s15 = 2nd adults, breed

Lefkovitch <- c(0,     fert2,   fert3,  fert4, fert4,  fert5,  fert6,  fert6,  fert7, fert7,  fert6,  fert6, fert7,  fert7,  fert7, fert7,
                s1,    0,       0,      0,     0,      0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0, 
                0,     s2,      0,      0,     0,      0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     s3,      0,      0,     0,      0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,        
                0,     s4,      0,      0,     0,      0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       s5,     0,     0,      0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       s6,     0,     0,      0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       s7,     0,     0,      0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       0,      s8,    s8,     0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       0,      s9,    s9,     0,      0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       0,      0,     0,      s10,    0,      0,      0,     0,      0,      0,     0,      0,      0,     0, 
                0,     0,       0,      0,     0,      s11,    0,      0,      0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       0,      0,     0,      0,      s12,    s12,    0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       0,      0,     0,      0,      s13,    s13,    0,     0,      0,      0,     0,      0,      0,     0,
                0,     0,       0,      0,     0,      0,      0,      0,      s14,   s14,    s12,    s12,   s14,    s14,    s14,   s14, 
                0,     0,       0,      0,     0,      0,      0,      0,      s15,   s15,    s13,    s13,   s15,    s15,    s15,   s15)
craw <- matrix(Lefkovitch,nrow=16,byrow=T)
craw_eigen <- eigen(craw)
stable_stage <- craw_eigen$vectors[,1]/sum(craw_eigen$vectors[,1]) 
stable_stage <-as.numeric(stable_stage)

#--------------------------------------------------------------------------------------------------#
##Population Viability Analysis - stage-based model for 50 years and 1000 iterations per scenario###
#--------------------------------------------------------------------------------------------------#
iter <- 1000  # Number of iterations
n.yrs <- 50   # Number of years to project pop for each iteration
breed.prob <- 0.8 # Probability of breeding in a given year
breed.state <- array(dim=c(n.yrs+1,iter))
breed.state[,] <- sapply(1:iter, function (i) rbinom(n.yrs+1,1,breed.prob))

# Set up main pop abundance (N) array
N<-array(dim=c(n.yrs+1,18,iter))   # Main abundance array for 18 classes (16 stage classes described in above matrix + 1 column for total adult breeders + 1 column for total abundance)
dim(N)    #Check if dimensions are right (1=years, 2=stages, 3=iterations, 4=mgmt scenarios)
lambda<-array(dim=c(n.yrs,iter))
crash <- array(dim=c(iter))
mean.lambda.iter<-array(dim=c(iter))
mean.lambda<-array(dim=1)
perc.change.iter<-array(dim=c(iter))
perc.change<-array(dim=1)
persist<-array(dim=c(iter))
persistence<-array(dim=c(iter))
persistence.prob<-array(dim=1)
N.tads<-array(dim=c(n.yrs+1,iter))
N.met<-array(dim=c(n.yrs+1,iter))
N.imms<-array(dim=c(n.yrs+1,iter))
phi.tad.DD<-array(dim=c(n.yrs+1,iter))

#### Population Viability Model
N.AB.0<-36    # Starting adult pop size
Total0 <- round(N.AB.0/sum(stable_stage[c(5,8,10,12,14,16)]),0)
stage.N <- round(stable_stage*Total0,0)
N0 <- c(stage.N,sum(stage.N[c(5,8,10,12,14,16)]),Total0)

# Begin PVA simulations
for (i in 1:iter){
  # Begin forecasting population abundance over 30 years
  N[1,,i] <- N0                                                                        #Initial abundances
  for (t in 1:n.yrs){
    # Incorporate parametric uncertainty around sensitive parameters (initial abundance, terrestrial survival rates)
    phi.S.juv <- rnorm(1,phi.juv,phi.juv.SD)
    phi.S.juv <- ifelse(phi.S.juv<phi.juv-2*phi.juv.SD,phi.juv-2*phi.juv.SD,ifelse(phi.S.juv>phi.juv+2*phi.juv.SD,phi.juv+2*phi.juv.SD,phi.S.juv))  # Bounded juvenile survival to +/- 2 SDs
    phi.S.adult.2 <- rnorm(1,phi.adult.2,phi.adult.2.SD)
    phi.S.adult.2 <- ifelse(phi.S.adult.2>0.85,0.85,ifelse(phi.S.adult.2<0.2,0.2,phi.S.adult.2))
    phi.S.adult.1 <- 0.9089*phi.S.adult.2 - 0.1924
    phi.S.adult.1 <- ifelse(phi.S.adult.1>0.7,0.7,ifelse(phi.S.adult.1<0.15,0.15,phi.S.adult.1))
    gam.BNB.S <- rnorm(1,gam.BNB,gam.BNB.SD)
    gam.NBNB.S <- rnorm(1,gam.NBNB,gam.NBNB.SD)
    
    N.tads[t+1,i] <- rbinom(1,N[t,1,i],phi.emb)                         # Eggs surviving to tadpoles
    #phi.tad.DD[t+1,i] <- rbeta(1,0.731538027,3250.548582)              # No density dependence at Cattail
    phi.tad.DD[t+1,i] <- rbeta(1,0.899,3000)                           # No density dependence at Cattail - new alpha/beta produce mean and median survival rates estimated from field data (0.0003 and 0.0002, respectively)
    N.met[t+1,i] <- rbinom(1,N.tads[t+1,i],phi.tad.DD[t+1,i])          # Metamorphs
    N[t+1,2,i] <- rbinom(1,N.met[t+1,i],phi.S.juv^0.75)                # Juvs (age 1) = first year survival through embryo, tadpole, and metamorph/juvenile stages
    
    N[t+1,3,i] <- rbinom(1,N[t,2,i],phi.S.juv*(1-p))    #Juvs (age 2) = juvs that didn't mature at age 2 (p) and survived
    N[t+1,4,i] <- rbinom(1,N[t,2,i],phi.S.juv*p*(gam.BNB.S))              #Juvs to ads, nb (age 2) = juvs that matured, survived, and didn't breed
    N[t+1,5,i] <- rbinom(1,N[t,2,i],phi.S.juv*p*(1-gam.BNB.S))                  #Juvs to ads, b (age 2) = juvs that matured, survived, and bred
    
    N[t+1,6,i] <- rbinom(1,N[t,3,i],phi.S.juv*(1-q))    #Juvs (age 3) = juvs that didn't mature at age 3 (q) and survived
    N[t+1,7,i] <- rbinom(1,N[t,3,i],phi.S.juv*q*(gam.BNB.S))              #Juvs to ads, nb (age 3) = juvs that matured, survived, and didn't breed
    N[t+1,8,i] <- rbinom(1,N[t,3,i],phi.S.juv*q*(1-gam.BNB.S))                  #Juvs to ads, b (age 3) = juvs that matured, survived, and bred
    N[t+1,9,i] <- rbinom(1,N[t,4,i],phi.S.adult.1*(gam.NBNB.S)) + rbinom(1,N[t,5,i],phi.S.adult.1*(gam.BNB.S))     #1st yr Ads, nb (age 3) = ads that survived and didn't breed
    N[t+1,10,i] <- rbinom(1,N[t,4,i],phi.S.adult.1*(1-gam.NBNB.S)) + rbinom(1,N[t,5,i],phi.S.adult.1*(1-gam.BNB.S))         #1st yr Ads, b (age 3) = ads that survived and bred
    
    N[t+1,11,i] <- rbinom(1,N[t,6,i],phi.S.juv*(gam.BNB.S))              #Juvs to ads, nb (age 4) = juvs that matured, survived, and didn't breed
    N[t+1,12,i] <- rbinom(1,N[t,6,i],phi.S.juv*(1-gam.BNB.S))                  #Juvs to ads, b (age 4) = juvs that matured, survived, and bred
    imms <- rnbinom(100, size=fit$estimate[1], mu=fit$estimate[2]) # Draw random numbers of breeding adult immigrants
    N.imms[t+1,i] <- imms[min(which(imms<=max(imm.dat)))] # Draw annual number of breeding adult immigrants, restricted to being less than or equal to the max observed number of immigrants
    N[t+1,13,i] <- rbinom(1,sum(N[t,c(7,11),i]),phi.S.adult.1*(gam.NBNB.S)) + rbinom(1,sum(N[t,c(8,12),i]),phi.S.adult.1*(gam.BNB.S))            #1st yr Ads, nb (age 4) = 1st yr ads that survived and didn't breed
    N[t+1,14,i] <- rbinom(1,sum(N[t,c(7,11),i]),phi.S.adult.1*(1-gam.NBNB.S)) + rbinom(1,sum(N[t,c(8,12),i]),phi.S.adult.1*(1-gam.BNB.S))        #1st yr Ads, b (age 4) = 1st yr ads that survived and bred
    N[t+1,14,i] <- N[t+1,14,i] + N.imms[t+1,i] #1st yr Ads, b (age 4) = 1st yr ads that survived and bred + immigrants
    N[t+1,15,i] <- rbinom(1,sum(N[t,c(9,13,15),i]),phi.S.adult.2*(gam.NBNB.S)) + rbinom(1,sum(N[t,c(10,14,16),i]),phi.S.adult.2*(gam.BNB.S))                #2nd+ yr Ads, nb (age 4) = 2nd+ yr ads that survived and didn't breed
    N[t+1,16,i] <- rbinom(1,sum(N[t,c(9,13,15),i]),phi.S.adult.2*(1-gam.NBNB.S)) + rbinom(1,sum(N[t,c(10,14,16),i]),phi.S.adult.2*(1-gam.BNB.S))                   #2nd+ yr Ads, b (age 4) = 2nd+ yr ads that survived and bred
    
    N[t+1,1,i] <- round(sum(N[t,c(5,8,10,12,14,16),i]) * clutch*prop.clutch.f * breed.state[t+1,i],0)  # nobreed is indicator variable simulating different lengths/frequencies of no recruitment years (e.g., droughts)
    
    # Track total abundance, lambda, time to extinction
    N[t+1,17,i] <- sum(N[t+1,c(5,8,10,12,14,16),i])   #Annual adult breeder abundance
    N[t+1,17,i]<-(N[t,17,i]>=2)*N[t+1,17,i]        #Makes N= 0 if adult abundance falls below 2 (quasi-extinction threshold)
    N[t+1,18,i] <- sum(N[t+1,1:16,i])                 #Annual total abundance
    N[t+1,18,i]<-(N[t,17,i]>=2)*N[t+1,18,i]        #Makes N= 0 if adult abundance falls below 2 (quasi-extinction threshold)
    lambda[t,i] <- N[t+1,17,i]/N[t,17,i]           #Annual population growth rate
  } #t
  lambda[lambda[,i]==0,i]<-NA                   #Replaces any 0s with NAs
  crash[i] <- ifelse(!is.na(lambda[n.yrs,i]), n.yrs+1, min(which(is.na(lambda[,i]))))          #Function to get time of crash (lambda=0 or NA) in each simulation
  
  mean.lambda.iter[i]<- ifelse(crash[i]>15,mean(lambda[16:(crash[i]-1),i],na.rm=T),next)  #Mean population growth rate per iteration
  persist[i]<-ifelse(N[n.yrs+1,18,i]>0,1,0)                     #Indicator (0 or 1) if population persisted at year 50
} #i

mean.lambda <- mean(mean.lambda.iter,na.rm=T)       #Mean population growth
perc.change.iter <- sapply(1:iter, function(i) ifelse(N[n.yrs+1,17,i]==0,-1,(N[n.yrs+1,17,i]-N[1,17,i])/N[1,17,i]))      #Percentage change in abundance from beginning of study
perc.change <- mean(perc.change.iter)      #Mean percentage change in abundance from beginning of study
persist[is.na(persist)] <- 0
persistence.prob <- mean(persist,na.rm=T)    #Mean persistence probability per strategy
#### End PVA Model

# Extract adult breeding abundance over time for 10 random iterations (for supplemental figure)
N.scen.dat_cat <- data.frame(N[,17,sample(1:iter,10,replace=F)])
N.scen.dat_cat$year <- 1:nrow(N.scen.dat_cat)
N.scen.dat_cat <- melt(N.scen.dat_cat,id="year")
colnames(N.scen.dat_cat)[2:3] <- c("iter","N")
N.scen.dat_cat$iter <- factor(rep(1:10,each=max(N.scen.dat_cat$year)))
N.scen.dat_cat$pond <- "Cattail Pond"

#Generate summary output at step 45 in all simulations
craw.out.cat <- data.frame(pond="Cattail",iter=1:iter,N.breeder=N[45,17,],s.tad=phi.tad.DD[45,],lambda=lambda[45,])
write.csv(craw.out.cat,"craw.out.cat_2.1.22.csv")

#-----------------------#
### CREATE RESULTS TABLE - CATTAIL ###
#-----------------------#
summary.fun = function(data) {  df = data.frame(site="Cattail")
df$mean = mean(data, na.rm=T)
df$median = median(data, na.rm=T)
df$low95 = quantile(data, prob=0.025, na.rm=T)
df$up95 = quantile(data, prob=0.975, na.rm=T)
df$low90 = quantile(data, prob=0.05, na.rm=T)
df$up90 = quantile(data, prob=0.9, na.rm=T)
return(df)
}

scen.lambda <- summary.fun(mean.lambda.iter)
scen.lambda[is.na(scen.lambda)] <- 0
scen.time.to.extinct <- summary.fun(crash)
crash.test <- crash
scen.abundance <- summary.fun(N[n.yrs+1,17,])   # Adult breeder abundance at end of model time frame
scen.percchange <- summary.fun(perc.change.iter)

craw.out.summary.cat <- data.frame(Site=c("Cattail"),N0.ad.b=N.AB.0,persist=persistence.prob,
                                   time.to.ext.mean=scen.time.to.extinct$mean, time.to.ext.low90=scen.time.to.extinct$low90, time.to.ext.up90=scen.time.to.extinct$up90,
                                   N.adult.median=scen.abundance$median,N.adult.low90=scen.abundance$low90,N.adult.up90=scen.abundance$up90,
                                   lambda.median=scen.lambda$median,lambda.low90=scen.lambda$low90,lamda.up90=scen.lambda$up90,
                                   percchange.med=scen.percchange$median,percchange.low90=scen.percchange$low90,percchange.up90=scen.percchange$up90)

#-----------------------#
### END CATTAIL MODEL ###
#-----------------------#

#-------------------------------#
### CREATE COMBINED RESULTS OUTPUT & SUMMARY TABLE ###
#-------------------------------#

craw.out <- rbind(craw.out.nat,craw.out.cat)
write.csv(craw.out,"craw.mod.output_2.1.22.csv")

plot(craw.out$s.tad,craw.out$lambda)

craw.out.summary <- rbind(craw.out.summary.nat,craw.out.summary.cat)
write.csv(craw.out.summary,"craw.mod.output.summary_2.1.22.csv")

###################################
#### Plot abundance over time #####
###################################
N.scen.dat <- rbind(N.scen.dat_nate,N.scen.dat_cat)

plot.N <- ggplot(N.scen.dat, aes(x=year,y=N, group=interaction(iter,pond))) + 
  geom_line(aes(linetype=pond),size=0.4) +
  geom_hline(yintercept=c(10,50),linetype="dashed",size=0.8) +
  labs(x= "Year of simulation",y="Number of breeding females") + 
  scale_x_continuous(expand=c(0,0),breaks=c(seq(1,n.yrs+1,by=5)),limits = c(1, n.yrs+1), labels=c(seq(0,n.yrs,by=5))) +
  scale_y_continuous(expand=c(0,0),breaks=c(seq(0,140,by=20))) +
  theme_classic() + theme(panel.background = element_rect(color="black"),
                          axis.text.y=element_text(size=11),axis.text.x=element_text(size=11),axis.title=element_text(size=12,face="bold"),
                          legend.position=c(0.1,0.9),legend.title=element_blank(), plot.margin=unit(c(0.3,0.6,0.4,0.4),"cm"))

plot.N

jpeg(file="Fig_Adult abundance over time for 10 iterations_imm_2.1.22.jpg",width = 7, height = 5, units = "in", res=600)
plot.N
dev.off()

# jpeg(file="Fig_Adult abundance over time for scenarios iter_1.7.21.jpg",width = 12, height = 8, units = "in", res=600)
# grid.arrange(arrangeGrob(plot.scen[[1]],plot.scen[[2]],plot.scen[[3]],plot.scen[[4]],
#                          layout_matrix = matrix(c(1,2,3,4), ncol=2, byrow=TRUE)))
# dev.off()

## Graph replationship between tadpole survival and istantaneous population growth

craw.mod.output <- read.csv("craw.mod.output_1.7.21.csv")
pA <- ggplot(data = craw.mod.output, aes(x = lambda)) +
  geom_histogram(aes(fill = pond), alpha =0.5)
pA
