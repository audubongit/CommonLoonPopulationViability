########################################################################################
# Integrated population model (IPM) for Wisconsin common loons, 1995 - 2019
# Code by Sarah Saunders, National Audubon Society, 2019 - 2020
# Data provided by Walter Piper, Chapman University

# Adapted from original scripts by Marc Kéry & Michael Schaub (2016)
# Modified by S. Saunders, 2019 - 2020
########################################################################################

# Load data and libraries
library(jagsUI)

nyears <- 25	  # Number of years in analysis

#Load function to create a m-array based on capture-recapture data (CH)
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

########################################################################
# Capture-recapture data: m-array of juveniles (HY) and adults (AHY)
########################################################################

#First read in capture histories for birds marked as chicks during 1995-2019
CH.J <- read.table("CH_HYmark95to19.txt")

#convert to matrix
CH.J <- data.matrix(CH.J)

#read in capture histories for birds marked as adults during 1995-2019
CH.A <- read.table("CH_AHYmark95to19.txt")

#convert to matrix
CH.A <- data.matrix(CH.A)

#create two m-arrays, one for juveniles and one for adults
cap <- apply(CH.J, 1, sum)
ind <- which(cap >= 2)
CH.J.R <- CH.J[ind,] # Juvenile CH recaptured at least once
CH.J.N <- CH.J[-ind,] # Juvenile CH never recaptured

# Remove first capture
first <- numeric()
for (i in 1:dim(CH.J.R)[1]){
  first[i] <- min(which(CH.J.R[i,]==1))
}
CH.J.R1 <- CH.J.R
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R1[i,first[i]] <- 0
}

# Add grown-up juveniles to adults and create m-array
CH.A.m <- rbind(CH.A, CH.J.R1)
CH.A.marray <- marray(CH.A.m)

# Create CH matrix for juveniles, ignoring subsequent recaptures
second <- numeric()
for (i in 1:dim(CH.J.R1)[1]){
  second[i] <- min(which(CH.J.R1[i,]==1))
}
CH.J.R2 <- matrix(0, nrow = dim(CH.J.R)[1], ncol = dim(CH.J.R)[2])
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R2[i,first[i]] <- 1
  CH.J.R2[i,second[i]] <- 1
}

# Create m-array for these
CH.J.R.marray <- marray(CH.J.R2)

# The last column should show the number of juveniles not recaptured again and should all be zeros, since all of them are released as adults
CH.J.R.marray[,dim(CH.J)[2]] <- 0

# Create the m-array for juveniles never recaptured and add it to the previous m-array
CH.J.N.marray <- marray(CH.J.N)
CH.J.marray <- CH.J.R.marray + CH.J.N.marray

#outputs: CH.A.marray and CH.J.marray
#convert outputs to names of m-arrays used in models

#delete last 2 rows of juvenile m-array (can't recap birds released in last 2 occassions)
marray.j <- CH.J.marray[-c(23:24),]
marray.a <- CH.A.marray

######################################################################################################################################################
# Population count data: all observed breeding pairs (1995 - 2019)
y <-  c(50,57,62,67,74,56,71,70,92,99,105,105,94,111,102,98,109,114,120,118,126,118,116,118,121)

# Productivity data (1995-2019)
J <- c(41,27,38,52,44,34,39,42,56,67,82,58,62,74,51,37,29,65,98,47,76,70,61,76,47) # Number of offspring/fledglings
R <-  c(50,57,62,67,71,52,67,70,92,99,105,105,94,111,102,98,109,114,120,118,126,118,116,118,119) # Number of surveyed pairs - these are virtually 'y'

#####################################################################################################################################################

## Read in scalar info to correct total pair data being fed into model---------------------------------
scalar <- c(0.47,0.63,0.65,0.77,0.71,0.53,0.75,0.75,0.88,1,1.1,1.08,0.82,1,1.07,1.02,1.18,1.11,1.16,1.24,1.43,1.27,1.35,1.36,1.31)
y.scaled <- y/scalar

#read in covariate data (top-supported model data only) ###########################
yeareffects = read.csv("Year_Covs.csv", header=T, sep=',', na.strings=T)
precip = yeareffects$Precip
nao1 = yeareffects$NAO_DJFM_PC1
nao2 = yeareffects$NAO_DJFM_PC2
dev5k = yeareffects$Dev_5km

#standardize effects
mprecip = mean(as.matrix(precip))
sdprecip = sd(as.vector(precip))
rain = as.vector((precip-mprecip)/sdprecip)

mnao1 = mean(as.matrix(nao1))
sdnao1 = sd(as.vector(nao1))
wnao1 = as.vector((nao1-mnao1)/sdnao1) #winter (Dec/Jan/Feb/March NAO)

mnao2 = mean(as.matrix(nao2))
sdnao2 = sd(as.vector(nao2))
wnao2 = as.vector((nao2-mnao2)/sdnao2) #winter (Dec/Jan/Feb/March NAO)

mdev5k = mean(as.matrix(dev5k))
sddev5k = sd(as.vector(dev5k))
dev5km = as.vector((dev5k-mdev5k)/sddev5k) #winter (Dec/Jan/Feb/March NAO)

#############################################################################
# Integrated population model (IPM) for Wisconsin common loons, 1995 - 2019
# Code by Sarah Saunders, National Audubon Society, 2020
# Data provided by Walter Piper, Chapman University
# Adapted from original scripts by Marc Kéry & Michael Schaub (2016)
# Modified by S. Saunders, 2019 - 2020
# See main text for full description of modeling framework
#
# Notations:
# nyears = 25
# marray.j is an m-array of capture histories for individuals first banded as
# chicks during 1995 - 2019
# marray.a is an m-array of capture histories for individuals first banded as
# adults during 1995 - 2019
# y = number of breeding pairs annually (scaled by number of lakes surveyed)
# b = number of surveyed broods annually
# J = number of fledglings observed annually
# wnao1 = winter North Atlantic Oscillation index with a 1-step lag (i.e. NAO
# immediately preceding the breeding season t)
# wnao2 = winter NAO with a 2-step lag (i.e. NAO from ~16 mo prior to
# breeding season t)
# rain = cumulative summer precipitation during year t
# dev5km = proportion of developed land cover within 5km, averaged across
# all study lakes annually
#############################################################################

sink("colo_ipm")
cat("
    model {
    #------------------------------------------------------------------------
    #  Integrated population model
    #  - Stage structured model with 3 stages: 
    #  <1 y olds (fledglings), new 3 y olds, 4+ y olds
    #  - Age at first breeding is fourth year
    #  - 1 and 2 y olds are 'invisible' (don't return to breeding grounds
    #  until age 3)
    #  - Prebreeding census, female-based
    #  - All vital rates are assumed to be time-dependent
    #  - Includes env. stochasticity thru random time effects for all params
    #  - Model includes the following covariates (resulting from model
    #  selection procedure described in the main text):
    #  - Summer rainfall, winter NAO from 2 winters prior, and their
    #  interaction on fecundity 
    #  - Preceding winter NAO and developed land cover within 5 km on adult
    #  survival
    #------------------------------------------------------------------------
    
    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    
    # Initial population sizes
    nadRet ~ dnorm(50, 0.001)I(0,)        # Adults 4+ y old returning

    N0[1] ~ dpois(0.5 * f[1] * NAdRet[1]) # Initial pop size of fledglings
    NAdRet[1] <- round(nadRet)   
    
    for (t in 1:3){
    nadNew[t] ~ dnorm(50, 0.001)I(0,)      # New 3 year olds; needs prior for 
    NAdNew[t] <- round(nadNew[t])          # first 3 years
    Ntot[t] <- NAdRet[t] + NAdNew[t]       # Ntot derived as returners and 
    }#t                                    # new 3 y olds (i.e. all adults)        

    # Mean demographic parameters (on appropriate scale)
    l.mphij ~ dnorm(0, 0.37)               # juvenile survival
    l.mphia ~ dnorm(0, 0.37)               # adult survival
    l.mfec ~ dnorm(0, 0.001)               # productivity
    l.p ~ dnorm(0, 0.37)                   # mean detection probability

    # Priors for beta coefficients
    beta.fec1 ~ dnorm(0, 0.01)
    beta.fec2 ~ dnorm(0, 0.01)
    beta.fec3 ~ dnorm(0, 0.01)
    beta.phia1 ~ dnorm(0, 0.01)
    beta.phia2 ~ dnorm(0, 0.01)

    # Precision of standard deviations of temporal variability
    sig.phij ~ dunif(0, 10)
    tau.phij <- pow(sig.phij, -2)
    sig.phia ~ dunif(0, 10)
    tau.phia <- pow(sig.phia, -2)
    sig.fec ~ dunif(0, 10)
    tau.fec <- pow(sig.fec, -2)
    sig.res ~ dunif(0, 10)
    tau.res <- pow(sig.res, -2)
    
    sig.obs ~ dunif(0.5, 50)   # residual variance
    tau.obs <- pow(sig.obs, -2)
    

    # Distribution of error terms (Bounded to help with convergence)
    for (t in 1:(nyears-1)){
    epsilon.phia[t] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.res[t] ~ dnorm(0, tau.res)T(-5,5)
    }
    
    for (t in 1:nyears){
    epsilon.fec[t] ~ dnorm(0, tau.fec)T(-5,5)
    }
    
    for (t in 1:(nyears-3)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)T(-5,5)
    }
    
    #---------------------------------------------
    # 2. Constrain parameters (temp variability)
    #---------------------------------------------
    for (t in 1:(nyears-1)){
    logit(phia[t]) <- l.mphia + beta.phia1*wnao1[t] + beta.phia2*dev5km[t] +
    epsilon.phia[t]  # epsilon.phia is random temporal effect for env. stoch.                           
    logit(p[t]) <- l.p + epsilon.res[t]       
    }
    
    for (t in 1:nyears){
    log(f[t]) <- l.mfec + beta.fec1*rain[t] + beta.fec2*wnao2[t] +
    beta.fec3*rain[t]*wnao2[t] + epsilon.fec[t]                  
    }
    
    for (t in 1:(nyears-3)){
    logit(phij[t]) <- l.mphij + epsilon.phij[t]         
    }
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    mphij <- exp(l.mphij)/(1+exp(l.mphij))        # Mean juvenile survival 
    mphia <- exp(l.mphia)/(1+exp(l.mphia))        # Mean adult survival 
    mfec <- exp(l.mfec)                           # Mean productivity

    # Population growth rate (total adult breeders [4+ y olds])
    for (t in 1:(nyears-1)){
    lambda[t] <- NAdRet[t+1] / (NAdRet[t] + 0.0001)   
    logla[t] <- log(lambda[t])
    }
    mlam <- exp((1/(nyears-1))*sum(logla[1:(nyears-1)]))   # Geo mean all yrs
    mlam.five <- exp((1/(nyears-20))*sum(logla[20:(nyears-1)]))    # Last 5 y
    mlam.ten <- exp((1/(nyears-15))*sum(logla[15:(nyears-1)]))    # Last 10 y
  
    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    # 4.1. Likelihood for population count data (state-space model)
    # 4.1.1 System process
    for (k in 2:nyears){
    mean1[k] <- 0.5 * f[k] * NAdRet[k]        # only NadRet can breed
    N0[k] ~ dpois(mean1[k])                   # Number of fledged females
    NAdRet[k] ~ dbin(phia[k-1], Ntot[k-1])
    }#k
    
    for (t in 4:nyears){  
    NAdNew[t] ~ dbin(phij[t-3],N0[t-3])                 
    Ntot[t] <- NAdNew[t] + NAdRet[t]        
    }#t

    # 4.1.2 Observation process
    for (t in 1:nyears){      
    y[t] ~ dnorm(NAdRet[t], tau.obs)          # Only NAdRet are counted
    }
    
    # 4.2 Likelihood for capture-recapture data: CJS model
    # Multinomial likelihood
    for (t in 1:(nyears-3)){                  
    marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t]) 
    }
    
    for (t in 1:(nyears-1)){
    marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])  
    }
    
    # m-array cell probabilities for juveniles
    for (t in 1:(nyears-3)){    
    # Main diagonal
    pr.j[t,(t+2)] <- phij[t]*p[t+2]
    # Above main diagonal
    for (j in (t+3):(nyears-1)){   
    pr.j[t,j] <- phij[t]*prod(phia[(t+3):j])*prod(q[(t+2):(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t+1)){  
    pr.j[t,j] <- 0   
    } #j
    # Last column
    pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
    } #t
    
    # m-array cell probabilities for adults  
    for (t in 1:(nyears-1)){      
    q[t] <- 1-p[t]     
    # Main diagonal
    pr.a[t,t] <- phia[t]*p[t]
    # above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    # Last column
    pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
    } #t
    
    # 4.3. Likelihood for productivity data: Poisson regression
    for (t in 1:nyears){   
    J[t] ~ dpois(rho[t])    # number young observed as fledged
    rho[t] <- b[t] * f[t]   # number monitored broods (b) and fecundity
    }
    }
    ",fill = TRUE)
sink()

###################################################################
# Bundle data
jags.data <- list(nyears = nyears, marray.j = marray.j, marray.a = marray.a, y = y.scaled, J = J, b = b, r.j = rowSums(marray.j), r.a = rowSums(marray.a), wnao1 = wnao1, rain = rain, wnao2 = wnao2, dev5km = dev5km) 

# Initial values
inits <- function(){list(l.mphij = rnorm(1, 0.2, 0.5), l.mphia = rnorm(1, 0.2, 0.5), l.mfec = rnorm(1, 0.2, 0.5), l.p = rnorm(1, 0.2, 1), sig.phij = runif(1, 0.1, 10), sig.phia = runif(1, 0.1, 10), sig.fec = runif(1, 0.1, 10), beta.fec1 = rnorm(1, 0, 1), beta.fec2 = rnorm(1, 0, 1), beta.fec3 = rnorm(1, 0, 1), beta.phia1 = rnorm(1, 0, 1), beta.phia2 = rnorm(1, 0, 1), nadRet = round(runif(1, 1, 50), 0), nadNew = round(runif(3, 5, 50), 0), sigma.obs = runif(1, 0, 1))}

# Parameters monitored
parameters <- c("phij", "phia", "f", "lambda", 
                "mphij", "mphia","mfec", "mlam", "mlam.five", "mlam.ten",
                "beta.fec1", "beta.fec2", "beta.fec3", 
                "beta.phia1", "beta.phia2",
                "l.mphij", "l.mphia","l.mfec",
                "p", "sig.phij", "sig.phia", "sig.fec", "sig.obs", 
                "N0", "NAdRet","NAdNew", "Ntot") 

# MCMC settings
ni <- 200000    
nt <- 5
nb <- 150000
nc <- 3

# Call JAGS from R
colo.out <- jags(jags.data, inits, parameters, "colo_ipm", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)

###############################################################################################################################
## Bayesian population viability analysis with top-supported model covariates
#############################################################################################################################

# future covariate values provided as data ------------------------------------------------------------
yeareffects = read.csv("Year_covs_Future.csv", header=T, sep=',', na.strings=T)
precip = yeareffects$Precip
dev5 = yeareffects$Dev_5km

#standardize effects
mprecip.fut = mean(as.matrix(precip[1:25]))
sdprecip.fut = sd(as.vector(precip[1:25]))
rain.fut = as.vector((precip-mprecip.fut)/sdprecip.fut)

mdev5.fut = mean(as.matrix(dev5[1:25]))
sddev5.fut = sd(as.vector(dev5[1:25]))
dev5k.fut = as.vector((dev5-mdev5.fut)/sddev5.fut)

###################################################
## Bootstrap method to obtain temporal variation
## in NAO values over time
##################################################

## Read in data
nao_all <- read.csv("Historic_NAO.csv", stringsAsFactors = FALSE)

## Bootstrap method for pulling 10 random values of historical NAO
## to calculate SD [temporal variability of NAO]
library(boot)

sdstats <- function(data, indices){
  dat <- data[indices,] #allows boot to select sample
  s1 <- sample(dat[,2], size=10, replace=FALSE)
  sd.s <- sd(s1)
  return(sd.s)
}
# bootstrapping wtih 1000 reps
set.seed(1234)
results <- boot(data=nao_all, statistic=sdstats,
                R=100000)

# view results
results
mnao.sd <- mean(results$t) #mean of bootstrapped samples: 1.12
sdnao.sd <- sd(results$t) #sd of bootstrapped samples: 0.24

#-----------------------------------------------------------------------------------------
# Model example below has prior/hyperprior structure for NAO with long-term NAO mean (no trend)
# Example structure for future scenarios A and D
#-----------------------------------------------------------------------------------------
sink("threestage.ipm.bpva.shared")
cat("
    model {
    
    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    
    # Initial population sizes
    nadRet ~ dnorm(50, 0.001)I(0,)        # Adults 4+ y old returning

    N0[1] ~ dpois(0.5 * f[1] * NAdRet[1]) # Initial pop size of fledglings based on fecundity est. and prior for breeding adults
    NAdRet[1] <- round(nadRet)   
    
    for (t in 1:3){
    nadNew[t] ~ dnorm(50, 0.001)I(0,)      # New 3 year olds 
    NAdNew[t] <- round(nadNew[t])  
    Ntot[t] <- NAdRet[t] + NAdNew[t]       # Ntot derived as returners (4+ y olds) and new 3 y olds (all adults)
    }#t                                                

    # Mean demographic parameters (on appropriate scale)
    l.mphij ~ dnorm(0, 0.37)               # juvenile survival 
    l.mphia ~ dnorm(0, 0.37)               # adult survival [>= 3 y old] (all ages that consistently migrate)
    l.mfec ~ dnorm(0, 0.001)               # productivity 
    l.p ~ dnorm(0, 0.37)                   # mean detection probability
    beta.fec1 ~ dnorm(0, 0.01)
    beta.fec2 ~ dnorm(0, 0.01)
    beta.fec3 ~ dnorm(0, 0.01)
    beta.phia1 ~ dnorm(0, 0.01)
    beta.phia2 ~ dnorm(0, 0.01)
    
    # Precision of standard deviations of temporal variability
    sig.phij ~ dunif(0, 10)
    tau.phij <- pow(sig.phij, -2)
    sig.phia ~ dunif(0, 10)
    tau.phia <- pow(sig.phia, -2)
    sig.fec ~ dunif(0, 10)
    tau.fec <- pow(sig.fec, -2)
    sig.res ~ dunif(0, 10)
    tau.res <- pow(sig.res, -2)
    
    sig.obs ~ dunif(0.5, 50)   # residual variance
    tau.obs <- pow(sig.obs, -2)
    
    sig.sig <- sdnao.sd  #sdnao.sd is 0.24 from bootstrap results 
    tau.sig <- pow(sig.sig, -2)
    sig.nao ~ dnorm(mnao.sd, tau.sig) #mnao.sd is 1.12 from bootstrap results
    tau.nao <- pow(sig.nao, -2)
    
    # Distribution of error terms (Bounded to help with convergence)
    for (t in 1:(nyears-1+K)){
    epsilon.phia[t] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.res[t] ~ dnorm(0, tau.res)T(-5,5)
    }
    
    for (t in 1:(nyears+K)){
    epsilon.fec[t] ~ dnorm(0, tau.fec)T(-5,5)  
    }
    
    for (t in 1:(nyears-3+K)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)T(-5,5)
    }
    
    for (t in nyears:(nyears-1+K)){         #pull future NAO values from prior
    wnao.pre[t] ~ dnorm(0.19, tau.nao)         
    wnao1.fut[t] <- (wnao.pre[t] - mnao1)/sdnao1
    wnao2.fut[t] <- (wnao.pre[t] - mnao2)/sdnao2
    }
    
    #---------------------------------------------
    # 2. Constrain parameters (temp variability)
    #---------------------------------------------
    # Past: same as best model
    for (t in 1:(nyears-1)){
    logit(phia[t]) <- l.mphia + beta.phia1*wnao1[t] + beta.phia2*dev5km[t] + epsilon.phia[t]              # Adult apparent survival (3+ y old survival)
    }
    
    for (t in 1:nyears){
    log(f[t]) <- l.mfec + beta.fec1*rain[t] + beta.fec2*wnao2[t] + beta.fec3*rain[t]*wnao2[t] + epsilon.fec[t]        # Productivity - additional estimate for final year
    }
    
    for (t in 1:(nyears-3+K)){                                             # same for past + future time periods
    logit(phij[t]) <- l.mphij + epsilon.phij[t]                            # Juv. apparent survival (from fledge to 3 y old)
    }
    
    for (t in 1:(nyears-1+K)){                       # same for past + future time periods
    logit(p[t]) <- l.p + epsilon.res[t]              # Recapture/resight probability (added epsilon here for temp variability)
    }
    
    # Future: use rain and dev5km future values; priors for NAO on fec and phia
    for (t in nyears:(nyears-1+K)){
    logit(phia[t]) <- l.mphia + beta.phia1*wnao1.fut[t] + beta.phia2*dev5k.fut[t] + epsilon.phia[t]              # Adult apparent survival (3+ y old survival)
    }
    
    for (t in (nyears+1):(nyears+K)){
    log(f[t]) <- l.mfec + beta.fec1*rain.fut[t] + beta.fec2*wnao2.fut[t-1] + beta.fec3*rain.fut[t]*wnao2.fut[t-1] + epsilon.fec[t]        # Productivity - additional estimate for final year
    }
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    mphij <- exp(l.mphij)/(1+exp(l.mphij))        # Mean juvenile survival probability
    mphia <- exp(l.mphia)/(1+exp(l.mphia))        # Mean adult survival probability
    mfec <- exp(l.mfec)                           # Mean productivity

    # Population growth rate (total adult breeders [4+ y olds])
    for (t in 1:(nyears-1+K)){
    lambda[t] <- NAdRet[t+1] / (NAdRet[t] + 0.0001)   
    logla[t] <- log(lambda[t])
    }
    mlam.hist <- exp((1/(nyears-1))*sum(logla[1:(nyears-1)]))          # Geometric mean growth (1995-2019)
    mlam.five <- exp((1/(nyears-20))*sum(logla[20:(nyears-1)]))        # Geometric growth (2015-2019)
    mlam.ten <- exp((1/(nyears-15))*sum(logla[15:(nyears-1)]))         # Geometric growth (2010-2019)
    mlam.tot <- exp((1/(nyears-1+K))*sum(logla[1:(nyears-1+K)]))       # Geometric mean (1995 - 2029)
    mlam.fut <- exp((1/(K-1))*sum(logla[nyears:(nyears-1+K)]))   	     # Geometric mean for 2020-2029
  
    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    # 4.1. Likelihood for population count data (state-space model)
    # 4.1.1 System process
    
    for (k in 2:(nyears+K)){
    mean1[k] <- 0.5 * f[k] * NAdRet[k]        # ONLY NadRet, i.e. 4+ y olds, are potential breeders
    N0[k] ~ dpois(mean1[k])                   # Number of fledged females 
    NAdRet[k] ~ dbin(phia[k-1], Ntot[k-1])
    }#k
    
    for (t in 4:(nyears+K)){  
    NAdNew[t] ~ dbin(phij[t-3],N0[t-3])                 # phij is survival & transitioning to age 3
    Ntot[t] <- NAdNew[t] + NAdRet[t]        
    }#t

    # 4.1.2 Observation process
    for (t in 1:nyears){      
    y[t] ~ dnorm(NAdRet[t], tau.obs)          
    }
    
    # 4.2 Likelihood for capture-recapture data: CJS model
    # Multinomial likelihood
    for (t in 1:(nyears-3)){                   
    marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t]) 
    }
    
    for (t in 1:(nyears-1)){
    marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])  
    }
    
    # m-array cell probabilities for juveniles
    for (t in 1:(nyears-3)){    
    # Main diagonal
    pr.j[t,(t+2)] <- phij[t]*p[t+2]
    # Above main diagonal
    for (j in (t+3):(nyears-1)){   
    pr.j[t,j] <- phij[t]*prod(phia[(t+3):j])*prod(q[(t+2):(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t+1)){  
    pr.j[t,j] <- 0   
    } #j
    # Last column
    pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
    } #t
    
    # m-array cell probabilities for adults  
    for (t in 1:(nyears-1)){      
    q[t] <- 1-p[t]     
    # Main diagonal
    pr.a[t,t] <- phia[t]*p[t]
    # above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    # Last column
    pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
    } #t
    
    # 4.3. Likelihood for productivity data: Poisson regression
    for (t in 1:nyears){   
    J[t] ~ dpois(rho[t])    
    rho[t] <- R[t] * f[t]   
    }
    }
    ",fill = TRUE)
sink()

###################################################################
# Bundle data
K <- 10 #number of years in prediction window
jags.data <- list(nyears = nyears, marray.j = marray.j, marray.a = marray.a, y = y.scaled, J = J, R = R, K = K,
                  r.j = rowSums(marray.j), r.a = rowSums(marray.a), 
                  dev5km = dev5km, rain = rain, wnao1 = wnao1, wnao2 = wnao2,
                  dev5k.fut = dev5k.fut, rain.fut = rain.fut,
                  mnao.sd = mnao.sd, sdnao.sd = sdnao.sd,
                  mnao1=mnao1, sdnao1=sdnao1, mnao2=mnao2, sdnao2=sdnao2) 

# Initial values
inits <- function(){list(l.mphij = rnorm(1, 0.2, 0.5), l.mphia = rnorm(1, 0.2, 0.5), l.mfec = rnorm(1, 0.2, 0.5), l.p = rnorm(1, 0.2, 1), 
                         sig.phij = runif(1, 0.1, 10), sig.phia = runif(1, 0.1, 10), sig.fec = runif(1, 0.1, 10),  
                         beta.fec1 = rnorm(1, 0, 1), beta.fec2 = rnorm(1, 0, 1), beta.fec3 = rnorm(1, 0, 1), 
                         beta.phia1 = rnorm(1, 0, 1), beta.phia2 = rnorm(1, 0, 1),
                         nadRet = round(runif(1, 1, 50), 0), 
                         nadNew = round(runif(3, 5, 50), 0), 
                         sigma.obs = runif(1, 0, 1))}    

# Parameters monitored
parameters <- c("phij", "phia", "f", "lambda", 
                "mphij", "mphia", "mfec", "mlam.hist", "mlam.five", "mlam.ten", "mlam.tot","mlam.fut",
                "wnao1.fut", "wnao2.fut", "wnao.pre",
                "beta.fec1", "beta.fec2", "beta.fec3",
                "beta.phia1", "beta.phia2",
                "p", "sig.phij", "sig.phia", "sig.fec", "sig.obs", 
                "N0", "NAdRet","NAdNew", "Ntot") 

# MCMC settings
ni <- 200000    
nt <- 5
nb <- 150000
nc <- 3

# Call JAGS from R
threestage.ipm.scaleddat95.bpva <- jags(jags.data, inits, parameters, "threestage.ipm.bpva.shared", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)
print(threestage.ipm.scaleddat95.bpva, digits=3)

#-----------------------------------------------------------------------------------------
# Model example below has prior/hyperprior structure for NAO with positive/negative trend of NAO
# Example structure used for future scenarios B, C, E, F
#-----------------------------------------------------------------------------------------

sink("threestage.ipm.bpva.shared.trend")
cat("
    model {
    
    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    
    # Initial population sizes
    nadRet ~ dnorm(50, 0.001)I(0,)        

    N0[1] ~ dpois(0.5 * f[1] * NAdRet[1]) 
    NAdRet[1] <- round(nadRet)   
    
    for (t in 1:3){
    nadNew[t] ~ dnorm(50, 0.001)I(0,)       
    NAdNew[t] <- round(nadNew[t])  
    Ntot[t] <- NAdRet[t] + NAdNew[t]       
    }#t                                                

    # Mean demographic parameters (on appropriate scale)
    l.mphij ~ dnorm(0, 0.37)               
    l.mphia ~ dnorm(0, 0.37)               
    l.mfec ~ dnorm(0, 0.001)               
    l.p ~ dnorm(0, 0.37)                   
    beta.fec1 ~ dnorm(0, 0.01)
    beta.fec2 ~ dnorm(0, 0.01)
    beta.fec3 ~ dnorm(0, 0.01)
    beta.phia1 ~ dnorm(0, 0.01)
    beta.phia2 ~ dnorm(0, 0.01)
    
    # Precision of standard deviations of temporal variability
    sig.phij ~ dunif(0, 10)
    tau.phij <- pow(sig.phij, -2)
    sig.phia ~ dunif(0, 10)
    tau.phia <- pow(sig.phia, -2)
    sig.fec ~ dunif(0, 10)
    tau.fec <- pow(sig.fec, -2)
    sig.res ~ dunif(0, 10)
    tau.res <- pow(sig.res, -2)
    
    sig.obs ~ dunif(0.5, 50)   # residual variance
    tau.obs <- pow(sig.obs, -2)
    
    sig.sig <- sdnao.sd  #sdnao.sd is 0.24 from bootstrap results 
    tau.sig <- pow(sig.sig, -2)
    sig.nao ~ dnorm(mnao.sd, tau.sig) #mnao.sd is 1.12 from bootstrap results
    tau.nao <- pow(sig.nao, -2)
    
    # Distribution of error terms (Bounded to help with convergence)
    for (t in 1:(nyears-1+K)){
    epsilon.phia[t] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.res[t] ~ dnorm(0, tau.res)T(-5,5)
    }
    
    for (t in 1:(nyears+K)){
    epsilon.fec[t] ~ dnorm(0, tau.fec)T(-5,5)  
    }
    
    for (t in 1:(nyears-3+K)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)T(-5,5)
    }
    
    mu.nao[24] <- 0.19                     #use long-term (0.08) or short-term (0.19) mean here
    for (t in nyears:(nyears-1+K)){         #pull future NAO values from prior
    mu.nao[t] <- mu.nao[t-1] - 0.24         # use small (0.02) LT trend or large ST trend (0.24)
    wnao.pre[t] ~ dnorm(mu.nao[t], tau.nao)         
    wnao1.fut[t] <- (wnao.pre[t] - mnao1)/sdnao1
    wnao2.fut[t] <- (wnao.pre[t] - mnao2)/sdnao2
    }
    
    #---------------------------------------------
    # 2. Constrain parameters (temp variability)
    #---------------------------------------------
    # Past: same as best model
    for (t in 1:(nyears-1)){
    logit(phia[t]) <- l.mphia + beta.phia1*wnao1[t] + beta.phia2*dev5km[t] + epsilon.phia[t]             
    }
    
    for (t in 1:nyears){
    log(f[t]) <- l.mfec + beta.fec1*rain[t] + beta.fec2*wnao2[t] + beta.fec3*rain[t]*wnao2[t] + epsilon.fec[t]        
    }
    
    for (t in 1:(nyears-3+K)){                                             
    logit(phij[t]) <- l.mphij + epsilon.phij[t]                            
    }
    
    for (t in 1:(nyears-1+K)){                       
    logit(p[t]) <- l.p + epsilon.res[t]              
    }
    
    # Future: use rain and dev5km future values; priors for NAO on fec and phia
    for (t in nyears:(nyears-1+K)){
    logit(phia[t]) <- l.mphia + beta.phia1*wnao1.fut[t] + beta.phia2*dev5k.fut[t] + epsilon.phia[t]              
    }
    
    for (t in (nyears+1):(nyears+K)){
    log(f[t]) <- l.mfec + beta.fec1*rain.fut[t] + beta.fec2*wnao2.fut[t-1] + beta.fec3*rain.fut[t]*wnao2.fut[t-1] + epsilon.fec[t]        # Productivity - additional estimate for final year
    }
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    mphij <- exp(l.mphij)/(1+exp(l.mphij))        # Mean juvenile survival probability
    mphia <- exp(l.mphia)/(1+exp(l.mphia))        # Mean adult survival probability
    mfec <- exp(l.mfec)                           # Mean productivity

    # Population growth rate (total adult breeders [4+ y olds])
    for (t in 1:(nyears-1+K)){
    lambda[t] <- NAdRet[t+1] / (NAdRet[t] + 0.0001)   
    logla[t] <- log(lambda[t])
    }
    mlam.hist <- exp((1/(nyears-1))*sum(logla[1:(nyears-1)]))          # Geometric mean growth (1995-2019)
    mlam.five <- exp((1/(nyears-20))*sum(logla[20:(nyears-1)]))        # Geometric growth (2015-2019)
    mlam.ten <- exp((1/(nyears-15))*sum(logla[15:(nyears-1)]))         # Geometric growth (2010-2019)
    mlam.tot <- exp((1/(nyears-1+K))*sum(logla[1:(nyears-1+K)]))       # Geometric mean (1995 - 2029)
    mlam.fut <- exp((1/(K-1))*sum(logla[nyears:(nyears-1+K)]))   	     # Geometric mean for 2020-2029
  
    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    # 4.1. Likelihood for population count data (state-space model)
    # 4.1.1 System process
    
    for (k in 2:(nyears+K)){
    mean1[k] <- 0.5 * f[k] * NAdRet[k]        
    N0[k] ~ dpois(mean1[k])                   
    NAdRet[k] ~ dbin(phia[k-1], Ntot[k-1])
    }#k
    
    for (t in 4:(nyears+K)){  
    NAdNew[t] ~ dbin(phij[t-3],N0[t-3])                 
    Ntot[t] <- NAdNew[t] + NAdRet[t]        
    }#t

    # 4.1.2 Observation process
    for (t in 1:nyears){      
    y[t] ~ dnorm(NAdRet[t], tau.obs)          
    }
    
    # 4.2 Likelihood for capture-recapture data: CJS model
    # Multinomial likelihood
    for (t in 1:(nyears-3)){                  
    marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t]) 
    }
    
    for (t in 1:(nyears-1)){
    marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])  
    }
    
    # m-array cell probabilities for juveniles
    for (t in 1:(nyears-3)){    
    # Main diagonal
    pr.j[t,(t+2)] <- phij[t]*p[t+2]
    # Above main diagonal
    for (j in (t+3):(nyears-1)){   
    pr.j[t,j] <- phij[t]*prod(phia[(t+3):j])*prod(q[(t+2):(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t+1)){  
    pr.j[t,j] <- 0   
    } #j
    # Last column
    pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
    } #t
    
    # m-array cell probabilities for adults  
    for (t in 1:(nyears-1)){      
    q[t] <- 1-p[t]     
    # Main diagonal
    pr.a[t,t] <- phia[t]*p[t]
    # above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    # Last column
    pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
    } #t
    
    # 4.3. Likelihood for productivity data: Poisson regression
    for (t in 1:nyears){   
    J[t] ~ dpois(rho[t])    
    rho[t] <- R[t] * f[t]   
    }
    }
    ",fill = TRUE)
sink()

###################################################################
# Bundle data
K <- 10
jags.data <- list(nyears = nyears, marray.j = marray.j, marray.a = marray.a, y = y.scaled, J = J, R = R, K = K,
                  r.j = rowSums(marray.j), r.a = rowSums(marray.a), 
                  dev5km = dev5km, rain = rain, wnao1 = wnao1, wnao2 = wnao2,
                  dev5k.fut = dev5k.fut, rain.fut = rain.fut,
                  mnao.sd = mnao.sd, sdnao.sd = sdnao.sd,
                  mnao1=mnao1, sdnao1=sdnao1, mnao2=mnao2, sdnao2=sdnao2) 

# Initial values
inits <- function(){list(l.mphij = rnorm(1, 0.2, 0.5), l.mphia = rnorm(1, 0.2, 0.5), l.mfec = rnorm(1, 0.2, 0.5), l.p = rnorm(1, 0.2, 1), 
                         sig.phij = runif(1, 0.1, 10), sig.phia = runif(1, 0.1, 10), sig.fec = runif(1, 0.1, 10),  
                         beta.fec1 = rnorm(1, 0, 1), beta.fec2 = rnorm(1, 0, 1), beta.fec3 = rnorm(1, 0, 1), 
                         beta.phia1 = rnorm(1, 0, 1), beta.phia2 = rnorm(1, 0, 1),
                         nadRet = round(runif(1, 1, 50), 0), 
                         nadNew = round(runif(3, 5, 50), 0), 
                         sigma.obs = runif(1, 0, 1))}    

# Parameters monitored
parameters <- c("phij", "phia", "f", "lambda", 
                "mphij", "mphia", "mfec", "mlam.hist", "mlam.five", "mlam.ten", "mlam.tot","mlam.fut",
                "wnao1.fut", "wnao2.fut", "wnao.pre",
                "beta.fec1", "beta.fec2", "beta.fec3",
                "beta.phia1", "beta.phia2",
                "p", "sig.phij", "sig.phia", "sig.fec", "sig.obs", 
                "N0", "NAdRet","NAdNew", "Ntot") 

# MCMC settings
ni <- 200000    
nt <- 5
nb <- 150000
nc <- 3

# Call JAGS from R
threestage.ipm.scaleddat95.bpva <- jags(jags.data, inits, parameters, "threestage.ipm.bpva.shared.trend", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)
print(threestage.ipm.scaleddat95.bpva, digits=3)

#-------------------------------------------------------------------------------------------------------------------
# Model example below has prior/hyperprior structure for NAO with long-term NAO mean (no trend)
# Has adjustment to lower juvenile survival. Example structure used for scenarios A and D under low juv survival
#------------------------------------------------------------------------------------------------------------------
sink("threestage.ipm.bpva.shared.lowphij")
cat("
    model {

    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    
    # Initial population sizes
    nadRet ~ dnorm(50, 0.001)I(0,)        

    N0[1] ~ dpois(0.5 * f[1] * NAdRet[1]) 
    NAdRet[1] <- round(nadRet)   
    
    for (t in 1:3){
    nadNew[t] ~ dnorm(50, 0.001)I(0,)      
    NAdNew[t] <- round(nadNew[t])  
    Ntot[t] <- NAdRet[t] + NAdNew[t]       
    }#t                                                

    # Mean demographic parameters (on appropriate scale)
    l.mphij ~ dnorm(0, 0.37)               
    l.mphia ~ dnorm(0, 0.37)               
    l.mfec ~ dnorm(0, 0.001)               
    l.p ~ dnorm(0, 0.37)                   
    beta.fec1 ~ dnorm(0, 0.01)
    beta.fec2 ~ dnorm(0, 0.01)
    beta.fec3 ~ dnorm(0, 0.01)
    beta.phia1 ~ dnorm(0, 0.01)
    beta.phia2 ~ dnorm(0, 0.01)
    
    # Precision of standard deviations of temporal variability
    sig.phij ~ dunif(0, 10)
    tau.phij <- pow(sig.phij, -2)
    sig.phia ~ dunif(0, 10)
    tau.phia <- pow(sig.phia, -2)
    sig.fec ~ dunif(0, 10)
    tau.fec <- pow(sig.fec, -2)
    sig.res ~ dunif(0, 10)
    tau.res <- pow(sig.res, -2)
    
    sig.obs ~ dunif(0.5, 50)   # residual variance
    tau.obs <- pow(sig.obs, -2)
    
    sig.sig <- sdnao.sd  #sdnao.sd is 0.24 from bootstrap results 
    tau.sig <- pow(sig.sig, -2)
    sig.nao ~ dnorm(mnao.sd, tau.sig) #mnao.sd is 1.12 from bootstrap results
    tau.nao <- pow(sig.nao, -2)
    
    # Distribution of error terms (Bounded to help with convergence)
    for (t in 1:(nyears-1+K)){
    epsilon.phia[t] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.res[t] ~ dnorm(0, tau.res)T(-5,5)
    }
    
    for (t in 1:(nyears+K)){
    epsilon.fec[t] ~ dnorm(0, tau.fec)T(-5,5)  
    }
    
    for (t in 1:(nyears-3+K)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)T(-5,5)
    }
    
    for (t in nyears:(nyears-1+K)){         #pull future NAO values from prior
    wnao.pre[t] ~ dnorm(0.19, tau.nao)         
    wnao1.fut[t] <- (wnao.pre[t] - mnao1)/sdnao1
    wnao2.fut[t] <- (wnao.pre[t] - mnao2)/sdnao2
    }
    
    #---------------------------------------------
    # 2. Constrain parameters (temp variability)
    #---------------------------------------------
    # Past: same as best model
    for (t in 1:(nyears-1)){
    logit(phia[t]) <- l.mphia + beta.phia1*wnao1[t] + beta.phia2*dev5km[t] + epsilon.phia[t]              
    }
    
    for (t in 1:nyears){
    log(f[t]) <- l.mfec + beta.fec1*rain[t] + beta.fec2*wnao2[t] + beta.fec3*rain[t]*wnao2[t] + epsilon.fec[t]        
    }
    
    for (t in 1:(nyears-3)){                                             
    logit(phij[t]) <- l.mphij + epsilon.phij[t]                            
    }
    
    for (t in 1:(nyears-1+K)){                       
    logit(p[t]) <- l.p + epsilon.res[t]              
    }
    
    # Future: use rain and dev5km future values; priors for NAO on fec and phia
    for (t in nyears:(nyears-1+K)){
    logit(phia[t]) <- l.mphia + beta.phia1*wnao1.fut[t] + beta.phia2*dev5k.fut[t] + epsilon.phia[t]              
    }
    
    for (t in (nyears+1):(nyears+K)){
    log(f[t]) <- l.mfec + beta.fec1*rain.fut[t] + beta.fec2*wnao2.fut[t-1] + beta.fec3*rain.fut[t]*wnao2.fut[t-1] + epsilon.fec[t]        # Productivity - additional estimate for final year
    }
    
    #define new mean phij for future
    phij.fut <- mean(phij[18:22])                              #use mean of last 5 years
    
    for (t in (nyears-2):(nyears-3+K)){                        #future time period uses juv survival that is lower than long-term mean                                
    logit(phij[t]) <- logit(phij.fut) + epsilon.phij[t]        
    }                                                                     
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    mphij <- exp(l.mphij)/(1+exp(l.mphij))        # Mean juvenile survival probability
    mphia <- exp(l.mphia)/(1+exp(l.mphia))        # Mean adult survival probability
    mfec <- exp(l.mfec)                           # Mean productivity

    # Population growth rate (total adult breeders [4+ y olds])
    for (t in 1:(nyears-1+K)){
    lambda[t] <- NAdRet[t+1] / (NAdRet[t] + 0.0001)   
    logla[t] <- log(lambda[t])
    }
    mlam.hist <- exp((1/(nyears-1))*sum(logla[1:(nyears-1)]))          # Geometric mean growth (1995-2019)
    mlam.five <- exp((1/(nyears-20))*sum(logla[20:(nyears-1)]))        # Geometric growth (2015-2019)
    mlam.ten <- exp((1/(nyears-15))*sum(logla[15:(nyears-1)]))         # Geometric growth (2010-2019)
    mlam.tot <- exp((1/(nyears-1+K))*sum(logla[1:(nyears-1+K)]))       # Geometric mean (1995 - 2029)
    mlam.fut <- exp((1/(K-1))*sum(logla[nyears:(nyears-1+K)]))   	     # Geometric mean for 2020-2029
  
    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    # 4.1. Likelihood for population count data (state-space model)
    # 4.1.1 System process
    
    for (k in 2:(nyears+K)){
    mean1[k] <- 0.5 * f[k] * NAdRet[k]        
    N0[k] ~ dpois(mean1[k])                   
    NAdRet[k] ~ dbin(phia[k-1], Ntot[k-1])
    }#k
    
    for (t in 4:(nyears+K)){  
    NAdNew[t] ~ dbin(phij[t-3],N0[t-3])                 
    Ntot[t] <- NAdNew[t] + NAdRet[t]        
    }#t

    # 4.1.2 Observation process
    for (t in 1:nyears){      
    y[t] ~ dnorm(NAdRet[t], tau.obs)          
    }
    
    # 4.2 Likelihood for capture-recapture data: CJS model
    # Multinomial likelihood
    for (t in 1:(nyears-3)){                  
    marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t]) 
    }
    
    for (t in 1:(nyears-1)){
    marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])  
    }
    
    # m-array cell probabilities for juveniles
    for (t in 1:(nyears-3)){    
    # Main diagonal
    pr.j[t,(t+2)] <- phij[t]*p[t+2]
    # Above main diagonal
    for (j in (t+3):(nyears-1)){   
    pr.j[t,j] <- phij[t]*prod(phia[(t+3):j])*prod(q[(t+2):(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t+1)){  
    pr.j[t,j] <- 0   
    } #j
    # Last column
    pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
    } #t
    
    # m-array cell probabilities for adults  
    for (t in 1:(nyears-1)){      
    q[t] <- 1-p[t]     
    # Main diagonal
    pr.a[t,t] <- phia[t]*p[t]
    # above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    # Last column
    pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
    } #t
    
    # 4.3. Likelihood for productivity data: Poisson regression
    for (t in 1:nyears){   
    J[t] ~ dpois(rho[t])    
    rho[t] <- R[t] * f[t]   
    }
    }
    ",fill = TRUE)
sink()

###################################################################
# Bundle data
K <- 10
jags.data <- list(nyears = nyears, marray.j = marray.j, marray.a = marray.a, y = y.scaled, J = J, R = R, K = K,
                  r.j = rowSums(marray.j), r.a = rowSums(marray.a), 
                  dev5km = dev5km, rain = rain, wnao1 = wnao1, wnao2 = wnao2,
                  dev5k.fut = dev5k.fut, rain.fut = rain.fut,
                  mnao.sd = mnao.sd, sdnao.sd = sdnao.sd,
                  mnao1=mnao1, sdnao1=sdnao1, mnao2=mnao2, sdnao2=sdnao2) 

# Initial values
inits <- function(){list(l.mphij = rnorm(1, 0.2, 0.5), l.mphia = rnorm(1, 0.2, 0.5), l.mfec = rnorm(1, 0.2, 0.5), l.p = rnorm(1, 0.2, 1), 
                         sig.phij = runif(1, 0.1, 10), sig.phia = runif(1, 0.1, 10), sig.fec = runif(1, 0.1, 10),  
                         beta.fec1 = rnorm(1, 0, 1), beta.fec2 = rnorm(1, 0, 1), beta.fec3 = rnorm(1, 0, 1), 
                         beta.phia1 = rnorm(1, 0, 1), beta.phia2 = rnorm(1, 0, 1),
                         nadRet = round(runif(1, 1, 50), 0), 
                         nadNew = round(runif(3, 5, 50), 0), 
                         sigma.obs = runif(1, 0, 1))}    

# Parameters monitored
parameters <- c("phij", "phia", "f", "lambda", 
                "mphij", "mphia", "mfec", "mlam.hist", "mlam.five", "mlam.ten", "mlam.tot","mlam.fut",
                "wnao1.fut", "wnao2.fut", "wnao.pre",
                "beta.fec1", "beta.fec2", "beta.fec3",
                "beta.phia1", "beta.phia2",
                "p", "sig.phij", "sig.phia", "sig.fec", "sig.obs", 
                "N0", "NAdRet","NAdNew", "Ntot") 

# MCMC settings
ni <- 200000    
nt <- 5
nb <- 150000
nc <- 3

# Call JAGS from R
threestage.ipm.scaleddat95.bpva <- jags(jags.data, inits, parameters, "threestage.ipm.bpva.shared.lowphij", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)
print(threestage.ipm.scaleddat95.bpva, digits=3)

#-------------------------------------------------------------------------------------------------------------
# Model example below has prior/hyperprior structure for NAO with positive/negative trend of NAO
# Has adjustment to low juvenile survival. Example structure for scenarios B, C, E, F under low juv survival
#--------------------------------------------------------------------------------------------------------------

sink("threestage.ipm.bpva.shared.trendlowphij")
cat("
    model {

    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    
    # Initial population sizes
    nadRet ~ dnorm(50, 0.001)I(0,)        

    N0[1] ~ dpois(0.5 * f[1] * NAdRet[1]) 
    NAdRet[1] <- round(nadRet)   
    
    for (t in 1:3){
    nadNew[t] ~ dnorm(50, 0.001)I(0,)      
    NAdNew[t] <- round(nadNew[t])  
    Ntot[t] <- NAdRet[t] + NAdNew[t]       
    }#t                                                

    # Mean demographic parameters (on appropriate scale)
    l.mphij ~ dnorm(0, 0.37)               
    l.mphia ~ dnorm(0, 0.37)               
    l.mfec ~ dnorm(0, 0.001)               
    l.p ~ dnorm(0, 0.37)                   
    beta.fec1 ~ dnorm(0, 0.01)
    beta.fec2 ~ dnorm(0, 0.01)
    beta.fec3 ~ dnorm(0, 0.01)
    beta.phia1 ~ dnorm(0, 0.01)
    beta.phia2 ~ dnorm(0, 0.01)
    
    # Precision of standard deviations of temporal variability
    sig.phij ~ dunif(0, 10)
    tau.phij <- pow(sig.phij, -2)
    sig.phia ~ dunif(0, 10)
    tau.phia <- pow(sig.phia, -2)
    sig.fec ~ dunif(0, 10)
    tau.fec <- pow(sig.fec, -2)
    sig.res ~ dunif(0, 10)
    tau.res <- pow(sig.res, -2)
    
    sig.obs ~ dunif(0.5, 50)   # residual variance
    tau.obs <- pow(sig.obs, -2)
    
    sig.sig <- sdnao.sd  #sdnao.sd is 0.24 from bootstrap results 
    tau.sig <- pow(sig.sig, -2)
    sig.nao ~ dnorm(mnao.sd, tau.sig) #mnao.sd is 1.12 from bootstrap results
    tau.nao <- pow(sig.nao, -2)
    
    # Distribution of error terms (Bounded to help with convergence)
    for (t in 1:(nyears-1+K)){
    epsilon.phia[t] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.res[t] ~ dnorm(0, tau.res)T(-5,5)
    }
    
    for (t in 1:(nyears+K)){
    epsilon.fec[t] ~ dnorm(0, tau.fec)T(-5,5)  
    }
    
    for (t in 1:(nyears-3+K)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)T(-5,5)
    }
    
    mu.nao[24] <- 0.19                     #use long-term (0.08) or short-term (0.19) mean here
    for (t in nyears:(nyears-1+K)){         #pull future NAO values from prior
    mu.nao[t] <- mu.nao[t-1] - 0.24         # use small (0.02) LT trend or large ST trend (0.24)
    wnao.pre[t] ~ dnorm(mu.nao[t], tau.nao)         
    wnao1.fut[t] <- (wnao.pre[t] - mnao1)/sdnao1
    wnao2.fut[t] <- (wnao.pre[t] - mnao2)/sdnao2
    }
    
    #---------------------------------------------
    # 2. Constrain parameters (temp variability)
    #---------------------------------------------
    # Past: same as best model
    for (t in 1:(nyears-1)){
    logit(phia[t]) <- l.mphia + beta.phia1*wnao1[t] + beta.phia2*dev5km[t] + epsilon.phia[t]              
    }
    
    for (t in 1:nyears){
    log(f[t]) <- l.mfec + beta.fec1*rain[t] + beta.fec2*wnao2[t] + beta.fec3*rain[t]*wnao2[t] + epsilon.fec[t]        
    }
    
    for (t in 1:(nyears-3)){                                             
    logit(phij[t]) <- l.mphij + epsilon.phij[t]                            
    }
    
    for (t in 1:(nyears-1+K)){                       
    logit(p[t]) <- l.p + epsilon.res[t]              
    }
    
    # Future: use rain and dev5km future values; priors for NAO on fec and phia
    for (t in nyears:(nyears-1+K)){
    logit(phia[t]) <- l.mphia + beta.phia1*wnao1.fut[t] + beta.phia2*dev5k.fut[t] + epsilon.phia[t]              
    }
    
    for (t in (nyears+1):(nyears+K)){
    log(f[t]) <- l.mfec + beta.fec1*rain.fut[t] + beta.fec2*wnao2.fut[t-1] + beta.fec3*rain.fut[t]*wnao2.fut[t-1] + epsilon.fec[t]        # Productivity - additional estimate for final year
    }
    
    #define new mean phij for future
    phij.fut <- mean(phij[18:22])                              #use mean of last 5 years
    
    for (t in (nyears-2):(nyears-3+K)){                        #future time period uses juv survival that is lower than long-term mean                                
    logit(phij[t]) <- logit(phij.fut) + epsilon.phij[t]        
    }  
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    mphij <- exp(l.mphij)/(1+exp(l.mphij))        # Mean juvenile survival probability
    mphia <- exp(l.mphia)/(1+exp(l.mphia))        # Mean adult survival probability
    mfec <- exp(l.mfec)                           # Mean productivity

    # Population growth rate (total adult breeders [4+ y olds])
    for (t in 1:(nyears-1+K)){
    lambda[t] <- NAdRet[t+1] / (NAdRet[t] + 0.0001)   
    logla[t] <- log(lambda[t])
    }
    mlam.hist <- exp((1/(nyears-1))*sum(logla[1:(nyears-1)]))          # Geometric mean growth (1995-2019)
    mlam.five <- exp((1/(nyears-20))*sum(logla[20:(nyears-1)]))        # Geometric growth (2015-2019)
    mlam.ten <- exp((1/(nyears-15))*sum(logla[15:(nyears-1)]))         # Geometric growth (2010-2019)
    mlam.tot <- exp((1/(nyears-1+K))*sum(logla[1:(nyears-1+K)]))       # Geometric mean (1995 - 2029)
    mlam.fut <- exp((1/(K-1))*sum(logla[nyears:(nyears-1+K)]))   	     # Geometric mean for 2020-2029
  
    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    # 4.1. Likelihood for population count data (state-space model)
    # 4.1.1 System process
    
    for (k in 2:(nyears+K)){
    mean1[k] <- 0.5 * f[k] * NAdRet[k]        
    N0[k] ~ dpois(mean1[k])                   
    NAdRet[k] ~ dbin(phia[k-1], Ntot[k-1])
    }#k
    
    for (t in 4:(nyears+K)){  
    NAdNew[t] ~ dbin(phij[t-3],N0[t-3])                 
    Ntot[t] <- NAdNew[t] + NAdRet[t]        
    }#t

    # 4.1.2 Observation process
    for (t in 1:nyears){      
    y[t] ~ dnorm(NAdRet[t], tau.obs)          
    }
    
    # 4.2 Likelihood for capture-recapture data: CJS model
    # Multinomial likelihood
    for (t in 1:(nyears-3)){                  
    marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t]) 
    }
    
    for (t in 1:(nyears-1)){
    marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])  
    }
    
    # m-array cell probabilities for juveniles
    for (t in 1:(nyears-3)){    
    # Main diagonal
    pr.j[t,(t+2)] <- phij[t]*p[t+2]
    # Above main diagonal
    for (j in (t+3):(nyears-1)){   
    pr.j[t,j] <- phij[t]*prod(phia[(t+3):j])*prod(q[(t+2):(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t+1)){  
    pr.j[t,j] <- 0   
    } #j
    # Last column
    pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
    } #t
    
    # m-array cell probabilities for adults  
    for (t in 1:(nyears-1)){      
    q[t] <- 1-p[t]     
    # Main diagonal
    pr.a[t,t] <- phia[t]*p[t]
    # above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    # Last column
    pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
    } #t
    
    # 4.3. Likelihood for productivity data: Poisson regression
    for (t in 1:nyears){   
    J[t] ~ dpois(rho[t])    
    rho[t] <- R[t] * f[t]   
    }
    }
    ",fill = TRUE)
sink()

###################################################################
# Bundle data
K <- 10
jags.data <- list(nyears = nyears, marray.j = marray.j, marray.a = marray.a, y = y.scaled, J = J, R = R, K = K,
                  r.j = rowSums(marray.j), r.a = rowSums(marray.a), 
                  dev5km = dev5km, rain = rain, wnao1 = wnao1, wnao2 = wnao2,
                  dev5k.fut = dev5k.fut, rain.fut = rain.fut,
                  mnao.sd = mnao.sd, sdnao.sd = sdnao.sd,
                  mnao1=mnao1, sdnao1=sdnao1, mnao2=mnao2, sdnao2=sdnao2) 

# Initial values
inits <- function(){list(l.mphij = rnorm(1, 0.2, 0.5), l.mphia = rnorm(1, 0.2, 0.5), l.mfec = rnorm(1, 0.2, 0.5), l.p = rnorm(1, 0.2, 1), 
                         sig.phij = runif(1, 0.1, 10), sig.phia = runif(1, 0.1, 10), sig.fec = runif(1, 0.1, 10),  
                         beta.fec1 = rnorm(1, 0, 1), beta.fec2 = rnorm(1, 0, 1), beta.fec3 = rnorm(1, 0, 1), 
                         beta.phia1 = rnorm(1, 0, 1), beta.phia2 = rnorm(1, 0, 1),
                         nadRet = round(runif(1, 1, 50), 0), 
                         nadNew = round(runif(3, 5, 50), 0), 
                         sigma.obs = runif(1, 0, 1))}    

# Parameters monitored
parameters <- c("phij", "phia", "f", "lambda", 
                "mphij", "mphia", "mfec", "mlam.hist", "mlam.five", "mlam.ten", "mlam.tot","mlam.fut",
                "wnao1.fut", "wnao2.fut", "wnao.pre",
                "beta.fec1", "beta.fec2", "beta.fec3",
                "beta.phia1", "beta.phia2",
                "p", "sig.phij", "sig.phia", "sig.fec", "sig.obs", 
                "N0", "NAdRet","NAdNew", "Ntot") 

# MCMC settings
ni <- 200000    
nt <- 5
nb <- 150000
nc <- 3

# Call JAGS from R
threestage.ipm.scaleddat95.bpva <- jags(jags.data, inits, parameters, "threestage.ipm.bpva.shared.trendlowphij", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)
print(threestage.ipm.scaleddat95.bpva, digits=3)

