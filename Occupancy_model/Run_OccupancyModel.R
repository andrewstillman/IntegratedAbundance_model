#### ------------------------------------------------------------------------------ ####
#### Specification of Bayesian temporal auto-logistic occupancy model for Black-backed Woodpeckers
# Random effect for fireID
# Habitat type modeled as a random intercept
# Fire-level intercept with covariates
# 13 covars on psi
# 3 covars on p
# Using R version 4.1.0 and JAGS version 4.3.0
#### ------------------------------------------------------------------------------- ####

library(R2jags)
library(dclone)

#### Load data and organize for the model
## Note: data are archived as .csv files to facilitate long-term storage
## Be sure that the current working directory is set to the "Occupancy_model" folder
occ.files <- list.files(pattern="*.csv")
for (i in 1:length(occ.files)) assign(occ.files[i], read.csv(occ.files[i]))

X.array <- simplify2array(list(as.matrix(X.array.yr1.csv),as.matrix(X.array.yr2.csv),
                            as.matrix(X.array.yr3.csv),as.matrix(X.array.yr4.csv),
                            as.matrix(X.array.yr5.csv),as.matrix(X.array.yr6.csv),
                            as.matrix(X.array.yr7.csv),as.matrix(X.array.yr8.csv),
                            as.matrix(X.array.yr9.csv),as.matrix(X.array.yr10.csv)))

X.naive <- apply(X.array, c(1, 3), max, na.rm = T) # max detection for each site
X.naive[X.naive == -Inf] <- 0
X.naive[is.na(X.naive)] <- 0
head(X.naive)

mod.data = list(
  X.array = X.array,
  X.naive = X.naive,
  ef = detection_covars.csv$ef,
  itype = detection_covars.csv$itype,
  S.elev = site_covars.csv$S.elev,
  S.lat = site_covars.csv$S.lat,
  bs30 = site_covars.csv$bs30,
  S.bs100 = site_covars.csv$S.bs100,
  S.d.patch = site_covars.csv$S.d.patch,
  S.pyro500 = site_covars.csv$S.pyro500,
  S.precc100 = site_covars.csv$S.precc100,
  S.fir100 = site_covars.csv$S.fir100,
  whr = as.factor(site_covars.csv$whr),
  S.fire.age = S.fire.age.csv,
  S.jday = S.jday.csv,
  site.code = site_covars.csv$site.code,
  fireID =as.factor(site_covars.csv$fireID),
  S.fire.size = fire_covars.csv$S.fire.size,
  fire.season = fire_covars.csv$fire.season
)

## Clean up the workspace
rm(list=setdiff(ls(), "mod.data"))
attach(mod.data)

## Create indicator variable for high severity
length(which(bs30 > 25)) #1258 sites
Ind <- rep(0, length(site.code))
Ind[which(bs30 > 25)] <- 1


#### -------------------------------- JAGS MODEL ------------------------------------------ ####

mod.function <- function(){
  
  #### PRIORS
  
  ## Priors on logit-linear model coefficients
  a.day ~ dnorm(0,0.2)    # coefficients on detection
  a.type ~ dnorm(0,0.2) 
  a.ef ~ dnorm(0,0.2)
  b.age ~ dnorm(0,0.2)    # coefficients on occupancy
  b.agesq ~ dnorm(0,0.2)
  b.elev ~ dnorm(0,0.2) 
  b.elevsq ~ dnorm(0,0.2)
  b.lat ~ dnorm(0,0.2)
  b.sev ~ dnorm(0,0.2)
  b.pyro ~ dnorm(0,0.2)
  b.patch ~ dnorm(0,0.2)
  b.precc ~ dnorm(0,0.2)
  b.fir ~ dnorm(0,0.2)
  b.elevlat ~ dnorm(0,0.2) # interaction effects
  b.agesev ~ dnorm(0,0.2)
  b.agefir ~ dnorm(0,0.2)
  phi ~ dnorm(0,0.2)       # auto-logistic component
  f0 ~ dnorm(0,0.001)      # fire-level coefficients on random intercept (linear model)
  f.size ~ dnorm(0,0.001)    
  f.season ~ dnorm(0,0.001)
  
  ## Prior on detection model intercept      
  p0 ~ dunif(0.01, 0.99)     # expected value in interval (0,1)
  a0 <- log(p0/(1-p0))       # logit-transform  
  
  ## Random effect for WHR-type
  whr.sigma ~ dunif(0.01, 3)    # sd between 0.01 and 3 on logit scale
  whr.tau <- 1/(whr.sigma * whr.sigma)
  for (m in 1:n.whr) {
    b.whr[m] ~ dnorm(0, whr.tau)	   # WHR-specific BLUP 
  }
  
  ## Hierarchical fire-level intercept
  fire.sigma ~ dunif(0.01, 3)    # sd between 0.01 and 3 on logit scale
  fire.tau <- 1/(fire.sigma * fire.sigma)
  # Hierarchical structure that governs fire-level intercept
  for (i in 1:n.fires){
    mu.fire[i] <- f0 + f.size*Size[i] + f.season*Season[i] # linear predictors on intercept
    b.fire[i] ~ dnorm(mu.fire[i], fire.tau)                # randomness added in
  }
  b0 <- mean(b.fire)  # monitor average intercept (across all fires) for predictions and plotting

  #### MODEL: Loop through all sites
  for(j in 1:n.sites){
    
    #### Model for first survey t = 1
    ## State process 
    ## MT: changed notation for b.whr. I need to change this below as well
    z[j,1] ~ dbern(psi[j,1])
    logit(psi[j,1]) <- b.fire[Fire[j]] + b.whr[WHR[j]] + 
                       b.age*Age[j,1] + b.agesq*(Age[j,1])^2 + 
                       b.elev*Elev[j] + b.elevsq*(Elev[j])^2 + b.lat*Lat[j] +
                       b.sev*Sev100[j] + b.pyro*Pyro500[j] + b.patch*Patch[j]*Ind[j] +
                       b.precc*Precc[j] + b.elevlat*(Elev[j]*Lat[j]) +
                       b.fir*Fir[j] + b.agefir*(Age[j,1]*Fir[j]) +
                       b.agesev*(Age[j,1]*Sev100[j])
    
    ## Observation process
    for(k in 1:n.surveys){
      y[j,k,1] ~ dbern(z.p[j,k,1])
      z.p[j,k,1] <- z[j,1]*p[j,k,1] 
      logit(p[j,k,1]) <- a0 + a.day*Day[j,1] + a.type*Type[k] + a.ef*Ef[k]
    }
    
    #### Model for subsequent surveys t > 1
    ## State process  
    for(t in 2:n.years){
      z[j,t] ~ dbern(psi[j,t])
      logit(psi[j,t]) <- b.fire[Fire[j]] + b.whr[WHR[j]] + 
                         b.age*Age[j,t] + b.agesq*(Age[j,t])^2 + 
                         b.elev*Elev[j] + b.elevsq*(Elev[j])^2 + b.lat*Lat[j] +
                         b.sev*Sev100[j] + b.pyro*Pyro500[j] + b.patch*Patch[j]*Ind[j] +
                         b.precc*Precc[j] + b.elevlat*(Elev[j]*Lat[j]) +
                         b.fir*Fir[j] + b.agefir*(Age[j,t]*Fir[j]) +
                         b.agesev*(Age[j,t]*Sev100[j]) + phi*(z[j,t-1])
      
      ## Observation process
      for(k in 1:n.surveys){
        y[j,k,t] ~ dbern(z.p[j,k,t])
        z.p[j,k,t] <- z[j,t]*p[j,k,t] 
        logit(p[j,k,t]) <- a0 + a.day*Day[j,1] + a.type*Type[k] + a.ef*Ef[k]
      }
    }
  }
}


#### -------------------------------- RUN JAGS ------------------------------------------ ####

## Data for JAGS
jags.data <- list(y = X.array,
                  Day = S.jday,
                  Type = itype,
                  Ef = ef,
                  Age = S.fire.age,
                  Elev = S.elev,
                  Lat = S.lat,
                  Sev100 = S.bs100,
                  Pyro500 = S.pyro500,
                  Patch = S.d.patch,
                  Precc = S.precc100,
                  WHR = as.numeric(whr),
                  Fir = S.fir100,
                  Ind = Ind,
                  Fire = as.numeric(fireID),
                  Size = S.fire.size,
                  Season = fire.season,
                  n.sites = dim(X.array)[1],
                  n.surveys = dim(X.array)[2],
                  n.years = dim(X.array)[3],
                  n.fires = length(unique(fireID)),
                  n.whr = length(unique(whr)))


params.save <- c("a0","a.day","a.type","a.ef","b.whr","b.age","b.agesq","b.elev","b.elevsq",
                 "b.lat","b.sev","b.pyro","b.patch","b.precc","b.elevlat","b.fir","b.agefir","b.agesev",
                 "f.size","f.season","phi","b0")  

## Initial values
inits <- function() {
  list(z = X.naive)  
}

## Values for MCMC
nc <- 3
n.adapt <- 1000
n.burn <- 30000
n.iter <- 50000
n.thin <- 100   # yields 1500 posterior samples

## Parallelize 
cl <- makePSOCKcluster(nc) # nc is number of cores (1 for each chain)
tmp <- clusterEvalQ(cl, library(dclone))  # Check that `dclone` is loaded on each of cl's workers. 
parLoadModule(cl, "glm")   # load the JAGS module 'glm' on each worker
parListModules(cl)         # make sure previous line worked.


#### Run the model
start.time<-Sys.time()
occ.run <- jags.parfit(cl, data=jags.data, params=params.save, model=mod.function, inits=inits, 
                       n.chains=nc, n.adapt=n.adapt, n.update=n.burn, n.iter=n.iter, thin=n.thin)
stopCluster(cl)    # close out the cluster.
elapsed.time = difftime(Sys.time(), start.time, units='mins')
elapsed.time

