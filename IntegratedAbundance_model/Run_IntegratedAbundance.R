#### ---------------------------------------------------------------------------------------- ####
#### Integrated space-use and occupancy model to predict Black-backed Woodpecker abundance and density
#
# 1. Organize data: Load data, create + scale variables, write function.
# 2. Predict abundance
# 3. Predict density 
# 4. Calculate uncertainty in density estimates
#
# This code was originally written for the RShiny app hosted by the Institute for Bird Populations (IBP).
# We recommend that users run the model using the app (available at the IBP website). 
# The code posted here provides an example using the Reading Fire (Lassen National Forest, 2012) 
# and pre-made spatial layers covering the fire area. 
#
# Code written by Andrew Stillman and Morgan Tingley. See Tingley et al. (Methods in Ecology and Evolution, 2016)
# for details on the model components and validation statistics. This version of the model includes 
# updated occupancy and snag model components based on additional surveys and new findings.  

#### ---------------------------------------------------------------------------------------- ####

#### Load required libraries and set folders
library(raster)
library(sf)

mod.folder <- getwd()
GIS.folder <- "../ReadingFire_data/"


#### ------------------------------------------- 1 ----------------------------------------- ####


#### ------------------ User inputs
setwd(GIS.folder)
fire_cc <- raster("Reading_cc.tif")
fire_outline <- read_sf("fire_outline.shp")
age <- 1                                  # Number of years post-fire           
input_season <- "January 1 to August 15"  # Ignition date
                                          # Two options: 1) "January 1 to August 15" 2) "August 16 to December 31"
input_size <- 28079                       # Area of fire, in acres.
                       


#### ------------------ Calculations on inputs to produce variables of interest

#### Burn severity. 9x9 cells approximates 100m buffer around each pixel
fire_cc100 <- focal(fire_cc, w = matrix(1/81, nr=9, nc=9))

#### Pyrodiversity - Simpsons diversity of severity bins within 500m
## First, reclassify to 11 categories
reclass_df <- c(-Inf, 0, 1,  
                0, 10, 2, 
                10, 20, 3,
                20, 30, 4,
                30, 40, 5,
                40, 50, 6,
                50, 60, 7,
                60, 70, 8,
                70, 80, 9,
                80, 90, 10,
                90, Inf, 11)
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)
fire_pyro <- reclassify(fire_cc, reclass_m)         

## Function to calculate pyrodiversity 
p <- 1:11 # number of categories
div.simp <- function(r, verbose=TRUE) {
  #
  # Create focal weights matrix.
  mat <- focalWeight(r, d=500, type='circle')
  #
  # Sum over squared proportions in each category
  entropy.r <- calc(r, function(x) 0)  # this is just a raster of 0s
  simp <- function(x) ifelse(x==0, 0, x^2)
  for (i in 1:length(p)) {
    cat("Computing for severity =", i, "... ")
    z <- system.time({
      r1 <- focal(r == i, mat, sum)    # proportion of category i in buffer
      r2 <- calc(r1, simp)             # if prop > 0, then square that cell
      entropy.r <- entropy.r + r2      # add the cell onto existing cells
    })
    cat(z[3], "seconds.\n")
  }
  return (entropy.r)
}

## Run function. Takes a minute
fire_pyro <- div.simp(fire_pyro) # returns D 
fire_pyro <- 1/fire_pyro         # Inverse simpsons diversity


#### Distance to patch edge
## reclassify >25% change cc = medium or high severity
reclass_df <- c(-Inf, 25, 1,  
                25, Inf, NA)
reclass_m <- matrix(reclass_df, ncol = 3, byrow = TRUE)
reclass_m

bi.fire_cc <- reclassify(fire_cc, reclass_m)   # binary raster: High severity is NA  
fire_patch <- distance(bi.fire_cc)             # takes a few seconds


#### Latitude
lat.pts <- as.data.frame(coordinates(fire_cc))
lat.pts <- st_as_sf(lat.pts, coords = c("x","y"), crs=st_crs(fire_cc))
lat.pts <- st_transform(lat.pts, crs = 4269)  # EPSG code for NAD83 lat/long
lat_coords <- st_coordinates(lat.pts)[,2]
fire_lat <- fire_cc
values(fire_lat) <- lat_coords


#### Tidy up the workspace
rm(lat.pts,lat_coords,reclass_m,reclass_df,p,bi.fire_cc)


#### ---------------------- Load pre-made layers

## NOTE: In the app, these layers are automatically produced for you. 
## See app documentation for information on data sources and considerations.

fire_dem <- raster("Reading_dem.tif")

fire_whr_crosswalk <- read.csv("veg_crosswalk.csv")
fire_whr <- raster("Reading_whr.tif")

fire_precc <- raster("Reading_precc.tif")

fire_fir <- raster("Reading_fir.tif")


#### ---------------------- Load model posteriors and scale variables

#### Load data and posteriors from three Bayesian models
## Data archived in .csv format in the IntegratedAbundance_model folder
setwd(mod.folder)
age.raw <- read.csv("fire.age_raw.csv")
fire.covars <- read.csv("fire.covars.csv")
site.covars.raw <- read.csv("site.covars_raw.csv")
snag.data <- read.csv("snag.data.csv")

post_hr <- read.csv("post.hr.csv")
post_psi <- read.csv("post.psi.csv")
post_snag <- read.csv("post.snag.csv")

#### Organize the data for modeling
occ_data <- list(
  cc.raw = site.covars.raw$cc.raw,
  pyro.raw = site.covars.raw$pyro.raw,
  patch.raw = site.covars.raw$patch.raw,
  age.raw = as.matrix(age.raw),
  precc.raw = site.covars.raw$precc.raw,
  elev.raw = site.covars.raw$elev.raw,
  lat.raw = site.covars.raw$lat.raw,
  whrtype = factor(site.covars.raw$whrtype),
  fir.raw = site.covars.raw$fir.raw,              
  firesize.log = fire.covars$firesize.log, 
  fireseason.raw = fire.covars$fireseason.raw,  
  scale.patch = list(mean=139.6914, sd=173.1197)) # data used to scale the d.patch variable
rm(fire.covars, site.covars.raw, age.raw)

snag_data <- as.list(snag.data$values)
names(snag_data) <- snag.data$statistic
rm(snag.data)

post_b.whr <- as.matrix(post_psi[,1:10])


#### Extracting prediction data from raster
## Extract index for pixels within the area for which you want to calculate total abundance
fire_outline_rast <- rasterize(x = fire_outline, y = fire_cc)  # create raster footprint of fire
in_fire           <- (1 - is.na(values(fire_outline_rast)))
pix               <- which(in_fire == 1)                       # identifies which pixels are inside of fire

## Extract values from rasters (for speed)
cc_val      <- values(fire_cc100)[pix]
pyro_val    <- values(fire_pyro)[pix]
patch_val   <- values(fire_patch)[pix]
lat_val     <- values(fire_lat)[pix]
elev_val    <- values(fire_dem)[pix]
precc_val   <- values(fire_precc)[pix]
fir_val     <- values(fire_fir)[pix]
veg_val_ID  <- values(fire_whr)[pix]

## Now it's safe to remove these rasters
rm(fire_cc100, fire_pyro, fire_patch, fire_lat, fire_dem, fire_precc, fire_whr, fire_fir)

#### Scaling of values for snag model
snag_cc_val <- values(fire_cc)[pix]
snag_cc_val <- (snag_cc_val - snag_data$mu.cc) / snag_data$sd.cc
snag_precc_val <- (precc_val - snag_data$mu.precc) / snag_data$sd.precc
snag_elev_val <- (elev_val - snag_data$mu.elev) / snag_data$sd.elev 
snag_lat_val <- (lat_val - snag_data$mu.lat) / snag_data$sd.lat
snag_size_val <- fire_whr_crosswalk$WHRSIZE[match(veg_val_ID, fire_whr_crosswalk$Value)]
# Two lines below force levels to 0 - 5
# X becomes 0. 5 and 6 are combined. Final categories are 0:5.
snag_size_val[which(snag_size_val=="X")] <- 0   # non-forest.
snag_size_val[which(snag_size_val=="6")] <- 5   # Multi-layer = 5
snag_size_val <- as.numeric(snag_size_val)

#### Scaling variables for occupancy model
cc_val    <- (cc_val - mean(occ_data$cc.raw)) / sd(occ_data$cc.raw)
pyro_val  <- (pyro_val - mean(occ_data$pyro.raw)) / sd(occ_data$pyro.raw)
patch_val <- (patch_val - occ_data$scale.patch$mean) / occ_data$scale.patch$sd         
lat_val <- (lat_val - mean(occ_data$lat.raw)) / sd(occ_data$lat.raw)
elev_val <- (elev_val - mean(occ_data$elev.raw)) / sd(occ_data$elev.raw)
precc_val <- (precc_val - mean(occ_data$precc.raw)) / sd(occ_data$precc.raw)
fir_val <- (fir_val - mean(occ_data$fir.raw)) / sd(occ_data$fir.raw)
age_val <- (age - mean(occ_data$age.raw, na.rm = TRUE)) / sd(occ_data$age.raw, na.rm = TRUE)
age_val <- rep(age_val, length.out = length(cc_val))

## WHR is a bit more tricky - match categories to IDs used in model (10 categories)
`%notin%` <- Negate(`%in%`)
veg_val <- data.frame(WHR = fire_whr_crosswalk$WHRTYPE[match(veg_val_ID, fire_whr_crosswalk$Value)])
veg_val$Life_form <- fire_whr_crosswalk$LIFE_FORM[match(veg_val_ID, fire_whr_crosswalk$Value)]

veg_val$WHR[veg_val$WHR %notin% c("EPN","JPN","LPN","PPN","RFR","SMC","WFR")] <- NA  # Forest types to save
veg_val$WHR[veg_val$Life_form %in% c("BARREN/OTHER","HERBACEOUS","SHRUB")] <- "NFR"  # Non-forest
veg_val$WHR[veg_val$Life_form == "HARDWOOD"] <- "HRD"                                # Hardwood
veg_val$WHR[is.na(veg_val$WHR)] <- "OTH"                                             # Other conifers
veg_val <- factor(veg_val$WHR, levels=levels(occ_data$whrtype))
veg_val <- as.numeric(veg_val)

## Scale fire size variable
fire_size <- input_size/2.471 # convert from acres to ha
size_val <- (log(fire_size) - mean(occ_data$firesize.log)) / sd(occ_data$firesize.log) # log it and standardize

## Scale the fire season variable: "ig" refers to ignition date
if(input_season == "January 1 to August 15"){
  ig_val <- 0
} else{
  ig_val <- 1
}


#### Last thing -- indicator variable for medium or high severity patches
cc30      <- values(fire_cc)[pix]
Ind <- rep(0, length(cc30))
Ind[which(cc30 > 25)] <- 1
rm(cc30)


#### -------------------- Write function for running the models!
`expit` <- function(x) {exp(x) / (1 + exp(x))}

mod_predict <- function(distrib_snag, distrib_psi, distrib_hr, distrib_firelevel, distrib_psi_b0){
  # Snag density
  snag_sim <- exp(distrib_snag[1] + 
                  distrib_snag[2] * snag_cc_val +
                  distrib_snag[3] * snag_cc_val*snag_cc_val + 
                  distrib_snag[4] * snag_precc_val + 
                  distrib_snag[5] * snag_precc_val * snag_cc_val + 
                  distrib_snag[6] * snag_size_val +  
                  distrib_snag[7] * snag_elev_val +
                  distrib_snag[8] * snag_lat_val +
                  distrib_snag[9] * snag_elev_val * snag_lat_val)
  snag_sim <- snag_sim*0.2295687*10  # Convert to m2/ha
    
  # Occupancy
  fire_b0 <- distrib_firelevel[1] +                       # Generate the fire-level intercept
    distrib_firelevel[2] * size_val + 
    distrib_firelevel[3] * ig_val 
  psi_sim <- expit(fire_b0 +                              # Fire-level intercept
                     distrib_psi_b0[veg_val] +            # WHR type intercept
                     distrib_psi[1] * age_val +  
                     distrib_psi[2] * age_val*age_val +
                     distrib_psi[3] * elev_val + 
                     distrib_psi[4] * elev_val*elev_val + 
                     distrib_psi[5] * lat_val + 
                     distrib_psi[6] * cc_val + 
                     distrib_psi[7] * pyro_val + 
                     distrib_psi[8] * patch_val*Ind + 
                     distrib_psi[9] * precc_val +  
                     distrib_psi[10] * elev_val*lat_val + 
                     distrib_psi[11] * age_val*cc_val +
                     distrib_psi[12] * fir_val +
                     distrib_psi[13] * age_val*fir_val) 

  # Home-range size
  hr_sim <- (exp(distrib_hr[1] + distrib_hr[2] * snag_sim))
  hr_sim[hr_sim < 20] <- 20  # Set minimum home range as 20 ha
  hr_sim[hr_sim > 825] <- 825  # Set maximum home range as 825 ha
  hr_sim <- 1 / hr_sim  # Density given occupancy, pairs per hectare
  hr_sim <- hr_sim * (((xres(fire_cc)) ^ 2) / (100 ^ 2))  # Density given occupancy, pairs per cell
  
  # Expected density
  dens_vals <- hr_sim * psi_sim
  dens_vals[is.na(dens_vals)] <- 0  # Force uninhabitable areas (NAs) to a density of 0
  
  # Output
  return(dens_vals)
}


#### ---------------------------------------- 2 ---------------------------------------- ####


#### ------------- Predict BBWO abundance across fire 

nsim      <- 1000  # an arbitrarily large number over which to sample from posterior. 
abund_est <- array(dim = nsim)  # output file

start.time<-Sys.time()
for(i in 1:nsim) {
  ## Pick parameter values for this sim from each posterior
  psi_i <- round(runif(1, 1, length(post_psi[[2]])))
  distrib_psi <- c(post_psi$b.age[psi_i], 
                   post_psi$b.agesq[psi_i], 
                   post_psi$b.elev[psi_i], 
                   post_psi$b.elevsq[psi_i], 
                   post_psi$b.lat[psi_i], 
                   post_psi$b.sev[psi_i],
                   post_psi$b.pyro[psi_i], 
                   post_psi$b.patch[psi_i], 
                   post_psi$b.precc[psi_i],
                   post_psi$b.elevlat[psi_i], 
                   post_psi$b.agesev[psi_i],
                   post_psi$b.fir[psi_i],
                   post_psi$b.agefir[psi_i])
  distrib_psi_b0 <- post_b.whr[psi_i, ]  # intercepts based on 10 WHR-types
  distrib_firelevel <- c(post_psi$f0[psi_i],
                         post_psi$f.size[psi_i],
                         post_psi$f.season[psi_i])
  hr_i <- round(runif(1, 1, length(post_hr[[1]])))
  distrib_hr <- c(post_hr$b0[hr_i], post_hr$b1[hr_i])
  snag_i <- round(runif(1, 1, length(post_snag[[1]])))
  distrib_snag <- c(post_snag$b0[snag_i], 
                    post_snag$b.bs[snag_i], 
                    post_snag$b.bs2[snag_i], 
                    post_snag$b.precc[snag_i], 
                    post_snag$b.bs.precc[snag_i], 
                    post_snag$b.size[snag_i], 
                    post_snag$b.elev[snag_i], 
                    post_snag$b.lat[snag_i], 
                    post_snag$b.elevlat[snag_i])  
  
  ## Use parameter values to predict density across the fire
  dens_vals <- mod_predict(distrib_snag = distrib_snag, distrib_psi = distrib_psi, 
                           distrib_hr = distrib_hr, distrib_firelevel = distrib_firelevel, 
                           distrib_psi_b0 = distrib_psi_b0)
 
  ## Sum density values to get abundance
  abund_est[i] <- sum(dens_vals)
  print(paste(round((i / nsim) * 100, 2), "%"))  # Calculation timer
}
end.time=Sys.time()
elapsed.time <- difftime(end.time, start.time, units='mins')
elapsed.time

## abund_est is a vector giving n=nsim estimates of potential abundance across the prediction area
hist(abund_est)

#### ------------------------------------------- 3 ---------------------------------------- ####


#### -------- Predict mean BBWO density across fire 

## Get posterior mean for each parameter
distrib_psi <- c(mean(post_psi$b.age), 
                 mean(post_psi$b.agesq), 
                 mean(post_psi$b.elev), 
                 mean(post_psi$b.elevsq), 
                 mean(post_psi$b.lat), 
                 mean(post_psi$b.sev),
                 mean(post_psi$b.pyro), 
                 mean(post_psi$b.patch), 
                 mean(post_psi$b.precc), 
                 mean(post_psi$b.elevlat), 
                 mean(post_psi$b.agesev),
                 mean(post_psi$b.fir),
                 mean(post_psi$b.agefir))
distrib_psi_b0 <- apply(post_b.whr, 2, mean)
distrib_firelevel <- c(mean(post_psi$f0),
                       mean(post_psi$f.size),
                       mean(post_psi$f.season))
distrib_hr <- c(mean(post_hr$b0), mean(post_hr$b1))
distrib_snag <- c(mean(post_snag$b0), 
                  mean(post_snag$b.bs), 
                  mean(post_snag$b.bs2), 
                  mean(post_snag$b.precc),
                  mean(post_snag$b.bs.precc),
                  mean(post_snag$b.size), 
                  mean(post_snag$b.elev), 
                  mean(post_snag$b.lat),
                  mean(post_snag$b.elevlat))

## Use mean parameter values to predict density across the fire
dens_vals <-  mod_predict(distrib_snag = distrib_snag, distrib_psi = distrib_psi, 
                          distrib_hr = distrib_hr, distrib_firelevel = distrib_firelevel, 
                          distrib_psi_b0 = distrib_psi_b0)


#### Post-processing of desired outputs
## Convert vector of predicted density into a raster
dens_mean_val <- rep(NA, length(in_fire))
dens_mean <- fire_outline_rast
values(dens_mean)[pix] <- 0    # we will map density values onto this empty raster

dens_mean_val[pix] <- dens_vals
values(dens_mean)  <- dens_mean_val

## dens_mean is a raster giving predicted density across the fire area (pairs per pixel)
plot(dens_mean)


#### -------------------------------------------- 4 -------------------------------------------- ####


#### ----------- Predict density uncertainty across the fire
if(input_size >= 250000){
  nsim <- 300
} else if(input_size < 250000 & input_size >= 100000){
  nsim <- 400
} else{
  nsim <- 500
}

dens_sim <- matrix(nrow = length(pix), ncol = nsim)  # output file

start.time<-Sys.time()
for(i in 1:nsim) {
  ## Pick parameter values for this sim from each posterior
  psi_i <- round(runif(1, 1, length(post_psi[[2]])))
  distrib_psi <- c(post_psi$b.age[psi_i], 
                   post_psi$b.agesq[psi_i], 
                   post_psi$b.elev[psi_i], 
                   post_psi$b.elevsq[psi_i], 
                   post_psi$b.lat[psi_i], 
                   post_psi$b.sev[psi_i],
                   post_psi$b.pyro[psi_i], 
                   post_psi$b.patch[psi_i], 
                   post_psi$b.precc[psi_i],
                   post_psi$b.elevlat[psi_i], 
                   post_psi$b.agesev[psi_i],
                   post_psi$b.fir[psi_i],
                   post_psi$b.agefir[psi_i])
  distrib_psi_b0 <- post_b.whr[psi_i, ]  # intercepts based on 10 WHR-types
  distrib_firelevel <- c(post_psi$f0[psi_i],
                         post_psi$f.size[psi_i],
                         post_psi$f.season[psi_i])
  hr_i <- round(runif(1, 1, length(post_hr[[1]])))
  distrib_hr <- c(post_hr$b0[hr_i], post_hr$b1[hr_i])
  snag_i <- round(runif(1, 1, length(post_snag[[1]])))
  distrib_snag <- c(post_snag$b0[snag_i], 
                    post_snag$b.bs[snag_i], 
                    post_snag$b.bs2[snag_i], 
                    post_snag$b.precc[snag_i], 
                    post_snag$b.bs.precc[snag_i], 
                    post_snag$b.size[snag_i], 
                    post_snag$b.elev[snag_i], 
                    post_snag$b.lat[snag_i], 
                    post_snag$b.elevlat[snag_i])  
  
 
  ## Use parameter values to predict density across the fire
  dens_vals <- mod_predict(distrib_snag = distrib_snag, distrib_psi = distrib_psi, 
                           distrib_hr = distrib_hr, distrib_firelevel = distrib_firelevel, 
                           distrib_psi_b0 = distrib_psi_b0)
  
  ## Add predicted denstiy to matrix
  dens_sim[,i] <- dens_vals
  print(paste(round((i / nsim) * 100, 2), "%"))  # Calculation timer
}
end.time=Sys.time()
elapsed.time <- difftime(end.time, start.time, units='mins')
elapsed.time


#### Post-processing of desired outputs
uncert_vals <- apply(dens_sim, 1, function(x) sd(x) / mean(x)) 
#uncert_vals <- apply(dens_sim, 1, function(x) sqrt(exp(sd(log(x))^2)-1)) # geometric CV
rm(dens_sim)  # save space

## Convert vector of density uncertainty into a raster
dens_uncert_val <- rep(NA, length(in_fire))
dens_uncert <- fire_outline_rast
values(dens_uncert)[pix] <- 0    # we will map density values onto this empty raster

dens_uncert_val[pix] <- uncert_vals
values(dens_uncert)  <- dens_uncert_val

## Reclassify coefficient of variation into 3 categories
#3 Good = < 0.4
#2 Fair = [0.4,1)
#1 Poor = > 1
reclass_df <- c(-Inf, 0.400, 3,    # Good 
                0.400, 1.00, 2,    # Fair
                1.00, Inf, 1)      # Poor
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)

dens_class <- reclassify(dens_uncert, reclass_m)

## dens_class gives density uncertainty in three categories.
## These categories are easily changed (e.g., from a cutoff at 0.400) based on desired thresholds.
plot(dens_class)



