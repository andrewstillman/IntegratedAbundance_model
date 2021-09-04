This file contains information and explanation for the data that accompany the following probject:

Stillman, A.N., D.R. Kaschube, R.L. Wilkerson, R.B. Siegel, S.C. Sawyer, and M.W. Tingley. 2021. Incorporating pyrodiversity into 
wildlife habitat assessments for post-fire management. Joint Fire Science Program.

This .README file accompanies the archived folder "IntegratedAbundance_model". Individual .csv files represent (1) data from the component models,
necessary to standardize input data for the fire of interest, and (2) posterior estimates from the three component models. The files are read 
into the script "Run_IntegratedAbundance" and reformatted to run the integrated space use and abundance model using an example fire contained
in the folder "ReadingFire_data". 

Data and code to run the snag model and occupancy model are provided in this repository. Data and code to run the home range model are
available from https://github.com/mtingley/BBWO_abundance in the folder "Model_HomeRange". 

__________________________________


Dataset: fire.age_raw.csv
Matrix with nrow = the number of survey points, and columns = 1 - 10 years post fire. Each row represents 1-10 years post-fire, 
with NA values in the years that the fire was not surveyed. 


Dataset: fire.covars.csv
Covariate values for the fire-level intercept in the occupancy model. 
firesize.log = natural log of the fire size, in hectares.
fireseason.raw = binary value for ignition date. 0 = fire started before or on August 15th, 1 = fire started after August 15th. 


Dataset: post.hr.csv
Posterior estimates for parameters from the home range model (a component of the integrated abundance model). Data and code to run the 
home range model are available from https://github.com/mtingley/BBWO_abundance in the folder "Model_HomeRange". 
b0 = posterior draws for the intercept parameter
b1 = posterior draws for the snag density parameter
deviance = estimated model deviance for each posterior draw


Dataset: post.psi.csv
Posterior estimates for parameters from the occupancy model (a component of the intergrated abundance model). 
b.whr.1 = posterior draws for the random intercept effect of CWHR habitat type. Type = "EPN".
b.whr.2 = posterior draws for the random intercept effect of CWHR habitat type. Type = "HRD".
b.whr.3 = posterior draws for the random intercept effect of CWHR habitat type. Type = "JPN".
b.whr.4 = posterior draws for the random intercept effect of CWHR habitat type. Type = "LPN".
b.whr.5 = posterior draws for the random intercept effect of CWHR habitat type. Type = "NFR" (non-forest).
b.whr.6 = posterior draws for the random intercept effect of CWHR habitat type. Type = "OTH" (other).
b.whr.7 = posterior draws for the random intercept effect of CWHR habitat type. Type = "PPN".
b.whr.8 = posterior draws for the random intercept effect of CWHR habitat type. Type = "RFR".
b.whr.9 = posterior draws for the random intercept effect of CWHR habitat type. Type = "SMC".
b.whr.10 = posterior draws for the random intercept effect of CWHR habitat type. Type = "WFR".
b.age = posterior draws for the fire age parameter.
b.agesev = posterior draws for the interaction parameter "fire age * burn severity".
b.agesq = posterior draws for the quadratic effect of fire age.
b.elev = posterior draws for the elevation parameter.
b.elevsq = posterior draws for the quadratic effect of elevation.
b.elevlat = posterior draws for the interaction parameter "elevation * latitude".
b.lat = posterior draws for the latitude parameter. 
b.patch = posterior draws for the distance to patch edge parameter.
b.precc = posterior draws for the pre-fire canopy cover parameter.
b.pyro = posterior draws for the simponson's diversity of burn severity parameter.
b.sev = posterior draws for the burn severity parameter.
b.fir = posterior draws for the fir basal area parameter.
b.agefir = posterior draws for the interaction parameter "fir basal area * fire age".
f.season = posterior draws for the fire-level ignition date parameter.
f.size = posterior draws for the fire-level fire size parameter. 
f0 = posterior draws for the fire-level intercept parameter.
b0 = posterior draws for the intercept parameter, taken as the mean of all 10 b.whr* estimates. 


Dataset: post.snag.csv
Posterior estimates for parameters from the snag density model (a component of the intergrated abundance model). 
b0 = posterior draws for the intercept parameter.
b.bs = posterior draws for the burn severity parameter.
b.bs2 = posterior draws for the quadratic effect of burn severity.
b.precc = posterior draws for the pre-fire canopy cover parameter.
b.elev = posterior draws for the elevation parameter.
b.lat = posterior draws for the latitude parameter. 
b.size = posterior draws for the CWHR tree size class parameter.
b.elevlat = posterior draws for the interaction parameter "elevation * latitude".
b.bs.precc = posterior draws for the interaction parameter "burn severity * pre-fire canopy cover".


Dataset: site.covars_raw.csv
Covariate values for the occupancy model. 
cc.raw = mean burn severity within 100m of the survey point, measured as the % change in canopy cover
	from before to immediately after fire.
pyro.raw = inverse Simpsons diversity of burn severity pixels (divided into 11 classes) within 500 m of the survey point.
patch.raw = distance from the survey point to the nearest pixel that burned at low severity (<25% change in canopy cover).
precc.raw = pre-fire canopy cover within 100 m of the survey point. 
elev.raw = elevation at the survey location (m).
lat.raw = latitude at the survey location (decimal degrees).
whrtype = california Wildlife Habitat Relationships (CWHR) classification at the survey point.
fir.raw = combined basal area of red and white fir within 100 m of the survey point. 


Dataset: snag.data.csv
This file contains mean and sd values necessary to standardize input data. 
mu.cc = mean of burn severity variable. 
sd.cc = sd of burn severity variable.
mu.precc = mean of pre-fire canopy cover variable. 
sd.precc = sd of pre-fire canopy cover variable.
mu.elev = mean of elevation variable. 
sd.elev = sd of elevation variable. 
mu.lat = mean of latitude variable. 
sd.lat = sd of latitude variable. 










