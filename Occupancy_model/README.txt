This file contains information and explanation for the data that accompany the following probject:

Stillman, A.N., D.R. Kaschube, R.L. Wilkerson, R.B. Siegel, S.C. Sawyer, and M.W. Tingley. 2021. Incorporating pyrodiversity into 
wildlife habitat assessments for post-fire management. Joint Fire Science Program.

This .README file accompanies the archived folder "Occupancy_model". Individual .csv files are read into the script "Run_OuccpancyModel" 
and reformatted for the model.
Continuous variables were standardized prior to modeling by substracting the mean and dividing by the standard deviation. Standardized 
variables are denoted by "S" before the variable name.  
__________________________________


Dataset: detection_covars.csv
ef = binary effort variable. 0 = 2 minute survey, 1 = 3 minute survey. 
itype = survey type. 0 = passive survey, 1 = broadcast survey. 

Dataset: fire_covars.csv
S.fire.size = the size of the fire area, in hectares. The natural log was applied before standardization.
fire.season = binary value for ignition date. 0 = fire started before or on August 15th, 1 = fire started after August 15th. 

Dataset: S.fire.age.csv
Matrix with nrow = the number of survey points, and columns = 1 - 10 years post fire. Each row represents 1-10 years post-fire. 

Dataset: S.jday.csv
The standardized ordinal date of the woodpecker survey at each site. nrow = the number of survey points, and columns = 1 - 10 years post fire.
Years with no woodpecker surveys are given the average value after standardization (~0). 

Dataset: site_covars.csv
Site-level covariate data at each survey point.
S.elev = Elevation at the survey location (m).
S.lat = Latitude at the survey location (decimal degrees).
bs30 = Burn severity of the pixel that overlaps the survey location, measured as the % change in canopy cover
	from before to immediately after fire.
S.bs100 = Mean burn severity within 100m of the survey point. 
S.d.patch = Distance from the survey point to the nearest pixel that burned at low severity (<25% change in canopy cover).
S.pyro500 = Inverse Simpsons diversity of burn severity pixels (divided into 11 classes) within 500 m of the survey point. 
S.precc100 = Pre-fire canopy cover within 100 m of the survey point. 
S.fir100 = Combined basal area of red and white fir within 100 m of the survey point. 
whr = California Wildlife Habitat Relationships classification at the survey point.
site.code = Unique ID assigned to each survey location.
fireID = Unique ID assigned to each fire where surveys took place. 

Datasets: X.array.yr*
These datasets provide the response variable for the model and represent the detection history for each point. Each row is a different survey point and 
each column is a different period of a single survey (see Methods). There are 10 .csv files, each representing a different year 1 - 10 years post-fire. 
NA = No survey conducted
0 = Black-backed Woodpecker not detected
1 = Black-backed Woodpecker detected
