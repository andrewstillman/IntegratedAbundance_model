This file contains information and explanation for the data that accompany the following probject:

Stillman, A.N., D.R. Kaschube, R.L. Wilkerson, R.B. Siegel, S.C. Sawyer, and M.W. Tingley. 2021. Incorporating pyrodiversity into 
wildlife habitat assessments for post-fire management. Joint Fire Science Program.

This .README file accompanies the archived folder "Snag_model". 
Continuous variables were standardized prior to modeling by substracting the mean and dividing by the standard deviation. Standardized 
variables are denoted by "S" before the variable name.  
__________________________________


Dataset: Snag_model_data.csv
Station_visit_ID = Unique ID for each survey, including the fire code, point name, and year. 
Snag_BA = Basal area of dead trees surrounding a survey point measured with a Crusier's Crutch at Basal Area Frequency = 10. 
WHRsize = California Wildlife Habitat Relationships classification of tree size at the survey point.
S.BS30 = Burn severity of the pixel that overlaps the survey location, measured as the % change in canopy cover
	from before to immediately after fire.
S.Precc30 = Pre-fire canopy cover at the pixel that overlaps the survey point. 
S.Elev = Elevation at the survey location (m).
S.Lat = Latitude at the survey location (decimal degrees).
log_Snag_BA = Natural log of snag basal area at the survey point.