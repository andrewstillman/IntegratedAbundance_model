# IntegratedAbundance_model
Data and code to run an integrated space use and occupancy model for Black-backed Woodpeckers

This repository contains information and explanation for the data and code that accompany the following probject:
Stillman, A.N., D.R. Kaschube, R.L. Wilkerson, R.B. Siegel, S.C. Sawyer, and M.W. Tingley. 2021. Incorporating pyrodiversity into 
wildlife habitat assessments for post-fire management. Joint Fire Science Program.

Code was originally written for the RShiny app hosted by the Institute for Bird Populations (IBP).
We recommend that users run the model using the app (available at the IBP website). The code posted here provides 
an example using the Reading Fire (Lassen National Forest, 2012) and pre-made spatial layers covering the fire area. 


/IntegratedAbundance_model: data and code to run the integrated space use and occupancy model to predict Black-backed Woodpecker 
abundance across an example fire (fire data contained in /ReadingFire_data). Predictions come from three component models: 
(1) A model of the relationship between snag density and home range size. Available at https://github.com/mtingley/BBWO_abundance/tree/master/Model_HomeRange
(2) A model to predict snag denstiy given remote sensing data.
(3) A temporal auto-logistic occupancy model based on data from a long-term population survey.

/Occupancy_model: data and code to run the occupancy model component. 

/ReadingFire_data: example input data for the integrated abundance model. 

/Snag_model: data and code to run the snag density model component. 
