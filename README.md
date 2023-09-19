# RoyleTurner2022
Code and data from this paper:
Royle, J.A. and Turner, H., 2022. Density estimation in terrestrial chelonian populations using spatial capture–recapture and search–encounter surveys. Journal of Herpetology, 56(3), pp.341-348.

There are three scripts that you need (only 2 uploaded so far, working on the 3rd!):
Script 1 processes all of the .gpx files and buffers each track. A number of objects are produced and saved. These are used when Script 2 is run.  It would be easy to automate and improve this script which I aim to do.

Script 2 summarizes all of the search tracks to compute a spatial effort covariate (number of searches through each grid cell basically). This can take awhile to run. 

Script 3: (not done yet)
Organizes the capture data and search effort data into the encounter data file (EDF) and trap deployment file (TDF) that SCR models require. Then some SCR models are fitted using the oSCR package. 


There are a number of .RData files included as well:
   1. utils.RData -- 3 utility functions used to make quick and dirty plots
   2. spatial_data.RData -- produced from Script 1 -- the buffered GPS tracks from the 2020 survey effort, including a few other objects created in Script 1.
   3. effortmatrix.RData -- this is produced by Script 2, it is the gridded effort covariate used in the paper. (basically it combines all buffered tracks, overlays the grid, addes them up).

      
