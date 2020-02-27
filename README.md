This is a python module for calculating various statistics from 
the radiosondes launched during EUREC4A. The soundings of the 
following platforms can be handled as long as the data is provided
in NetCDF files: 

Barbados Cloud Observatory (BCO) 
Ronald H. Brown (RHB)
Meteor (M161)
Atalante (ATL)
Maria S Merian (MSM)

The naming convention for the NetCDF files is of the form:
'{platform_abbreviation}_SoundingAscentProfile_*_{YYYYMMDD}_*.nc' 
'{platform_abbreviation}_SoundingDescentProfile_*_{YYYYMMDD}_*.nc'

The file utils.py is a module that contains functions to calculate 
various atmospheric properties based on the provided radiosonde data.
This module is used in the 3 provided jupyter notebooks to plot and to calculate
different kinds of statistics of radiosonde data, which shall here 
be briefly described:


soundings_one_platform_stats.ipynb:
-----------------------------------
* For a given platform and time period, calculate various atmospheric 
  properties for each individual sounding. 
* Calculate the mean, standard deviation, minimum, maximum and 
  quantiles over the individual atmospheric properties calculated 
  for each sounding. 
* Plot a set of vertical profiles for all the soundings of the 
  most prominent atmospheric quantities.
  
soundings_compare_platforms.ipynb:
-----------------------------------
* For a given day and time period, plot a map with the positions of each platform and the trajectories of the soundings.
* Plot vertical profiles of different quantities for all platforms.

soundings_temporal_evolution.ipynb:
-----------------------------------
* For a given platform and time period plot the temporal evolution of different quantities (relative humidity, wind speed, ...)  with height.
 
