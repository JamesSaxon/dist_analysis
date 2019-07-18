# Voting Data

The voting data for this project were drawn from a number of sources:
  the excellent and ambitious project by Jonathan Rodden at Stanford (mostly for earlier data),
  the Census's 2010 Voter Tabulation Districts (VTDs/precincts), 
    and data from individiual states and even counties.
States were chosen simply based on the availability of data, 
  prioritizing states where multiple years were possible.
Shapefiles and voting returns often use slightly different names,
  so these were stitched together with manual edits and human-verified fuzzy matches.
Because I used Census tracts for the clustering,
  the process ended with a merge from precincts to clusters:
  precinct were first merged to Census tracts that contained their centroids,
  and a nearest neighbor match was used to match any unmatched precincts.
  
(It is worth pointing out that Minnesota Wisconsin have phenomenal data/shapefiles.
  The data for Texas and North Carolina require a bit more work, but are very good.)
  
The code for each state is in this directory.
The `votes/` directory holds Democratic and Republican votes by Census tract.
The `mapped/` directory contains geojson files of the precinct returns
  for each state/year, which can be viewed natively on github.
(In Maryland 2016, this works a bit less well, since I had polling places 
  instead of precinct boundaries, and the points do not display as nicely.)
  
### Florida (2008)

### Ilinois (2008, 2016)

### Louisiana (2012, 2016)

### Maryland (2008, 2016)

### Minnesota (2008-2016)

### North Carolina (2012, 2016)

### Pennsylvania (2000-2012)

### Tennessee (2016)

### Texas (2000-2016)

(Check that the online version is just Presidential, esp. for 2008 and 2012.)

### Virginia

Due to a bug (since fixed) at the time of submission, Virginia was not included in this paper.

### Wisconsin (2004-2015)
