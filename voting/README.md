# Voting Data

The voting data for this project were drawn from a number of sources:
  the excellent and ambitious project by Ansolabehere, Rodden et al. (mostly for earlier data),
  [Voter Tabulation Districts (VTDs) from the Census](https://www2.census.gov/geo/tiger/TIGER2010/VTD/2010/)
    and data from individiual states and even counties.
States were chosen simply based on the availability of data, 
  prioritizing states where multiple years were possible.
All data were downloaded in mid-2017 and some links and data may, in the interim, have changed.
I have tried to update links where possible.
The original, raw data can be made available upon request.

(As an aside: Minnesota and Wisconsin have phenomenally clean data/shapefiles,
  for which I think the state GIS teams should really be commended.
 Texas is the runner-up, and North Carolina and Louisiana trail behind.)
  
Shapefiles and voting returns often use slightly different names,
  so these were stitched together with manual edits and human-verified fuzzy matches.
Because I used Census tracts for the clustering,
  the process ended with a merge from precincts to clusters:
  precinct were first merged to Census tracts that contained their centroids,
  and a nearest neighbor match was then used to match any unmatched precincts.
In the state-leve descriptions below, I will call this the _standard match_.
  
In North Carolina and Louisiana, provisional and absentee ballots
  are reported at the county level.
I reallocate Democractic and Republican votes back 
  to precincts according to each precinct's share of _that party's_
  county-wide polling place total.
In Pennsylvania, some reallocation was previously done by Rodden et co.
These all result in fractional vote tallies.

The code for cleaning each state is in this directory.
The `votes/` directory holds Democratic and Republican votes by Census tract.
The `mapped/` directory contains geojson files of the precinct returns
  for each state/year, which can be viewed natively on github.
(In Maryland 2016, this works a bit less well, since I had polling places 
  instead of precinct boundaries, and the points do not display as nicely.)

There are (brief!) notes on the strategy for each state, at the outset of each notebook.
This mirrors exactly the data description from the appendix of the paper.
* Florida (2008)
* Ilinois (2008, 2016)
* Louisiana (2012, 2016)
* Maryland (2008, 2016)
* Minnesota (2008-2016)
* North Carolina (2012, 2016)
* Pennsylvania (2000-2012)
* Tennessee (2016)
* Texas (2000-2016)
* Virginia (2016)
* Wisconsin (2004-2016)
