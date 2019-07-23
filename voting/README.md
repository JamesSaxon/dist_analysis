# Voting Data

The voting data for this project were drawn from a number of sources:
  the excellent and ambitious project by [Ansolabehere, Rodden et al.](https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/21919) (mostly for earlier data),
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
  precincts were first merged to Census tracts that contained their centroids,
  and a nearest neighbor match was then used to match any unmatched precincts
  (see [`../dist_tools.py`](https://github.com/JamesSaxon/district_analysis/blob/master/dist_tools.py#L100)).
  
In Louisiana, North Carolina, and Tennessee, provisional and absentee ballots
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
I end each notebook with a call of
  [map_sanity_check()](https://github.com/JamesSaxon/district_analysis/blob/master/dist_tools.py#L117]),
  which reads in the completed data along with the Census tract geometry.
While developing these scripts, I had additional checks to ensure that no votes were lost in the process.

### Included States and Sources

There are (brief!) notes on the strategy for each state, at the outset of each notebook.
This mirrors exactly the data description from the appendix of the paper.

* Florida (2008): [Ansolabehere and Rodden, 2011, "Florida Data Files"](https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/16797)
* Ilinois (2008, 2016)
  * Data for 2008 is [Ansolabehere, Rodden 2011](https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/15845)
  * Voting returns for 2016 come from the state.  I downloaded these in May 2017.  The link has since changed, and the precinct-level presidential returns now seem to be [here](https://www.elections.il.gov/Downloads/ElectionOperations/ElectionResults/ByOffice/51/51-120-PRESIDENT%20AND%20VICE%20PRESIDENT-2016GE.csv), but data may have changed.
  * Geographies for 2016 are an amalgamation of the the 2010 VTDs with updated data for "Chicagoland" counties 
    [Chicago](https://data.cityofchicago.org/Facilities-Geographic-Boundaries/Precincts-current-/uvpq-qeeq) and 
    the [balance of Cook County](https://datacatalog.cookcountyil.gov/GIS-Maps/Historical-ccgisdata-Election-Precinct-Data-2015-t/mtie-43p4), 
    [DuPage](http://gisdata-dupage.opendata.arcgis.com/datasets/election-precincts), and 
    [Lake](http://data-lakecountyil.opendata.arcgis.com/datasets/voting-precincts-1).
    Each of these links have changed, and several of them (Chicago, DuPage, and Lake)
    seem now to point to current precincts rather than those used in 2016.
    These counties account for the vast majority of precinct changes, and about 7/12 of Illinois's population.
    The rest of the state is matched by precinct and county name. (Notes in file.)
* Louisiana (2012, 2016)
  * Presidential Votes by Precinct:
    * https://voterportal.sos.la.gov/static/2012-11-06/resultsRegion/46257
    * https://voterportal.sos.la.gov/static/2016-11-08/resultsRegion/53898
  * Watching the data on the network shows that json is retrievable through calls like this, for each county (1-64):
    * https://voterportal.sos.la.gov/ElectionResults/ElectionResults/Data?blob=20161108/VotesRaceByPrecinct/Votes_53898_01.htm
  * Precinct-Level Shapefiles:
    * http://house.louisiana.gov/H_Redistricting2011/default_LouisianaPrecinctShapefiles.htm
* Maryland (2008, 2016)
  * 2008 Vote Totals are from [Ansolabehere, Palmer, and Lee, 2014](https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/21919); shapes are from Census VTDs.  There is some silliness as to "Maryland Counties", which are not FIPS codes.
  * 2016 Data:
    * [Polling place addresses](https://elections.maryland.gov/elections/2016/2016_Precincts_and_polling_places_GENERAL.xls) are from the state.  These are geocoded (with some manual assists) through Google.
    * [General election precinct returns](https://elections.maryland.gov/elections/2016/election_data/All_By_Precinct_2016_General.csv)
* Minnesota (2008-2016): This is simply the best data in the country.  Votes and precincts are merged cleanly, and the links are reliable.
  * Clearinghouse: https://www.gis.leg.mn/html/download.html
  * ftp://ftp.commissions.leg.state.mn.us/pub/gis/shape/elec2016.zip
  * ftp://ftp.commissions.leg.state.mn.us/pub/gis/shape/elec2012.zip
  * ftp://ftp.commissions.leg.state.mn.us/pub/gis/shape/elec08.zip
* North Carolina (2012, 2016)
  * Shapefiles for [2012](https://s3.amazonaws.com/dl.ncsbe.gov/PrecinctMaps/SBE_PRECINCTS_09012012.zip) [2016](https://s3.amazonaws.com/dl.ncsbe.gov/PrecinctMaps/SBE_PRECINCTS_20161004.zip): http://dl.ncsbe.gov/index.html?prefix=PrecinctMaps/ 
  * Precinct voting returns: https://er.ncsbe.gov/downloads.html?election_dt=11/08/2016
* Pennsylvania (2000-2012)
  * Votes are from 2000 to 2012 are from [Ansolabehere, Palmer, and Lee, 2014](https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/21919).
  * Shapefiles are the Census VTDs, which require some name cleaning in 2000-2008, and much more in 2012.  The latter was accomplished with some explicit fixes, followed by human-verified fuzzy matches (see dist_tools for Jaro-Winkler one-liner).
* Tennessee (2016)
  * Precinct returns: http://sos-tn-gov-files.s3.amazonaws.com/StateGeneralbyPrecinctNov2016.xlsx
  * 3D KMZ files are split by counties.  The location has changed since I downloaded these.  The landing page is now: https://apps.cot.tn.gov/DPAMaps/Redistrict/Counties.  A script in the repo shows the download.
* Texas (2000-2016)
  * 2008-2008 uses election returns by [Ansolabehere, et al. (2015)](https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/21919)
  * For 2012 and 2016, the [shapefiles](ftp://ftpgis1.tlc.state.tx.us/2011_Redistricting_Data/VTDs/Geography/) and 
  [election returns](ftp://ftpgis1.tlc.state.tx.us/elections) are from the Texas Secretary of State.
* Virginia (2016)
  * Initial shapefile from the [VA Public Access Project](https://github.com/vapublicaccessproject/va-precinct-maps-2016).  These were faulty/mismatched for Roanoke, and new data was obtained from the county GIS office.
  * Election _returns_ are from the state's excellent [historical election viewer](https://historical.elections.virginia.gov/elections/view/80871/), which hides (very slightly) an [API](http://historical.elections.virginia.gov/elections/download/80871/precincts_include:1/).
  * Note that these data were not used for the project, because I messed up the merge, before submitting the first time.  (This has since been fixed, but Virginia is still not in the paper.)
* Wisconsin (2004-2016): Nearly-Minnesota quality data (thanks, _Gill v. Whitford!_), so not much to do. The addresses have changed since I originally downloaded these, but the formats seem largely the same.
  * http://legis.wisconsin.gov/ltsb/gis/data/
  * https://data-ltsb.opendata.arcgis.com/search?tags=election%20data
  * https://data-ltsb.opendata.arcgis.com/datasets/2002-2010-wi-election-data-with-2017-wards
  * https://data-ltsb.opendata.arcgis.com/datasets/2012-2020-wi-election-data-with-2017-wards
