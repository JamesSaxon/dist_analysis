# C4 Political Analysis: Legislative Avenues for Gerrymandering Reform

## Overview of the Analysis

There are four components to the C4 Analysis.  These are divided between two repositories.
1. *Geographies*: these are shapefiles drawn from the US Census, preprocessed using postgres/postgis.  The scripts used to build this database can be found in the [c4 repository in the `db/` directory](https://github.com/JamesSaxon/C4/tree/master/db).  A [README](https://github.com/JamesSaxon/C4/blob/master/db/README.md) in that folder describes those scripts.  To avoid dependence on this private database, the data are also cached and saved in [`c4/shapes/`](https://github.com/JamesSaxon/C4/tree/master/shapes) and [`c4/demographics/`](https://github.com/JamesSaxon/C4/tree/master/demographic).
   * *Replication suggestion*: check the import scripts and verify that the cached data in `shapes/` are sensible.
2. *Election data*: these data are drawn from the Census Voter Tabulation Districts from 2010, state and county election websites, and the work by [Ansolabehere, Rodden et al.](https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/21919).  The scripts for munging these input data are in this repository, in the [`voting/`](https://github.com/JamesSaxon/district_analysis/tree/master/voting) directory.  The script for each individual state describes the sources and the strategy for that state (they are pretty consistent).
   
   This generates two outputs: [`mapped/`](https://github.com/JamesSaxon/district_analysis/tree/master/voting/mapped) precinct-level GeoJSONs for each race, and state [`votes/`](https://github.com/JamesSaxon/district_analysis/tree/master/voting/votes) CSVs.  The GeoJSON files display natively on github -- they are zoomable and colored, and lend credence that the basic process has worked (cities are blue, rural areas are red).  The vote CSVs are copied to the `c4` software (`demographic/`) to calculate winners (where possible), but because these could be (and were) updated after running simulation, the actual results of the paper use the original `voting/votes/??.csv` files to recalculate winners and competitiveness by simply merging them with the tract cluster assignments from the clustering simulation.
   * *Replication suggestions*: See any of the geojson maps, comparing e.g. [IL](https://github.com/JamesSaxon/district_analysis/blob/master/voting/mapped/il_2016.geojson) or [LA](https://github.com/JamesSaxon/district_analysis/blob/master/voting/mapped/la_2016.geojson) in 2016, or MN across the years ([2008](https://github.com/JamesSaxon/district_analysis/blob/master/voting/mapped/mn_2008.geojson), [2012](https://github.com/JamesSaxon/district_analysis/blob/master/voting/mapped/mn_2012.geojson), [2016](https://github.com/JamesSaxon/district_analysis/blob/master/voting/mapped/mn_2016.geojson)).
3. *Clustering software (c4)*: the "main contribution" of this project is the [c4](https://github.com/JamesSaxon/C4/) software, which occupies its own repository, and has its own README.  The software is written in c++, bound to python through cython.  Preprocessing of the data is performed through python, using geopandas and pysal, which in turn depend on gdal.  The software _should_ be relatively easy to install -- it is meant to be accessible -- but may be a bit tricky.  The software ran within a Docker container ([AWS Dockerfile](https://github.com/JamesSaxon/C4/blob/master/DockerfileAWS)), with jobs [launched to AWS ec2/batch](https://github.com/JamesSaxon/C4/blob/master/aws_launch.sh), and output staged to S3 storage.  I then downloaded these data for analysis.
   * *Replication suggestion 1*: First, build the package or pull the non-Amazon docker container to run the code and generate districts.
     Only the 11 states with voting data (FL, IL, LA, MD, MN, NC, PA, TX, TN, VA, WI) will generate winners and losers. 
     Other states should all work, though California may get slow.
     Instructions for this are in the [c4 repo](https://github.com/JamesSaxon/C4#running-c4-as-a-docker-container-simple),
     but running from the image on [DockerHub](https://cloud.docker.com/repository/docker/jamessaxon/c4) is a one-liner:
     ```
     docker run -e STATE=pa -e SEED=2 -e METHOD=POWER -e SHADING=all \
                -v $(pwd)/res/:/C4/res/ \
                --rm -it  jamessaxon/c4:replication
     ```
     The results will be written to `res/` on your machine.
   * *Replication suggestion 2*: Then head to this [interactive map](https://saxon.harris.uchicago.edu/redistricting_map/), where many maps (but not all) can be viewed.  Please consider, however, that this is an outreach project and not an "official" part of the analysis.  Where there are discrepancies, the official analysis has priority.  In particular, the map has a fault in the Virginia voting data (excluded from the paper, though fixed in this repo).
4. *Analysis*: This is the work of most of the scripts in _this_ directory/repo, which picks up after districts have been simulated on AWS.
   * *Replication suggestion*: run the code to reproduce figures, as described below.
   
## Data and Processing Code in this Repository

### Extracting state json files from bulk sim
The C4 program [outputs](https://github.com/JamesSaxon/c4#outputs) 
  json, csv, and geojson data for generated districts.
The csv files are simple tract to seat assignments,
  but the json and geojson contain fairly broad descriptions of district demographics,
    partisanship, spatial scores, population deviations, and so forth.

The full simulation came to around 30k maps, for each of the 11 core states.
It was computationally inefficient to load each of these individually.
I therefore used the excellent json manipulation library
  [`jq`](https://stedolan.github.io/jq/manual/)
  to combine these into single files.
The jq scripts are as follows: 
* The main jq manipulation is `extract.jq` which is run by `extract.sh`.
  This is the code that pulls out all voting returns,
  and creates a single output file for the main tables of the paper.
  These files are not included; instead, the files with final votes are (next section).
* For Appendix G (post-selection of minority majority districts),
  I select plans with at least two "minority" seats 
  (defined as both over half Democratic and more than 30\% Black VAP),
  and record the two Democratic seats with the highest share Black.
  These are scripts `nc_race_extract.jq` and `nc_race_shares.jq`.

### Merging district simulation and voting data

I updated the vote aggregation process several times, 
  by improving the spatial correspondence between precincts and tracts.
I also corrected a bug that affected two Texas elections, at an earlier stage of this review.

These updated vote tallies of course affect partisan outcomes.
The [`update_election_json.ipynb`](update_election_json.ipynb) notebook performs this update.
The code is closely commented and should be self-explanatory.
It writes files for analysis to `data/{date}/{st}_redux.json`.
The committed version of the data is the one created
  while preparing the replication materials, [`data/190717/`](data/190717/).

### Evaluating the spatial properties of historical districts

For appendix J (correlations and principal component analysis of enacted maps),
  the compactness scores of enacted districts are needed.
These scores are evaluated using geopandas instead of C4.
The functions are defined in `dist_tools.py`.
They are called in `run_historic_tracts.py`.
This evaluates spatial scores at the same granularity as C4,
  although some approximations differ.
For example, C4 replaces cell geometries with centroids in many cases (like convex hull) 
    where `dist_tools.py` uses the actual geometry.
The data are saved to `data/decennial_census.*` for inspection.

The similar script, `run_historic_blocks.py`
  was used to calculate spatial and demographics scores at finer 
  granularity, for the [webmap](https://saxon.harris.uchicago.edu/redistricting_map/).
In that case, the exact populations were necessary,
  since I also aggregated the demographic characteristics and calculated the population deviance.
But the same functions are called.

### Power Diagram Districts

* state maps for race!!

## Generating Tables and Figure

The scripts in this section are the ones
  that the replicator is expected to run.

### State-Level Seat Share Tables
* Code: [`main_tables.ipynb`](main_tables.ipynb)

### Competitiveness Table
* Code: [`main_tables.ipynb`](main_tables.ipynb)

### Minority Seat Share
* Code: [`minority_seats.ipynb`](minority_seats.ipynb)

### Pennsylania Plans
### Appendix: Non-Compact Districts
### Appendix: North Carolina Race 
### Appendix: County Splits

* polsby_w

### Appendix: PCA of Historic Seats: `app_historic_interp_corr.ipynb`

