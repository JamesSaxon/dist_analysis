# C4 Political Analysis: Legislative Avenues for Gerrymandering Reform

## Overview of the Analysis

There are four components to the C4 Analysis.  These are divided between two repositories.
1. *Geographies*: these are shapefiles drawn from the US Census, preprocessed using postgres/postgis.  The scripts used to build this database can be found in the [c4 repository in the `db/` directory](https://github.com/JamesSaxon/C4/tree/master/db).  A [README](https://github.com/JamesSaxon/C4/blob/master/db/README.md) in that folder describes those scripts.  To avoid dependence on this private database, the data are also cached and saved in [`c4/shapes/`](https://github.com/JamesSaxon/C4/tree/master/shapes) and [`c4/demographics/`](https://github.com/JamesSaxon/C4/tree/master/demographic).
   * *Replication suggestion*: check the import scripts and verify that the cached data in `shapes/` are sensible.
2. *Election data*: these data are drawn from the Census Voter Tabulation Districts from 2010, state and county election websites, and the work by [Ansolabehere, Rodden et al.](https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/21919).  The scripts for munging these input data are in this repository, in the [`voting/`](https://github.com/JamesSaxon/district_analysis/tree/master/voting) directory.  The script for each individual state describes the sources and the strategy for that state (they are pretty consistent).  This generates two outputs: [`mapped/`](https://github.com/JamesSaxon/district_analysis/tree/master/voting/mapped) precinct-level geojsons for each race, and state [`votes/`](https://github.com/JamesSaxon/district_analysis/tree/master/voting/votes) csvs.  The geojson files display natively on github -- they are zoomable and colored, and lend credence that the basic process has worked (cities are blue, rural areas are red).  The vote csvs are in the `c4` software (`demographic/`) to calculate winners (where possible), but because these could be (and were) updated after running simulation, the actual results use those files to recalculate winners and competitiveness by simply merging tracts with `/voting/votes/??.csv`.
   * *Replication suggestions*: See any of the geojson maps, comparing e.g. [IL](https://github.com/JamesSaxon/district_analysis/blob/master/voting/mapped/il_2016.geojson) or [LA](https://github.com/JamesSaxon/district_analysis/blob/master/voting/mapped/la_2016.geojson) in 2016, or MN across the years ([2008](https://github.com/JamesSaxon/district_analysis/blob/master/voting/mapped/mn_2008.geojson), [2012](https://github.com/JamesSaxon/district_analysis/blob/master/voting/mapped/mn_2012.geojson), [2016](https://github.com/JamesSaxon/district_analysis/blob/master/voting/mapped/mn_2016.geojson)).
3. *Clustering software (c4)*: the "main contribution" of this project is the [c4](https://github.com/JamesSaxon/C4/) software, which occupies its own repository, and has its own README.  The software is written in c++, bound to python through cython.  Preprocessing of the data is performed through python, using geopandas and pysal, which in turn depend on gdal.  The software _should_ be relatively easy to install -- it is meant to be accessible -- but may be a bit tricky.  The software ran within a Docker container ([AWS Dockerfile](https://github.com/JamesSaxon/C4/blob/master/Dockerfile), with jobs launched to [run on aws ec2/batch](https://github.com/JamesSaxon/C4/blob/master/aws_launch.sh) and write output to the S3 storage).  I then downloaded these data for analysis.
   * *Replication suggestion 1*: First, build the package or pull the non-Amazon docker container to run the code and generate districts.
     Instructions for this are in the other repo, but it basically comes to one line of code:
     ```
     docker run -v $(pwd)/res/:/C4/res/ \
                -e STATE=pa -e SEED=2 -e METHOD=POWER \
                --rm -it  jamessaxon/c4:replication
     ```
     Only the 11 states with voting data (FL, IL, LA, MD, MN, NC, PA, TX, TN, VA, WI) will generate winners and losers. 
     Other states should all work, though California may get slow.
   * *Replication suggestion 2*: Then head to this [interactive map](https://saxon.harris.uchicago.edu/redistricting_map/), where many maps (but not all) can be viewed.  Please consider, however, that this is an outreach project and not an "official" part of the analysis.  Where there are discrepancies, the official analysis has priority.  In particular, the map has a fault in the Virginia voting data (excluded from the paper, though fixed in this repo).
4. *Analysis*: This is the work of most of the scripts in _this_ directory/repo, which picks up after districts have been simulated on AWS.
   * *Replication suggestion*: run the code to reproduce figures, as described below.
   
## Processing Code in this Repository

### Extracting state json files from bulk sim
### Merging simulation and voting data
### Evaluating the spatial properties of historical districts

## Generating Tables and Figure

### State-Level Seat Share Tables
### Competitiveness Table
### Minority Seat Share
### Appendix: County Splits
### Appendix: PCA of Historic Seats
