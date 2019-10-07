# C4 Political Analysis: Legislative Avenues for Gerrymandering Reform

## Overview of the Analysis

There are four components to the C4 Analysis.  These are divided between two repositories.
1. *Geographies*: these are shapefiles drawn from the US Census, preprocessed using postgres (9.5.17, with postgis 2.2.1, on Ubuntu 16.04.6 LTS).
   The scripts used to build this database can be found in the [c4 repository in the `db/` directory](https://github.com/JamesSaxon/C4/tree/master/db).  A [README](https://github.com/JamesSaxon/C4/blob/master/db/README.md) in that folder describes those scripts.  To avoid dependence on this private database, the data are also cached and saved in [`c4/shapes/`](https://github.com/JamesSaxon/C4/tree/master/shapes) and [`c4/demographics/`](https://github.com/JamesSaxon/C4/tree/master/demographic).
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
   * *Replication suggestion*: run the code described below to reproduce the figures and, above all, the [main tables](https://github.com/JamesSaxon/district_analysis#state-level-seat-share-tables).

## Required Software

The analysis scripts use the following standard python libraries,
which can be installed most-easily through Anaconda, and I recommend this course strongly.
I note the versions on my own machine.
* **python (3.5.6)**
* jupyter (4.2.3)       - jupyter notebooks -- self-commenting python code.
* matplotlib (2.2.2)    - general plotting
* seaborn (0.7.1)       - plot styling...
* geopandas (0.4.1)     - mapping (require gdal)
* pandas (0.24.2)       - general analysis
* numpy (1.15.2)        - numerical simulations
* statsmodels (0.9.0)   - probit model for race
* scipy (1.1.0)         - general statistics
* scikit-learn (0.18.1) - PCA of historic districts, in appendix.
* xlrd                  - loading excel, for minority seats
* descartes             - for plotting polygons from geopandas
* pyproj (1.9.5.1)      - map projecting (c4 only)
* pysal (1.14.4)        - to generate contiguity matrix (c4 only)
* cython (0.28.5)       - wrapping c++ to python (c4 only)

In a recent check on Sept 17, 2019 for the Political Analysis replication, I was able to make an appropriate environment 
  with simply:

```
conda create --yes --channel conda-forge --name PAenv \
      matplotlib seaborn geopandas pandas numpy scipy \
      scikit-learn jellyfish psycopg2 jupyter
```

Because this does not freeze versions, this will not necessarily remain valid permanently.
GDAL's linking is notoriousy unstable in Anaconda, so the following may be more reliable:

```
conda env create --name pa_rep -f pa_rep.yaml
```

where `pa_rep.yaml` is the file included in the base of this directory.

Most of these scripts are implemented as python notebooks.
This allows intermediate outputs are comments to be displayed inline.
Jupyter is a standard python package, included with Anaconda installs.
It can be launched on a remote machine, and code can be examined
   and run through the web browser.
If users do not have it installed, and do not wish to install it,
  there are a number of docker containers for threading it directly to users' browser, as for instance, here:

https://hub.docker.com/r/continuumio/anaconda3/

Users preferring python notebooks can always convert them to python:
```
jupyter nbconvert --to python <notebook>.ipynb
```

or execute them directly from the command line:
```
jupyter nbconvert --to notebook --inplace --execute <notebook>.ipynb
```

In order to avoid timeouts and to specify which kernel/version of python to use, additional command line options can be specified:
```
jupyter nbconvert --ExecutePreprocessor.timeout=600 \
                  --ExecutePreprocessor.kernel_name=python3 \
                  --to notebook --inplace --execute main_tables.ipynb
```

## Computational Requirements

This analysis was performed on Ubuntu 16.04.6 LTS, with 16 GB of RAM and a modest CPU i5-6260U @ 2.9 GHz.  Run times for these settings are listed below.  None of this analysis was multithreaded.

The C4 map generation of course, was parallellized over thousands of cores, on AWS.  This is described on the [C4](http://github.com/JamesSaxon/c4) page.

## Data and Processing Code in this Repository

### Extracting state json files from bulk sim
The C4 program [outputs](https://github.com/JamesSaxon/c4#outputs) 
  json, csv, and geojson data for generated districts.
The csv files are simple tract to seat assignments,
  but the json and geojson contain fairly broad descriptions of district demographics,
    partisanship, spatial scores, population deviations, and so forth.

The full simulation came to around 30k maps, for each of the 11 core states;
  adding in the power diagrams for the full country, it is about 55 GB of raw data.
It was computationally inefficient to load each of these individually.
I therefore used the excellent json manipulation tool
  [`jq`](https://stedolan.github.io/jq/manual/)
  to combine these into single files.
The jq scripts are as follows: 
* The main jq manipulation is `extract.jq` which is run by `extract.sh`.
  This is the code that pulls out all voting returns,
  and creates a single output file for the main tables of the paper.
  These files are not included; instead, the files with final votes are (next section).
* For the power diagram districting and the implications for minorities, 
  minority VAP fractions for each district are extracted using `power_minority.jq` and `power_minority.sh`.
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

## Generating Tables and Figure

The scripts in this section are the ones
  that the replicator is expected to run.

A note on the data.
This was a bit of a balance to strike, 
  between the very low-density json data,
  and providing data that is _so_ processed
    that the replication itself is uninformative.
I have tried to skew towards the raw side (jq output) for the main results.
For the tables of the appendices,
  I did not always pre-consider how easy it would be to provide data
  in my responses to the referees.
The county-split table in particular is challenging
  without my local database and _all_ of the simulated data for MD, NC, and PA,
  and the table is only peripheral to the main argument.
For this table, I saved all partitionings of MD, NC, and PA.
Even compressed, these files cannot be uploaded to GitHub
  (though they could go to LFS, I don't want to pay for this),
  so they are only on the dataverse.

In all cases, the code to calculcate the derived data is provided,
  and the raw data can be made available on request.
The files are a bit large for the dataverse -- 
  _uncompressed_ the data come to around 55 GB;
  dropping geojson and tarring comes 20 GB.
I worked on compressing and reducing at some length, 
  and got the raw data under 10 GB.
So if necessary, I could split up states and get this onto the dataverse.
(The maximum file size is 2 GB, 
  and it seems that the dataset limit is 10 GB.)

### State-Level Seat Share Tables

* Code: [`main_tables.ipynb`](main_tables.ipynb)
* Run time: 1m45s real, 1m41s user
* Output: 
  * LaTeX files for the tables are written to `tex/`.
    * Table 2 -- seat shares for Pennsylvania -- is `tex/PA_table.tex`.  
    * Appendix Tables A.1-A.9 -- seat shares for the other states -- are named analogously.  
    * Appendix D -- Texas senate table is `TX_senate_table.tex`.
    * Appendix F -- NC race table is `NC_race.tex`.
    * Table 3 -- competitiveness is `tex/competitiveness_table.tex`
  * These tables "place" histograms generated by the scripts,
    themselves written to `mini_hist/` (seat shares) and `comp_hist/` (competitiveness).
    The naming conventions should be transparent, and are consistent between the tables and figures.

If C4 is the main "product" of the paper,
  the competitiveness and party share tables 
  are its main _results_.
These tables are fundamentally really simple:
  histograms and means for a bunch of different elections
    and partitioning methods (enacted and simulated).
The code should thus be straightforward.
It loads the jq'ed data,
  runs the plotting functions and calculates means or integrals,
  and finally runs a bunch of string manipulation to prettify the LaTeX outputs.

Note that the tables are built in a slightly unorthodox way:
  each of the "mini-histograms" in the table is
  a figure within a cell.
So the output of the code is
  (a) a latex file and (b) lots of mini-histograms.

### Competitiveness Table
* Code: [`main_tables.ipynb`](main_tables.ipynb)
* (See above for total run time, and output locations.)

The competitiveness table runs in the same notebook as the seat shares,
  and is built in basically the same way.

### Minority Seat Share
* Code: [`minority_seats.ipynb`](minority_seats.ipynb)
* Run time: 12s real, 14.5s user
* Output: Figure 4
  * `paper_figs/black_vap_representation.pdf`
  * `paper_figs/hispanic_vap_representation.pdf`

Minority seat shares first loads jq'ed power diagram simulation data.
It orders the seat shares and weights to create the "survival function"
  of minority seats against minority share.
It then runs a simple probit model of minority representation.
All of this should be clear from the code and its comments.

This model relies on the ethnic/racial identities of Representatives
  from the House's official history pages for [Blacks](http://history.house.gov/Exhibitions-and-Publications/BAIC/Historical-Data/Black-American-Representatives-and-Senators-by-Congress/)
  and [Hispanics](http://history.house.gov/Exhibitions-and-Publications/HAIC/Historical-Data/Hispanic-American-Representatives,-Senators,-Delegates,-and-Resident-Commissioners-by-Congress/),
  and analogous sources from the Press Gallery ([Blacks])https://pressgallery.house.gov/member-data/demographics/african-americans), [Hispanics](https://pressgallery.house.gov/member-data/demographics/hispanic-americans)).
There are a few discrepancies between these -- members in the press gallery but not in the History page, when I retrived it.
Those members are: Trent Franks R-AZ-8. John Garamendi D-CA-3, and Brian Mast R-FL-18.
I have include all three of them.

As for replacements, Chaka Fattah was replaced by Dwight Evants, and Xavier Becerra was replaced by Jimmy Gomez.
In each case, the ethnicity/race stayed the same, and so I included the districts just once.

Pew Research puts out pretty much the
  [same article on diversity](http://www.pewresearch.org/fact-tank/2017/01/24/115th-congress-sets-new-high-for-racial-ethnic-diversity/),
  for each Congress, which is also relevant.

### Pennsylania Plans
* Code: [`geomap.py`](geomap.py)
* Run time: 4.5s real, 4.4s user
* Figure 3: one pdf for each compactness method -- is written to `paper_figs/pa_ex/*pdf`

This code simply plots the "representative" plots of Pennsylvania,
  which are stored `data/pa_ex/`, based on seed 280.
The code simply plots the state colorfully,
  superimposes a circle over Pittsburgh,
  and writes the output to `paper_figs/pa_ex/`,

### Appendix: The Least Compact Districts 
* Code: [`cd_printer.py`](cd_printer.py)
* Run time: 5.1s real / 4.9s user
* <b>Not for replication</b>: this script is dependent on a private database.  It will not run, but it is simply a convenience script, to print a clean PDF of any US Congressional district.
* Outputs: Figure I.1 written to `paper_figs/bad_districts/`, of the form `[usps]_[cd].pdf`, e.g., `fl_5.pdf` for Florida's 5th congressional district.

This is a convenience script to plot
  a clean pdf of any US Congressional District.
These can be looked up with the Census, on Wikipedia, or on any number of Congressional widgets.
There is no "replication" to do.

The choice of districts is based on the PCA described below.

### Appendix: North Carolina Race 

* Code: [`app_post_selection_nc.ipynb`](app_post_selection_nc.ipynb)
* Run time: 2.8 s real / 2.2s user
* Outputs: Figure F: `paper_figs/nc_dseats_black_shares.pdf`

This code simply uses the most and second-most Black districts,
  as extracted with the jq files.
This creates violin plots of Black fraction,
  for the discussion of post-selection or post-analysis of 
  and "ancillary" properties of simulated districts.

### Appendix: County Splits

* Code: [`app_county_splits.ipynb`](app_county_splits.ipynb)
* Run time: approx 40 minutes.  As noted, this is by far the most awkward to replicate, and it is only tangential to the argument...
* Output: Table E: `tex/splits_table.tex` and corresponding mini histograms, written to `splits/*pdf`.

The question here is whether the compactness-based simulation
  "endangers" other districting principles.
To evaluate this, I calculated the 
  likelihood of a county resident having the same representative
  as someone else in his or her county.
This is corrected for the fact that some counties
  contain many Congressional districts, and in this case, 
  co-representation is of course impossible.
The evaluated quantity is thus pop(CD ∩ county) / min(pop(CD), pop(county)).

I prepared this table in reponse to a reviewer's question,
  and it is the one table where I did not preconsider the ease of replication.
I have therefore copied a lot of data to allow this one to be calculated from the dataverse files.
Sorry that it is somewhat slow and memory intensive.

### Appendix: PCA of Historic Seats

* Code: [`app_historic_interp_corr.ipynb`](app_historic_interp_corr.ipynb)
* Run time: 6.7s real / 7.5s user
* Output: Figure I.2: 
  * Pearson: `paper_figs/historic_correlation.pdf`
  * Spearman: `paper_figs/historic_spearmans.pdf`


This appendix considers a PCA of districts
  from multimember states over the past three districting cycles (107th, 11th, and 114th Congresses).
The aim is to understand the degree to which _compactness_ is a single concept
  (think: Massey/Denton segregation, but for shapes).
The result is that between 60 and 70% of the loading is on a single dimension,
  though this would of course change as a function of the sample
  and the specific compactness definitions included in the first place.

The relevance to the main argument is that interpersonal distance
  (which is closely related to the power diagram algorithm),
  is strongly correlated with the first component.
This is part of the justification for using it for the whole country,
  in the analysis of the impacts on minority representation.
  
This procedure is also used to identify the "worst" 5 districts of the
  114th Congress (based on the first component of the PCA),
  which are also plotted.


