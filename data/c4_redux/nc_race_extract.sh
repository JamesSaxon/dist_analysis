#!/bin/bash 

reduxdir="/home/jsaxon/proj/dist_analysis/data/c4_redux"
s3json="/media/jsaxon/brobdingnag/data/s3/res/json"

# output list too long, if expanding full paths of inputs
cd $s3json 

# North Carolina race selection
jq -f ${reduxdir}/nc_race_extract.jq -c nc_*_s2[7-9]*.json nc_*_s3*.json nc_split_s001.json > $reduxdir/nc_race_redux.json

jq -f ${reduxdir}/nc_race_shares.jq -c nc_power_s2[7-9]*.json  nc_power_s3*.json  > $reduxdir/nc_power_race.json
jq -f ${reduxdir}/nc_race_shares.jq -c nc_polsby_s2[7-9]*.json nc_polsby_s3*.json > $reduxdir/nc_polsby_race.json

cd $reduxdir

