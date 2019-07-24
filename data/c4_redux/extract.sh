#!/bin/bash 

reduxdir="/home/jsaxon/proj/dist_analysis/data/c4_redux"
s3json="/media/jsaxon/brobdingnag/data/s3/res/json"

# output list too long, if expanding full paths of inputs
cd $s3json 

# json extraction for each state, using jq.
# It is unworkably inefficient to open each file in python...
for st in fl il la md mn nc pa tx tn va wi; do 

  jq -f ${reduxdir}/extract.jq -c ${st}_*_s2[7-9]*.json ${st}_*_s3*.json ${st}_split_s001.json > $reduxdir/${st}_redux.json

done 

cd $reduxdir

