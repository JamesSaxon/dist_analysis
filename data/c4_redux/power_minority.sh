#!/bin/bash 

reduxdir="/home/jsaxon/proj/dist_analysis/data/c4_redux"
s3json="/media/jsaxon/brobdingnag/data/s3/res/json"

# output list too long, if expanding full paths of inputs
cd $s3json 

# jq is used to do the extract.  then round with awk and drop quotes with sed.
jq -r -f ${reduxdir}/power_minority.jq -c *_power_s26*.json |
   awk -F "," -v OFS=, '{ print($1,$2,sprintf("%.5f",$3),sprintf("%.5f",$4)); }' |
   sed "s/\"//g" >  $reduxdir/power_minority_redux.csv

jq -r -f ${reduxdir}/power_minority.jq -c *_power_s2[7-9]*.json |
   awk -F "," -v OFS=, '{ print($1,$2,sprintf("%.5f",$3),sprintf("%.5f",$4)); }' |
   sed "s/\"//g" >> $reduxdir/power_minority_redux.csv

## Append on the same thing, but with seeds in the "300s."
## Otherwise, we cannot list them all!
jq -r -f ${reduxdir}/power_minority.jq -c *_power_s3*.json |
   awk -F "," -v OFS=, '{ print($1,$2,sprintf("%.5f",$3),sprintf("%.5f",$4)); }' |
   sed "s/\"//g" >> $reduxdir/power_minority_redux.csv

cd $reduxdir

sort power_minority_redux.csv > tmp
mv tmp power_minority_redux.csv

