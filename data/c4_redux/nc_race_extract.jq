# jq -f nc_race_extract.jq -c /media/jsaxon/brobdingnag/data/s3/res/json/nc_*_s2[7-9]*.json /media/jsaxon/brobdingnag/data/s3/res/json/nc_*_s3*.json /media/jsaxon/brobdingnag/data//s3/res/json/nc_split_s001.json > nc_race_redux.json 

select([.Districts[] | {DemFrac : ([.Elections[].DemFrac] | add / length), BlackVAPFrac : .Populations.BlackVAPFrac, HispanicVAPFrac : .Populations.HispanicVAPFrac} | if (.DemFrac > 0.5) and (.BlackVAPFrac > 0.30) then 1 else 0 end] | add | . >= 2) | 
{
  UID : .UID, 
  PopulationDeviation : .PopulationDeviation, 
  Score : .Score, 
  DemSeats : (.Elections | with_entries(.value |= .DemSeats)), 
  RepFrac : [.Districts[] | .Elections[] | .RepFrac ],
  BlackSeats : ([.Districts[] | {DemFrac : ([.Elections[].DemFrac] | add / length), BlackVAPFrac : .Populations.BlackVAPFrac, HispanicVAPFrac : .Populations.HispanicVAPFrac} | if (.DemFrac > 0.5) and (.BlackVAPFrac > 0.3) then 1 else 0 end] | add)
}

