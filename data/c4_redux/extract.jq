# jq -f extract.jq -c /media/jsaxon/brobdingnag/data/s3/res/json/tx_*_s2[7-9]*.json ../s3/res/json/tx_*_s3*.json ../s3/res/json/tx_split_s001.json

{
  UID : .UID, 
  PopulationDeviation : .PopulationDeviation, 
  Score : .Score, 
  DemSeats : (.Elections | with_entries(.value |= .DemSeats)), 
  RepFrac : [.Districts[] | .Elections[] | .RepFrac ]
}
