# jq -f nc_race_shares.jq -c /media/jsaxon/brobdingnag/data/s3/res/json/nc_power_s[23]*.json > nc_power_race.json

{
  DemSeats : .Elections["2012"].DemSeats,
  BlackShare1 : [.Districts[] | select(.Elections["2012"].Party == "D") | .Populations.BlackVAPFrac] | sort | reverse | .[0],
  BlackShare2 : [.Districts[] | select(.Elections["2012"].Party == "D") | .Populations.BlackVAPFrac] | sort | reverse | .[1], 
} 

