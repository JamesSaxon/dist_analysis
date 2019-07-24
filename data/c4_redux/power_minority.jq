# select plans with no deviations larger than 2%
# Note the state USPS and the simulation UID.
# Then grab the Black and Hispanic VAP Fractions per district.
# Output as a CSV.

select(.PopulationDeviation < 0.02) | 
(.UID | gsub(".*power/"; "")) as $UX | 
.USPS as $USPS | 
(.Districts[].Populations | 
 [$USPS, $UX, .BlackVAPFrac, .HispanicVAPFrac]) | 
@csv
