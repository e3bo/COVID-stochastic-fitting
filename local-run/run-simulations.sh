#!/usr/bin/env bash

set -e


# Problem states: Florida Massachusetts Lousianna "Nebraska" "New Jersey" Washington
fdt="${fdt:-2020-11-16}"
declare -a locs=(
"Wisconsin"
"West Virginia"
"Wyoming"
  "Missouri"
  "Minnesota"
  "Michigan"
  "Maine"
  "Maryland"
  "Kentucky"
  "Kansas"
  "Indiana"
  "Illinois"
  "Idaho"
  "Iowa"
  "Hawaii"
  "Delaware"
  "Connecticut"
  "Colorado"
  "California"
  "Arizona"
  "Alaska"
  "Alabama"
  "Arkansas"
  "Georgia"
  "Mississippi"
"Montana"
"North Carolina"
"North Dakota"
"New Hampshire"
"New Mexico"
"Nevada"
"New York"
"Ohio"
"Oklahoma"
"Oregon"
"Pennsylvania"
"Rhode Island"
"South Carolina"
"South Dakota"
"Tennessee"
"Texas"
"Utah"
"Virginia"
"Vermont"
)
                
for l in "${locs[@]}"; do
  fdt=$fdt locname=$l ./weekly-forecast.R
done
