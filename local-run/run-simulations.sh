#!/usr/bin/env bash

set -e


#  Lousianna
fdt="${fdt:-2020-11-16}"
declare -a locs=(
"Washington"
"New Jersey"
"Nebraska"
"Louisiana"
"Massachusetts"
"Florida"
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
