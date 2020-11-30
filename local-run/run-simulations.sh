#!/usr/bin/env bash

set -e


# Problem states: Florida
fdt="${fdt:-2020-11-16}"
declare -a locs=(
  "Hawaii"
  "Delaware"
  "Connecticut"
  "Colorado"
  "California"
  "Arizona"
  "Alaska"
  "Alabama"
  "Arkansas"
  "Minnesota" 
  "Georgia"
                )
                
for l in "${locs[@]}"; do
  fdt=$fdt locname=$l ./weekly-forecast.R
done
