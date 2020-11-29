#!/usr/bin/env bash

set -e

fdt="${fdt:-2020-11-16}"
declare -a locs=("Minnesota" 
                "Georgia"
                )
                
for l in "${locs[@]}"; do
  fdt=$fdt locname=$l ./weekly-forecast.R
done
