#!/usr/bin/env bash

set -e

fdt="${fdt:-2020-11-16}"

for i in $(seq 25 50); do
  if [[ "$i" == '2' || "$i" == '4' || "$i" == '21' || "$i" == '26' || "$i" == '27' || "$i" == '39' || "$i" == '48' || "$i" == '50' ]]; then
    continue
  fi
  fdt=$fdt state_index=$i ./weekly-forecast.R
done
