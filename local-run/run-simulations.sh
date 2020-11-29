#!/usr/bin/env bash

set -e

fdt="${fdt:-2020-11-16}"

fdt=$fdt state_index=1 ./weekly-forecast.R