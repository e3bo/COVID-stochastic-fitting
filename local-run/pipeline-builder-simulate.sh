#!/usr/bin/env bash

set -e

fdt="${fdt:-2020-11-16}"

dvc run \
    -d ../code \
    -d ./header \
    -d weekly-forecast.R \
    -d archive/$fdt \
    -o weekly-forecast-simulations/$fdt \
    -n simulate-$fdt \
    --force \
    podman run -v "$PWD/..":/root -e fdt=$fdt -w /root/local-run --rm eamon/2020mif:v20201129 bash -x ./run-simulations.sh