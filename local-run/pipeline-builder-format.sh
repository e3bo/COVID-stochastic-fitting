#!/usr/bin/env bash

set -e

fdt="${fdt:-2020-11-16}"

dvc run \
    -d format-simulations.R \
    -d weekly-forecast-simulations/$fdt \
    -o CEID-compart_mif_sq/$fdt-CEID-compart_mif_sq.csv \
    -o CEID-compart_mif_rn/$fdt-CEID-compart_mif_rn.csv \
    -o CEID-compart_mif_li/$fdt-CEID-compart_mif_li.csv \
    -n format-$fdt \
    --force \
    podman run -v "$PWD":/root -e fdt=$fdt -w /root --rm eamon/2020mif:v20201129 ./format-simulations.R