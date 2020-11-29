#!/usr/bin/env bash

set -e

fdt="${fdt:-2020-11-16}"

dvc run \
    -d pull-current-output-from-date.R \
    -o archive/$fdt \
    -n archive-$fdt \
    --force \
    podman run -v "$PWD":/root -e fdt=$fdt -w /root --rm eamon/2020mif:v20201129 ./pull-current-output-from-date.R