#!/usr/bin/env bash

set -e

dvc run \
    -d CEID-compart_mif_li \
    -d CEID-compart_mif_rn \
    -d CEID-compart_mif_sq \
    -d evaluate-forecasts.R \
    -o evaluate-forecasts.html \
    -n evaluate \
    --force \
    podman run -v "$PWD":/root -w /root --rm eamon/2020mif:v20201129 'Rscript -e "knitr::spin(\"evaluate-forecasts.R\")" '