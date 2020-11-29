#!/usr/bin/env bash

dvc run -d ../data \
  -d ../code \
  -d ../output \
  -d ./weekly-forecast-prep.R \
  -o header \
  -n make-header-dir \
  podman run -v "$PWD/..":/root -e DISABLE_AUTH=true -w /root/local-run --rm eamon/2020mif:v20201129 Rscript ./weekly-forecast-prep.R