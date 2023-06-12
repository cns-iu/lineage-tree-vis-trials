#!/bin/bash
source constants.sh
set -ev

# Grab the "first" dataset
dataset=$(ls ./datasets | head -1)

time ./src/map4sci.sh $dataset ./scripts/99x-run-all-datasets.sh
