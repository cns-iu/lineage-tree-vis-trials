#!/bin/bash
source constants.sh
set -ev

RUNNER=${RUNNER:-docker}
CWL_FILE=https://raw.githubusercontent.com/cns-iu/map4sci/main/map4sci.cwl

SCRIPT="./scripts/99x-run-all-datasets.sh"
if [ "$2" != "" ]; then
  SCRIPT="$2"
fi

if [ "$RUNNER" == "singularity" ]; then
  export DATASETS_DIR=$(pwd)/datasets
  export RAW_DATA_DIR=$(pwd)/raw-data
  export SITE_DIR=$(pwd)/site

  export NODE_OPTIONS=--max-old-space-size=58000
  export CURRENT_DATASET=$1
  export CURRENT_VERSION=${CURRENT_VERSION:-v1}

  singularity exec docker://ghcr.io/cns-iu/map4sci:main /bin/bash \
    -c "cd /workspace/data-processor && ${SCRIPT}"
elif [ "$RUNNER" == "cwl" ]; then
  cwl-runner $CWL_OPTS $CWL_FILE \
    --dataset $1 --version $CURRENT_VERSION \
    --datasets_dir ./datasets --rawdata_dir ./raw-data --site_dir ./site \
    --script_cmd "cd /workspace/data-processor && ${SCRIPT}"
else
  docker compose -f src/docker-compose.yml run --rm \
    -e CURRENT_DATASET=$1 -e CURRENT_VERSION=${CURRENT_VERSION} map4sci ${SCRIPT}
fi
