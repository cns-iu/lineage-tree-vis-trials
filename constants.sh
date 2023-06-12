shopt -s expand_aliases

export CURRENT_VERSION=v1

# Options are docker or cwl
export RUNNER=docker

# If using singularity, add --singularity to this environment variable
export CWL_OPTS=""

# Load environment
if [ -e env.sh ]; then
  source env.sh
fi

# Create directories
mkdir -p ./datasets ./raw-data ./site
