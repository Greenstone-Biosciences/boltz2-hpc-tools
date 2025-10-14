#!/bin/bash
# boltz2-hpc-tools setup script.  Run this once after clonging the repository and following other instuctions. 
#Written by Jeremy Leitz (he is the one to blame!). 

read -p "Path to this repository [$PWD?]: " NUTZ_PATH
NUTZ_PATH=${NUTZ_PATH:-"$HOME/Software/boltz-hpc-tools"}

read -p "Conda environment name [boltz]: " BOLTZ_CONDA_ENV
BOLTZ_CONDA_ENV=${BOLTZ_CONDA_ENV:-"boltz"}

read -p "Path to the affinity_template yaml [$PWD]: " AFFINITY_TEMPLATE
AFFINITY_TEMPLATE=${AFFINITY_TEMPLATE:-"$HOME/Software/boltz-hpc-tools"}

# Replace in all scripts
for script in prep_batch_boltz Nutz.sh analyze_batch_boltz; do
    sed -i "s|NUTZ_PATH|$NUTZ_PATH|g" $script
    sed -i "s|BOLTZ_CONDA_ENV|$BOLTZ_CONDA_ENV|g" $script
    sed -i "s|AFFINITY_TEMPLATE|$AFFINITY_TEMPLATE|g" $script
done

echo "Setup complete!"
