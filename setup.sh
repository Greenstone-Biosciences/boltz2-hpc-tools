#!/bin/bash
# boltz2-hpc-tools setup script.  Run this once after clonging the repository and following other instuctions. 
#Written by Jeremy Leitz (he is the one to blame!). 

read -p "Path to this repository [$PWD?]: " NUTZ_PATH
NUTZ_PATH=${NUTZ_PATH:-"$HOME/Software/boltz-hpc-tools"}

read -p "Conda environment name [boltz]: " BOLTZ_CONDA_ENV
BOLTZ_CONDA_ENV=${CONDA_ENV:-"boltz"}

read -p "AFFINITY_TEMPLATE $NUTZ_PATH/affinity_template.yaml" AFFINITY_TEMPLATE

# Replace in all scripts
for script in boltz_batch_prepare Nutz.sh Analyze_boltz; do
    sed -i "s|NUTZ_PATH=.*|NUTZ_PATH=\"$NUTZ_PATH\"|" $script
    sed -i "s|BOLTZ_CONDA_ENV=.*|BOLTZ_CONDA_ENV=\"$BOLTZ_CONDA_ENV\"|" $script
    sed -i "s|AFFINITY_TEMPLATE=.*|AFFINITY_TEMPLATE=\"$AFFINITY_TEMPLATE\"|" $script
done

echo "Setup complete!"
