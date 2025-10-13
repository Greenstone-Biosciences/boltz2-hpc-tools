#!/bin/bash
# boltz2-hpc-tools setup script.  Run this once after clonging the repository and following other instuctions. 
#Written by Jeremy Leitz (he is the one to blame!). 

read -p "Path to Boltz installation [$HOME/Software/boltz]: " BOLTZ_PATH
BOLTZ_PATH=${BOLTZ_PATH:-"$HOME/Software/boltz"}

read -p "Conda environment name [boltz]: " CONDA_ENV
CONDA_ENV=${CONDA_ENV:-"boltz"}

read -p "Slurm partition [defq-gpu]: " PARTITION
PARTITION=${PARTITION:-"defq-gpu"}

# Replace in all scripts
for script in boltz_batch_prepare Nutz.sh Analyze_boltz; do
    sed -i "s|BOLTZ_PATH=.*|BOLTZ_PATH=\"$BOLTZ_PATH\"|" $script
    sed -i "s|BOLTZ_CONDA_ENV=.*|BOLTZ_CONDA_ENV=\"$CONDA_ENV\"|" $script
    sed -i "s|SLURM_PARTITION=.*|SLURM_PARTITION=\"$PARTITION\"|" $script
done

echo "Setup complete!"
