#!/bin/bash
# boltz2-hpc-tools setup script.
# Run this once after cloning the repository and following the other instructions.
# Written by Jeremy Leitz @ Greenstone Biosciences
# https://github.com/Greenstone-Biosciences/boltz2-hpc-tools

read -p "Path to this repository [$PWD]: " NUTZ_PATH
NUTZ_PATH=${NUTZ_PATH:-"$PWD"}

read -p "Conda environment name [boltz]: " BOLTZ_CONDA_ENV
BOLTZ_CONDA_ENV=${BOLTZ_CONDA_ENV:-"boltz"}

# Replace placeholders in all scripts
for script in prep_batch_boltz Nutz.sh prep_inverse_screen InverseNutz.sh analyze_batch_boltz; do
    if [[ -f "$script" ]]; then
        sed -i "s|NUTZ_PATH|$NUTZ_PATH|g"           "$script"
        sed -i "s|BOLTZ_CONDA_ENV|$BOLTZ_CONDA_ENV|g" "$script"
    fi
done

echo "Setup complete!"
echo "  Repository path : $NUTZ_PATH"
echo "  Conda env       : $BOLTZ_CONDA_ENV"
