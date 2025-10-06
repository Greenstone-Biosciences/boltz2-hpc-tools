#!/bin/bash
#SBATCH -J boltz2
#SBATCH -p defq-gpu           # or your GPU partition
##SBATCH --gres=gpu:a100:1   # or a30:1
#SBATCH -N 1
#SBATCH -c 8
##SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o boltz2.%j.out
#SBATCH -e boltz2.%j.err

#module purge
module load cuda/12.1            # match your cluster’s CUDA
# (optional) module load cudnn    # if you use it

Input=$1 

source ~/anaconda3/bin/activate boltz

# Sanity check
#nvidia-smi
python -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())"

# Run Boltz-2 affinity prediction
srun boltz predict $Input --use_msa_server

