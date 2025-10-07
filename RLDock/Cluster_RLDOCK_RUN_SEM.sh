#!/bin/bash
## This script is to be run as an array.  
## This script is used to dock ligands to RNA

#SBATCH --error=array_%A_RNA_Dock.err
#SBATCH --output=array_%A_RNA_Dock.out
#SBATCH --export=ALL
#SBATCH --export-file=/data/Shared/Docking_Scripts/ConvertNDock_RNA.sh
#SBATCH --export=HOME

LigandFile_orig=$1
OutDIR=$2
Receptor=$3

export $HOME
export /data/Shared/Docking_Scripts/ConvertNDock_RNA.sh
echo $PATH
#source /home/jleitz/env/bin/activate
export /data/Shared/Docking_Scripts

source $HOME/anaconda3/etc/profile.d/conda.sh
conda activate dock

ml load gcc

#Make the Directories
if [[ ! -d $OutDIR ]]; then mkdir $OutDIR; fi
if [[ ! -d $OutDIR/Ligands ]]; then mkdir $OutDIR/Ligands; fi
if [[ ! -d $OutDIR/Docked ]]; then mkdir $OutDIR/Docked; fi
#Move into the docked directory in the output
cd $OutDIR/Docked

#If this is a synthemol file, convert it into a compatable format and rename to "SynthemolName.smi"
if [[ "${isSynMol}" == "Y" ]];
	then
	sed 's/,/\ /g' $LigandFile_orig > $OutDIR/Ligands/Synthemol_OG.smi
	len=$(cat $LigandFile_orig | wc -l)
	for i in $(seq 1 $len); do
		echo "SyntheMol_"$i >> $OutDIR/Ligands/SynthemolName.smi
	done
	cat $OutDIR/Ligands/Synthemol_OG.smi | awk '{ print $1 }' > $OutDIR/Ligands/SynthemolSmile.smi
	paste $OutDIR/Ligands/SynthemolSmile.smi $OutDIR/Ligands/SynthemolName.smi > $OutDIR/Ligands/SynthemolOut.smi
	LigandFile_orig=$OutDIR/Ligands/SynthemolOut.smi
fi

#How many jobs to run, assumes 10cpu per job.
jobs=$(($SLURM_CPUS_PER_TASK/10))

## Split Ligand Smiles into groups
LFName=$(basename $LigandFile_orig)
# Calculate lines per split file
total_lines=$(wc -l < "$LigandFile_orig")
lines_per_split=$(( (total_lines + $SLURM_ARRAY_TASK_COUNT - 1) / $SLURM_ARRAY_TASK_COUNT ))

# Split by lines, not bytes
split -l $lines_per_split -d "$LigandFile_orig" "$OutDIR/Ligands/$LFName.split."
wait

for filelist in $OutDIR/Ligands/*split*;
do
	Stack=(${Stack[@]} "$filelist")
done

echo $SLURM_ARRAY_TASK_ID
LigandFile=${Stack[$SLURM_ARRAY_TASK_ID]}

#Actual meat of the program, run ConvertNDock_RNA
for line in $(seq 1 $(cat $LigandFile | wc -l)); do	
	smile=$(sed "${line}q;d" $LigandFile | awk '{ print $1 }');
	name=$(sed "${line}q;d" $LigandFile | awk '{ print $2 }' | tr -d $'\r');
	#Check the file doesn't already exist
	if [ ! -f $OutDIR/Docked/$name"_out.mol2" ]; then
	echo 'converting '$name' from '$LigandFile
	echo "sem -j $jobs --id $$ -u timeout 30m /data/Shared/Docking_Scripts/ConvertNDock_RNA.sh $line $LigandFile $OutDIR $Receptor"
	sem -j $jobs --id $$ -u timeout 30m /data/Shared/Docking_Scripts/ConvertNDock_RNA.sh $line $LigandFile $OutDIR $Receptor;
	else 
		continue
	fi

		
done
sem --id $$ wait

conda deactivate
