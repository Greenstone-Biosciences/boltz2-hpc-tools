#!/bin/bash
## This script is to be run as an array.  
## Prior to running, you need in create a conda environment 'dock' using the environment (.yaml) file in /data/Shared/Docking_Scripts/ 
#SBATCH --error=array_%A_Dock.err
#SBATCH --output=array_%A_Dock.out
#SBATCH --export=ALL
#SBATCH --export-file=/data/Shared/Docking_Scripts/ConvertNDock.sh
#SBATCH --export=HOME

LigandFile_orig=$1
OutDIR=$2
Exhaustiveness=$3
VPDIR=$4
isSynMol=$5 # "Y/N"

export $HOME
export /data/Shared/Docking_Scripts/ConvertNDock.sh
echo $PATH
#source /home/jleitz/env/bin/activate
export /data/Shared/Docking_Scripts

source $HOME/anaconda3/etc/profile.d/conda.sh
#conda activate dock
conda activate pdbfixer310

if [[ ! -d $OutDIR ]]; then mkdir $OutDIR; fi
if [[ ! -d $OutDIR/Ligands ]]; then mkdir $OutDIR/Ligands; fi
if [[ ! -d $OutDIR/Docked ]]; then mkdir $OutDIR/Docked; fi
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

jobs=$(($SLURM_CPUS_PER_TASK/10))
##### IF array job or not #####

if [[ -z "$SLURM_ARRAY_TASK_ID" ]]; then
	echo "Running on single node."
	LigandFile="$LigandFile_orig"
else
	echo "Running in array mode. Splitting ligand file..."
	LFName=$(basename $LigandFile_orig)
	split -dn $SLURM_ARRAY_TASK_COUNT $LigandFile_orig $OutDIR/Ligands/$LFName.split
	wait

	for filelist in $OutDIR/Ligands/*split*;
	do
		Stack=(${Stack[@]} "$filelist")
	done

	echo "SLURM ARRAY TASK ID = $SLURM_ARRAY_TASK_ID"
	LigandFile=${Stack[$SLURM_ARRAY_TASK_ID]}
fi

for line in $(seq 1 $(cat $LigandFile | wc -l)); do	
	smile=$(sed "${line}q;d" $LigandFile | awk '{ print $1 }');
	if [[ "${smile,,}" != "smiles" ]]; then
		name=$(sed "${line}q;d" $LigandFile | awk '{ print $2 }' | tr -d $'\r');
		CenterX=$(sed '2q;d' $VPDIR/*.coord | awk '{ print $2 }')
		CenterY=$(sed '2q;d' $VPDIR/*.coord | awk '{ print $3 }')
		CenterZ=$(sed '2q;d' $VPDIR/*.coord | awk '{ print $4 }')
		SizeX=$(sed '2q;d' $VPDIR/*.coord | awk '{ print $5 }')
		SizeY=$(sed '2q;d' $VPDIR/*.coord | awk '{ print $6 }')
        	SizeZ=$(sed '2q;d' $VPDIR/*.coord | awk '{ print $7 }')
		Recept=$VPDIR/*.pdbqt
		if [ ! -f $OutDIR/Docked/$name"_out.pdbqt" ]; then
			echo 'converting '$name' from '$LigandFile
			echo "sem -j $jobs --id $$ -u timeout 30m /data/Shared/Docking_Scripts/ConvertNDock.sh $line $LigandFile $OutDIR $Exhaustiveness $CenterX $CenterY $CenterZ $SizeX $SizeY $SizeZ $Recept"
			sem -j $jobs --id $$ -u timeout 30m ~/Software/VerdeLeitz_scripts/SlurmScripts/ConvertNDock_Meeko.sh $line $LigandFile $OutDIR $Exhaustiveness $CenterX $CenterY $CenterZ $SizeX $SizeY $SizeZ $Recept;
		else 
			continue
		fi
	else 
		continue
	fi

		
done
sem --id $$ wait

conda deactivate
