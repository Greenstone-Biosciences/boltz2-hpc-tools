#!/bin/bash

#SBATCH --error=array_%A_Dock.err
#SBATCH --output=array_%A_Dock.out
#SBATCH --export=ALL
#SBATCH --export-file=/home/jleitz/Software/VerdeLeitz_scripts/SlurmScripts/ConvertNDock.sh
#SBATCH --export=HOME

LigandFile_orig=$1
OutDIR=$2
Exhaustiveness=$3
VPDIR=$4
JSplit=$5

export HOME='/home/jleitz'
export '$HOME/Software/VerdeLeitz_scripts/SlurmScripts/ConvertNDock.sh'
echo $PATH
source /home/jleitz/env/bin/activate
export /home/jleitz/Software/

if [[ ! -d $OutDIR ]]; then mkdir $OutDIR; fi
if [[ ! -d $OutDIR/Ligands ]]; then mkdir $OutDIR/Ligands; fi
if [[ ! -d $OutDIR/Docked ]]; then mkdir $OutDIR/Docked; fi

jobs=$(($SLURM_CPUS_PER_TASK/10))

## Split Ligand Smiles into groups

LFName=$(basename $LigandFile_orig)
split -dn $JSplit $LigandFile_orig $OutDIR/Ligands/$LFName.split
for filelist in $OutDIR/Ligands/*split*;
do
	Stack=(${Stack[@]} "$filelist")
done

echo $SLURM_ARRAY_TASK_ID
LigandFile=${Stack[$SLURM_ARRAY_TASK_ID]}

for line in $(seq 1 $(cat $LigandFile | wc -l)); do 
	smile=$(sed "${line}q;d" $LigandFile | awk '{ print $1 }');
	name=$(sed "${line}q;d" $LigandFile | awk '{ print $2 }' | tr -d $'\r');
	CenterX=$(sed '2q;d' $VPDIR/*.coord | awk '{ print $2 }')
	CenterY=$(sed '2q;d' $VPDIR/*.coord | awk '{ print $3 }')
	CenterZ=$(sed '2q;d' $VPDIR/*.coord | awk '{ print $4 }')
	SizeX=$(sed '2q;d' $VPDIR/*.coord | awk '{ print $5 }')
        SizeY=$(sed '2q;d' $VPDIR/*.coord | awk '{ print $6 }')
        SizeZ=$(sed '2q;d' $VPDIR/*.coord | awk '{ print $7 }')
	Recept=$VPDIR/*.pdbqt
	echo 'converting '$name' from '$LigandFile
	sem -j $jobs --id $$ -u timeout 30m  ~/Software/VerdeLeitz_scripts/SlurmScripts/ConvertNDock.sh $line $LigandFile $OutDIR $Exhaustiveness $CenterX $CenterY $CenterZ $SizeX $SizeY $SizeZ $Recept;
	wait
		
done
