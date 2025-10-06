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
conda activate dock

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

LigandFile=$LigandFile_orig

## This does the bulk of the work. Gets the smile, name and coordinates from the ligand file and the coords file, respectively.
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
	if [ ! -f $OutDIR/Docked/$name"_out.pdbqt" ]; then
	echo 'converting '$name' from '$LigandFile
	echo "sem -j $jobs --id $$ -u timeout 30m /data/Shared/Docking_Scripts/ConvertNDock.sh $line $LigandFile $OutDIR $Exhaustiveness $CenterX $CenterY $CenterZ $SizeX $SizeY $SizeZ $Recept"
	sem -j $jobs --id $$ -u timeout 30m /data/Shared/Docking_Scripts/ConvertNDock.sh $line $LigandFile $OutDIR $Exhaustiveness $CenterX $CenterY $CenterZ $SizeX $SizeY $SizeZ $Recept;
	else 
		continue
	fi

		
done
sem --id $$ wait

conda deactivate
