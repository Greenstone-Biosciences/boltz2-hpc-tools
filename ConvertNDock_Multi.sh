#!/bin/bash

#This script is to convert a single molecule from smiles to pdbqt and dock it using autodockvina.

line=$1
LigandFile=$2
OutDIR=$3
Exhaustiveness=$4
CenterX=$5
CenterY=$6
CenterZ=$7
SizeX=$8
SizeY=$9
SizeZ=${10}
Recept=${11}
SiteDIR=${12}

check=$(cat -v $LigandFile | grep '\^M' | wc -l)
if [[ $check -gt 0 ]]; then sed -i -e "s/\r//g" $LigandFile; fi

smile=$(sed "${line}q;d" $LigandFile | awk '{ print $1 }');
name=$(sed "${line}q;d" $LigandFile | awk '{ print $2 }' | tr -d $'\r');
echo "running "$name
python3 ~/Software/VerdeLeitz_scripts/SmileTo3D_RDKIT.py $smile $name $OutDIR/Ligands
wait;
echo 'done converting '$name
obabel $OutDIR/Ligands/$name".mol" -opdbqt -O $OutDIR/Ligands/$name".pdbqt" -h
wait;
echo 'done converting '$name' to pdbqt';
echo 'running in Vina'
/home/jleitz/Software/vina_1.2.3_linux_x86_64 --receptor $Recept --ligand $OutDIR/Ligands/$name".pdbqt" $OutDIR/Ligands/$name".pdbqt" --out $OutDIR/Docked/$name'_out.pdbqt' --cpu 20 --exhaustiveness $Exhaustiveness --center_x $CenterX --center_y $CenterY --center_z $CenterZ --size_x $SizeX --size_y $SizeY --size_z $SizeZ --num_modes 10; 

