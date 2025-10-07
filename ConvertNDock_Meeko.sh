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

smile=$(sed "${line}q;d" $LigandFile | awk '{ print $1 }');
name=$(sed "${line}q;d" $LigandFile | awk '{ print $2 }' | tr -d $'\r');
dt=$(date '+d%/%m/%Y %H:%M:%S')
nodename=$(hostname)
echo "running "$name 'at ' $dt 'on ' $nodename >> $OutDIR/Ligands/Log
python3 ~/Software/VerdeLeitz_scripts/SmileTo3D_RDKIT_WithQED.py $smile $name $OutDIR/Ligands
wait;
echo 'done converting '$name
mk_prepare_ligand.py -i $OutDIR/Ligands/$name.mol -o $OutDIR/Ligands/$name.pdbqt
wait;
echo 'done converting '$name' to pdbqt';
echo 'running in Vina'
/home/jleitz/Software/vina_1.2.3_linux_x86_64 --receptor $Recept --ligand $OutDIR/Ligands/$name'.pdbqt' --out $OutDIR/Docked/$name'_out.pdbqt' --cpu 10 --exhaustiveness $Exhaustiveness --center_x $CenterX --center_y $CenterY --center_z $CenterZ --size_x $SizeX --size_y $SizeY --size_z $SizeZ --num_modes 10; 

