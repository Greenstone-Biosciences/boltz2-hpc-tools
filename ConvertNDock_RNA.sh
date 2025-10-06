#!/bin/bash

#This script is to convert a single molecule from smiles to pdbqt and dock it using autodockvina.

line=$1
LigandFile=$2
OutDIR=$3
receptor=$4

smile=$(sed "${line}q;d" $LigandFile | awk '{ print $1 }');
name=$(sed "${line}q;d" $LigandFile | awk '{ print $2 }' | tr -d $'\r');
dt=$(date '+d%/%m/%Y %H:%M:%S')
nodename=$(hostname)
echo "running "$name 'at ' $dt 'on ' $nodename >> $OutDIR/Ligands/Log
python3 /data/Shared/Docking_Scripts/SmileTo3D_RDKIT_WithQED.py $smile $name $OutDIR/Ligands
wait;
echo 'done converting '$name
obabel $OutDIR/Ligands/$name".mol" -omol2 -O $OutDIR/Ligands/$name".mol2" -h
wait;
echo 'done converting '$name' to mol2';
echo 'running in RLDOCK'
/home/jleitz/Software/RLDOCK/RLDOCK/bin/RLDOCK -i $receptor -l $OutDIR/Ligands/$name".mol2" -o $name -c 10 -n 10 -s /home/jleitz/Software/RLDOCK/RLDOCK/src/sphere.dat
