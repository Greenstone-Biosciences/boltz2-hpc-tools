#!/bin/bash
#This script is used to analyze the docking folder output.  Last modified 05/28/2024 by J.MFing.Leitz

#SBATCH --job-name=AnalyzeDocking
#SBATCH --error=AnalD_%A_%a.err
#SBATCH --output=AnalD_%A_%a.out
#SBATCH --cpus-per-task=5
#SBATCH --ntasks=1

Input=$1 #Input folder containing pdbqts
Output=$2 #Output folder. Output is single file
fname=$3 #Name of file
QEDfile=$4 #File containing QED measure.
OrigSmi=$5
## argument parsing

#while getopts r:t:n:c: flag
#do
#    case "${flag}" in
#        i) Input=${OPTARG};;
#        o) Output=${OPTARG};;
#        n) fname=${OPTARG};;
#        q) QEDfile=${OPTARG};;
#    esac
#done

if [[ -z $fname ]]; then fname='Docking_Results_'; echo 'saving as '$fname; else echo 'saving as '$fname; fi
if [[ ! -d $Output ]]; then mkdir $Output; fi

if compgen -G "${Output}/$fname"* > /dev/null; 
then
	N=$( ls $Output/$fname* | sort -V | tail -n1 | rev | cut -d'_' -f 1 | rev);
	next=$(($N + 1));
else
	next=1
fi

for file in $Input/*.pdbqt; 
do 
energy=$(sed '2q;d' $file | awk '{ print $4 }'); 
name=$(sed '7q;d' $file | awk '{ print $4 }'); 
model=$(grep 'MODEL' $file | wc -l); 
ave=$(grep 'VINA' $file | awk '{ print $4 }' | awk '{ sum+=$1 }END { print sum/NR }'); 
QED=$(grep -w "${name}" $QEDfile | awk '{ print $3 }')
smile=$(grep -w "${name}" $OrigSmi | awk '{ print $1 }')
echo $fname'_'$next
echo $name
echo $smile $name $energy $model $ave $QED $file >> $Output/$fname 
done
