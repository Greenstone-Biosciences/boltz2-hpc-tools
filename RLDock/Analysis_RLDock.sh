#!/bin/bash
#Anal_RLDock.sh
#Written by: Jeremy Leitz
#This is a script to analyze RLDock Output data specifically the cluster.mol2 files
#Collects name based on filename, best binding score, average binding score, number of poses.

set -e #Exit on error

#Usage
show_usage() {
	echo "usage: $0 -i <input directory> [-s file_containing_smiles]  [-o <path_to_Output_directory>] [--overwrite]"
	exit ${1:-0}
}

#Default values
Output="${PWD}/RLDock_results.csv"
overwrite=false

#Parse arguments
while [[ $# -gt 0 ]]; do
	key="$1"
	case $key in
		-i|--inputDIR)
			inputDIR="$2"
			shift
			shift
			;;
		-o|--Output)
			Output="$2"
			shift
			shift
			;;
		-s|--smiles)
			smiles_file="$2"
			shift
			shift
			;;
		--overwrite)
			overwrite=true
			shift
			;;
		*)
			show_usage 1
			;;
	esac
done

#Figure out if Output is a directory or a file, make the directory if needed
if [[ -d "$Output" ]]; then
	#Output is an existing directory
	result="${Output}/RLDock_results.csv"
	echo "Writing default RLDock_results.csv to ${Output}"
elif [[ "$Output" == */ ]]; then
	#Path ends with /, treat as directory
	mkdir -p "$Output"
	result="${Output}/RLDock_results.csv"
	echo "Created directory and writing RLDock_results.csv to ${Output}"
else
	#Treat it as a full file path
	result="$Output"
	out_dir=$(dirname "$result")
	if [[ ! -d "$out_dir" ]]; then
		mkdir -p "$out_dir" 
		echo "Created Output directory: $out_dir"
	fi
	echo "Writing results to $result"
fi

#Check if the results file already exits
if [[ -f "$result" ]]; then
	echo "${result} already exists." 
	if [[ "${overwrite}" == true  ]]; 
	then
		echo "Overwritting $result"
	else 
		echo "To overwrite include the '--overwrite' flag"
		exit 1
	fi
else
	echo "${result} does not exist, creating results file..."
fi


#Make header for results
echo "name,smiles,best,average,pose_number" > "$result"

for file in $inputDIR/*_cluster.mol2;
do
	if [[ ! -f $file ]]; then
		echo "Warning: No _cluster.mol2 files were found in ${inputDIR}"
		break
	fi

	name=$(basename "$file" "_cluster.mol2")
	best=$(awk '/^# Total_Energy:/ {print $3; exit}' "$file")
	average=$(awk '/^# Total_Energy:/ {sum+=$3; count++} END {if(count>0) print sum/count}' "$file")
	pose_num=$(grep -c 'Total_Energy' $file)
	if [[ -n "$smiles_file" && -f "$smiles_file" ]]; then
		smile=$(grep -w "${name}" "$smiles_file" | awk '{ print $1}')
	else
		smile=""
	fi

	echo "${name},${smile},${best},${average},${pose_num}" >> $result
done

echo "All Done Boss!!!"
echo "Your results are in ${result}"
