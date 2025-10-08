#!/bin/bash
#
# valid_boltz.sh
#
# Description: This script is intended to filter and validate that ligands are binding to correct locations on a receptor or protein from the outputs of Boltz. The script will take .cif files, align them, find centers of ligands, and then compare ligand positions last. Outputting if any have a distance from the aligned binding pockets.

# Input: directory containing .cifs with ligands as an output from boltz2
# Usage: ./valid_boltz.sh -i input_dir
# Requires: boltz env and USalign
#
#
#
# Default Values
INPUT_DIR=""
OUTPUT_DIR="${OUTPUT_DIR:-./aligned_cifs}" # Set the default values to ./aligned_cifs

# Help message
shoe_help() {
	cat << EOF
Usage: ${0##*/} [OPTIONS]

Aligns .cifs and then checks centers of ligands for binding positions and then compares.
To be generally used with receptor .cifs with 'bound' ligands either through boltz2 pipeline or other.

Required Arguments:
	-i, --input DIR		Input directory containing .cifs (/your/path)

Optional Arguments:
	-o, --output DIR
	-h, --help			Show this help message and exit
	-v, --verbose			Verbose output with dates

Examples:
	${0##*/} -i /data/Shared/Docking_Scripts/my_unaligned_cifs
	${0##*/} -input /data/Shared/Docking_Scripts/my_unaligned_cifs

EOF
	exit 0
}

# Parse Arguments
while [[ $# -gt 0 ]]; do
       case $1 in
       		-i|--input)
			if [[ -z "$2" || "$2" == -* ]]; then
				echo "Error: -i|--input requires a value" >&2
				exit 1
			fi
			INPUT_DIR="$2"
			shift 2
			;;
		-o|--output)
			if [[ -z "$2" || "$2" == -* ]]; then
				echo "Error: -o|--output requires a directory path" >&2
				exit 1
			fi
			OUTPUT_DIR="$2"
			shift 2
			;;
		-h|--help)
			show_help
			;;
		*)
			echo "Error: Unknown option: $1" >&2
			echo "Use -h or --help for usage information" >&2
			exit 1
	esac
done



# Validate required arguments
if [[ -z "$INPUT_DIR" ]]; then
	echo "Error: Missing required argument -i|--input" >&2
	echo "Use -h or --help for usage information" >&2
	exit 1
fi


# Validate input directory exists
if [[ ! -d "$INPUT_DIR" ]]; then
	echo "Error: Input directory '$INPUT_DIR' not found" >&2
	exit 1
fi

# Validate for output directory, if not create one to put aligned .cifs
if [[ ! -d "$OUTPUT_DIR" ]]; then
	mkdir -p $OUTPUT_DIR || {
		echo "Error: Cannot create output directory 'OUTPUT_DIR'" >&2
		exit 1
	}
fi

# Validate if output directory is even writable
if [[ ! -w "$OUTPUT_DIR" ]]; then
	echo "Error: Output directory '$OUTPUT_DIR' is not writable" >&2
	exit 1
fi

# Check for .cif files
CIF_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -name "*.cif" | wc -l)
if [[ $CIF_COUNT -eq 0 ]]; then
	echo "Error: No .cif files found in '$INPUT_DIR'" >&2
	exit 1
fi

# Source the conda setup and activate the correct conda env
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate boltz


# Navigate to the input directory

cd "$INPUT_DIR"

# Print configuration
echo "========================================"
echo "Boltz2 Docking Analysis"
echo "========================================"
echo "Input directory:   $INPUT_DIR"
echo "Output directory:  $OUTPUT_DIR"
echo "CIF files found:   $CIF_COUNT"
echo "========================================"
echo ""

## (p1)start the meat of the alignment

