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
# Default values
INPUT_DIR=""
OUTPUT_DIR=""
VERBOSE=false

# Help message
show_help() {
	cat << EOF
Usage: ${0##*/} [OPTIONS]

Analyzes Boltz2 docking results for kinases with cofolded ligands by:
  1. Aligning all protein structures using US-align
  2. Parsing aligned .cif files to extract ligand (HETATM) coordinates
  3. Calculating ligand centroids
  4. Computing pairwise distances between ligand positions

Required Arguments:
    -i, --input DIR           Directory containing .cif files from Boltz2

Optional Arguments:
    -o, --output DIR          Output directory (default: ./analysis_output)
    -v, --verbose             Verbose output
    -h, --help                Show this help message

Requirements:
    - US-align must be in PATH
    - Python 3 with NumPy

Input Requirements:
    - Directory must contain .cif files with kinase + ligand (HETATM)
    - All structures should be of similar kinases with different ligand poses

Output:
    - aligned_cifs/           Aligned .cif files (input files unchanged)
    - centroids.csv           Ligand centroid coordinates
    - distances.csv           Pairwise distances between ligands
    - binding_pocket.txt      Predicted binding pocket center
    - analysis_summary.txt    Summary report

Examples:
    ${0##*/} -i ./boltz2_output/
    ${0##*/} --input ./kinase_results/ --output ./analysis/

EOF
    exit 0
}


# Parse Arguments
while [[ $# -gt 0 ]]; do
       case $1 in
       		-i|--input)
			if [[ -z "$2" || "$2" == -* ]]; then
				echo "Error: -i|--input requires a directory path" >&2
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
		-v|--verbose)
			VERBOSE=true
			shift
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


# Source the conda setup and activate the correct conda env
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate boltz


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


# Check for .cif files
CIF_FILES=($(find "$INPUT_DIR" -maxdepth 1 -name "*.cif" -type f | sort))
CIF_COUNT=${#CIF_FILES[@]}

if [[ $CIF_COUNT -eq 0 ]]; then
	echo "Error: No .cif files found in '$INPUT_DIR'" >&2
	exit 1
fi

# Set default output directory if not provided
OUTPUT_DIR="${OUTPUT_DIR:-./analysis_output}"

# Create output directory
if [[ ! -d "$OUTPUT_DIR" ]]; then
	mkdir -p "$OUTPUT_DIR" || {
		echo "Error: Cannot create output directory '$OUTPUT_DIR'" >&2
		exit 1
	}
fi

# Verify output directory is writable
if ! touch "$OUTPUT_DIR/.write_test" 2>/dev/null; then
	echo "Error: Cannot write to output directory '$OUTPUT_DIR'" >&2
        exit 1
fi

rm "$OUTPUT_DIR/.write_test"

# Create subdirectories
mkdir -p "$OUTPUT_DIR/aligned_cifs"
mkdir -p "$OUTPUT_DIR/usalign_logs"

# Check if US-align is available
if ! command -v USalign &> /dev/null; then
	echo "Error: US-align not found in PATH" >&2
	echo "Please install US-align or ensure it's in your conda environment" >&2
	exit 1
fi

# Log function
log() {
	if [[ "$VERBOSE" == true ]]; then
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
	fi
}



# Print configuration
echo "========================================"
echo "Boltz2 Docking Analysis"
echo "========================================"
echo "Input directory:   $INPUT_DIR"
echo "Output directory:  $OUTPUT_DIR"
echo "CIF files found:   $CIF_COUNT"
echo "========================================"
echo ""

log "Starting analysis pipeline of .cifs..."


## Step 1:
