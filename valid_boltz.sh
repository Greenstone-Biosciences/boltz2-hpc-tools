#!/bin/bash
################################################################################
# Script Name:    valid_boltz.sh
# Description:    Align Boltz2 structures and calculate ligand centroids, binding pockets
# Usage:          ./valid_boltz.sh -i INPUT_DIR [-o OUTPUT_DIR]
# Written by:     Chris Yan
################################################################################

# Default values
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INPUT_DIR=""
OUTPUT_DIR=""
VERBOSE=false
CLUSTER_THRESHOLD="${CLUSTER_THRESHOLD:-5.0}"
DRY_RUN=false

# Help message
show_help() {
	cat << EOF
Usage: ${0##*/} [OPTIONS]

Aligns .cifs protein structures and calculates ligand centroids
(e.g.Boltz2 outputs cofolded single receptor with multiple ligands)
  1. Aligning all protein structures using gemmi
  2. Parsing aligned .cif files to extract ligand and calculate ligand centroids
  3. Calculated average position of all ligands and individual ligand distances
  4. Clusters ligand centroids to identify and assign binding pockets

Required Arguments:
    -i, --input DIR           Directory containing .cif files

Optional Arguments:
    -o, --output DIR          Output directory (default: ./analysis_output)
    -t, --threshold FLOAT     Clustering threshold in Å (default: 5.0)
                              Set via CLUSTER_THRESHOLD environment variable
    --dry-run                 Show what files would be processed without running analysis
    -v, --verbose             Verbose output
    -h, --help                Show this help message

Requirements:
    - gemmi Python package
    - bc calculator
    - identify_pockets.py (in same directory)
    - (boltz) conda environment
    - USalign

Input Requirements:
    - Directory must contain .cif files

Output:
    - aligned_cifs/           Aligned .cif files (protein + ligand)
    - ligand_centers.csv      Ligand center of mass coordinates
    - ligand_deviations.csv   Distance from average position
    - pocket_assignments.csv  Pocket IDs per ligand"
    - cluster_statistics.csv  Per-pocket cluster statistics"

Examples:
    ${0##*/} -i ./my_unaligned_cifs/
    ${0##*/} --input ./my_unaligned_cifs/ --output ./analyzed_.cifs/

EOF
    exit 0
}


# Parse Arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            # Input directory validation
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output)
            # Output directory
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threshold)
            # Clustering threshold in Angstroms
            CLUSTER_THRESHOLD="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
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
            ;;
    esac
done


# Source the conda setup and activate the correct conda env
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate boltz


# Validate input directory requirements
if [[ -z "$INPUT_DIR" ]]; then
	echo "Error: Missing required argument -i|--input" >&2
	echo "Use -h or --help for usage information" >&2
	exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
	echo "Error: Input directory '$INPUT_DIR' not found" >&2
	exit 1
fi


# Check for .cif files - conditional detection for Boltz2 output structures
# mapfile -t CIF_FILES < <(find "$INPUT_DIR" -maxdepth 4 -name "*.cif" -type f | sort)
if find "$INPUT_DIR" -type d -name "Ligand_yamls" -maxdepth 3 2>/dev/null | grep -q .; then
    echo "Detected Boltz2 output structure, using predictions pattern..."
    mapfile -t CIF_FILES < <(find "$INPUT_DIR" -path "*/predictions/*/*.cif" -type f | sort)
else
    echo "No Ligand_yamls folder detected, searching for all CIF files..."
    mapfile -t CIF_FILES < <(find "$INPUT_DIR" -name "*.cif" -type f | sort)
fi

CIF_COUNT=${#CIF_FILES[@]}

if [[ $CIF_COUNT -eq 0 ]]; then
	echo "Error: No .cif files found in '$INPUT_DIR'" >&2
	exit 1
fi

# Dry run: show what would be processed and exit
if [[ "$DRY_RUN" == true ]]; then
    echo ""
    echo "========================================"
    echo "DRY RUN MODE - No files will be processed"
    echo "========================================"
    echo "Would process $CIF_COUNT CIF files"
    echo ""
    echo "First 10 files:"
    printf '  %s\n' "${CIF_FILES[@]}" | head -10
    if [[ $CIF_COUNT -gt 10 ]]; then
        echo " ... and $((CIF_COUNT - 10)) more files"
    fi
    echo ""
    echo "Exiting (no processing performed)"
    exit 0
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
mkdir -p "$OUTPUT_DIR/logs"

# Check gemmi is available in environment
if ! python3 -c "import gemmi" 2>/dev/null; then
    echo "Error: gemmi not found. Install with: pip install gemmi" >&2
    exit 1
fi

# Check scikit-learn (for clustering)
if ! python3 -c "import sklearn" 2>/dev/null; then
    echo "Warning: scikit-learn not found. Pocket identification will be skipped." >&2
    echo "Install with: pip install scikit-learn" >&2
fi

# Check bc is available
if ! command -v bc &> /dev/null; then
    echo "Error: bc calculator not found" >&2
    exit 1
fi

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
echo ""


##################################################
## Step 1: Retrieve .cifs and alignment
##################################################
# Run alignment and extract ligand centers of mass
python3 << EOF
import gemmi
import sys
import os
from pathlib import Path

input_dir = "$INPUT_DIR"
output_dir = "$OUTPUT_DIR"
verbose = "$VERBOSE"

print("Step 1: Retrieve .cifs and alignment")

# Get all CIF files - match bash logic for Boltz2 structure detection
# DEPRECATED cif_files = sorted([str(f) for f in Path(input_dir).glob("*/predictions/*/*.cif")])
input_path = Path(input_dir)

# Check if Boltz2 directory structure exists (look for Ligand_yamls directory)
ligand_yamls_dirs = list(input_path.glob("**/Ligand_yamls"))

if ligand_yamls_dirs:
    # Boltz2 directory structure: search within Ligand_yamls for predictions pattern
    if verbose:
        print("Detected Boltz2 structure, using predicions pattern")
    cif_files = []
    for lyd in ligand_yamls_dirs:
        cif_files.extend(lyd.glob("*/predictions/*/*.cif"))
    cif_files = sorted([str(f) for f in cif_files])
else:
    # Flat directory structure: find any CIF files
    if verbose:
        print("No Ligand_yamls detected, searching for all CIF files")
    cif_files = sorted([str(f) for f in input_path.glob("**/*.cif")])

if not cif_files:
    print("Error: No .cif files found", file=sys.stderr)
    sys.exit(1)


# Use first file as reference
ref_file = cif_files[0]
ref_basename = os.path.basename(ref_file)
ref_name = os.path.splitext(ref_basename)[0]  # Remove .cif extension

# Load reference structure and extract protein polymer
print("Loading reference structure")
print(f"Reference: {ref_basename}\n")

# Load reference structure
ref_st = gemmi.read_structure(ref_file)
ref_model = ref_st[0]

# Get reference polymer (protein chain for alignment)
ref_polymer = None
for chain in ref_model:
	poly = chain.get_polymer()
	if poly:
		ref_polymer = poly
		if verbose:
			print(f"Found polymer in chain {chain.name}: {len(poly)} residues")
		break

# Save reference structure to output_dir
ref_output = os.path.join(output_dir, "aligned_cifs", f"{ref_name}_aligned.cif")
ref_st.make_mmcif_document().write_file(ref_output)
print(f"Wrote reference file to {ref_output}")


print(f"\n{'Structure':<50} {'RMSD (Å)'}")
print("-" * 65)
print(f"{ref_basename:<50} {'0.0000'}")

for mobile_file in cif_files[1:]:
	mobile_basename = os.path.basename(mobile_file)
	mobile_name = os.path.splitext(mobile_basename)[0]
#	print(mobile_basename)
#	print(mobile_name)

	try:
		mobile_st = gemmi.read_structure(mobile_file)
		mobile_model = mobile_st[0]
		mobile_polymer = None
		for chain in mobile_model:
			poly = chain.get_polymer()
			if poly:
				mobile_polymer = poly
				break

		if not mobile_polymer:
			print(f"{mobile_basename:<50} No polymer found", file=sys.stderr)
			continue

		# Main Calculation of mobile polymer superposition relative to reference polymer
		# ptype = gemmi.PolymerType.PeptideL
		sup = gemmi.calculate_superposition(
			ref_polymer,
			mobile_polymer,
			gemmi.PolymerType.PeptideL,
			gemmi.SupSelect.CaP
		)
	
		# The alignment step where it applies the changes to mobile polymer
		for chain in mobile_model:
			for residue in chain:
				for atom in residue:
					atom.pos = gemmi.Position(sup.transform.apply(atom.pos))
	
		print(f"{mobile_basename:<50} {sup.rmsd:.3f}")
#		print(f"Aligned {sup.count} matching CA atoms")
			
	
		# Writes out the new aligned .cif
		output_file = os.path.join(output_dir, "aligned_cifs", f"{mobile_name}_aligned.cif")
		mobile_st.make_mmcif_document().write_file(output_file)
	
	
	except Exception as e:
		print(f"{mobile_basename:<50} Error: {e}", file=sys.stderr)
		continue

print(f"\n✓ Alignment complete")


EOF

if [[ $? -ne 0 ]]; then
	echo "Error: Alignment failed :(" >&2
	exit 1
fi



##################################################
## Step 2:  Extract Ligand Centroids of Aligned .cifs
##################################################

echo ""
echo "Step 2: Extracting ligand centers of mass from aligned structures"
echo ""

CENTROID_FILE="$OUTPUT_DIR/ligand_centers.csv"
echo "Structure,X,Y,Z,N_atoms" > "$CENTROID_FILE"

printf "%-50s %35s %10s\n" "Structure" "Ligand Centroid (x, y, z)" "# Atoms"
echo "--------------------------------------------------------------------------------"

for CIF in "$OUTPUT_DIR/aligned_cifs"/*.cif; do
	BASENAME=$(basename "$CIF")

	# Extract HETATM lines for LIG1 ligand and get coordinates
	# Only works if col 11, 12, and 13 of .cif are the x, y, z
	COORDS=$(awk '
        /^ATOM/ || /^HETATM/ {
            # Check if this line contains LIG1
                if ($0 ~ /LIG1/) {
	                    # CIF format: fields are space-separated
	                    # Print columns 11 (x), 12 (y), 13 (z)
	                print $11, $12, $13
                }
        }
	' "$CIF")
    
	if [[ -z "$COORDS" ]]; then
	    echo "$BASENAME: No LIG1 ligand found" >&2
	     continue
	fi
	    
	# Calculate centroid using bc
	N_ATOMS=$(echo "$COORDS" | wc -l)
	SUM_X=0
	SUM_Y=0
	SUM_Z=0
	
	while read -r X Y Z; do
		SUM_X=$(echo "$SUM_X + $X" | bc)
		SUM_Y=$(echo "$SUM_Y + $Y" | bc)
		SUM_Z=$(echo "$SUM_Z + $Z" | bc)
	done <<< "$COORDS"
	
	CENT_X=$(echo "scale=3; $SUM_X / $N_ATOMS" | bc)
	CENT_Y=$(echo "scale=3; $SUM_Y / $N_ATOMS" | bc)
	CENT_Z=$(echo "scale=3; $SUM_Z / $N_ATOMS" | bc)
	
	# Save to file
	echo "$BASENAME,$CENT_X,$CENT_Y,$CENT_Z,$N_ATOMS" >> "$CENTROID_FILE"
	
	printf "%-50s (%8.3f, %8.3f, %8.3f) %10s\n" "$BASENAME" "$CENT_X" "$CENT_Y" "$CENT_Z" "$N_ATOMS"
done

echo ""
echo "✓ Ligand centers extracted and saved to ligand_centers.csv"



################################################################################
# STEP 3: CALCULATE AVERAGE LIGAND POSITION AND DEVIATIONS
################################################################################

echo ""
echo "Step 3: Analyzing ligand position variability"
echo ""


# Calculate average position (centroid of centroids)
SUM_X=0
SUM_Y=0
SUM_Z=0
COUNT=0

while IFS=',' read -r FILE X Y Z N; do
	[[ "$STRUCTURE" == "Structure" ]] && continue
	SUM_X=$(echo "$SUM_X + $X" | bc)
	SUM_Y=$(echo "$SUM_Y + $Y" | bc)
	SUM_Z=$(echo "$SUM_Z + $Z" | bc)
	COUNT=$((COUNT + 1))
done < "$CENTROID_FILE"


if [[ $COUNT -eq 0 ]]; then
	echo "Error: No ligand centroids found" >&2
	exit 1
fi


AVG_X=$(echo "scale=3; $SUM_X / $COUNT" | bc)
AVG_Y=$(echo "scale=3; $SUM_Y / $COUNT" | bc)
AVG_Z=$(echo "scale=3; $SUM_Z / $COUNT" | bc)



echo "Average ligand position (binding site center):"
echo "  X: $AVG_X Å"
echo "  Y: $AVG_Y Å"
echo "  Z: $AVG_Z Å"
echo "  Based on $COUNT ligand poses"
echo ""

# Calculate deviations from average
DEVIATIONS_FILE="$OUTPUT_DIR/ligand_deviations.csv"
echo "Structure,Distance_from_average,Direction_X,Direction_Y,Direction_Z" > "$DEVIATIONS_FILE"

printf "%-50s %15s %30s\n" "Structure" "Distance (Å)" "Direction Vector (x, y, z)"
echo "--------------------------------------------------------------------------------"


while IFS=',' read -r FILE X Y Z N; do
	[[ "$STRUCTURE" == "Structure" ]] && continue

	# Calculate distance from average (magnitude)
	DX=$(echo "$X - $AVG_X" | bc)
	DY=$(echo "$Y - $AVG_Y" | bc)
	DZ=$(echo "$Z - $AVG_Z" | bc)

	DIST=$(echo "scale=3; sqrt(($DX)^2 + ($DY)^2 + ($DZ)^2)" | bc)

	# Normalize direction vector (unit vector)
	if (( $(echo "$DIST > 0" | bc -l) )); then
		NORM_X=$(echo "scale=3; $DX / $DIST" | bc)
		NORM_Y=$(echo "scale=3; $DY / $DIST" | bc)
		NORM_Z=$(echo "scale=3; $DZ / $DIST" | bc)
	else
		NORM_X=0.000
		NORM_X=0.000
		NORM_X=0.000
	fi

	# Safe to file
	echo "$FILE,$DIST,$NORM_X,$NORM_Y,$NORM_Z" >> "$DEVIATIONS_FILE"

	printf "%-50s %15.3f (%6.3f, %6.3f, %6.3f)\n" "$FILE" "$DIST" "$NORM_X" "$NORM_Y" "$NORM_Z"

done < "$CENTROID_FILE"


echo ""
echo "✓ Deviations saved to ligand_deviations.csv"

################################################################################
# STEP 4: DETERMINE BINDING POCKETS BY CLUSTERING
################################################################################

echo ""
echo "Step 4: Identifying binding pockets via clustering"
echo ""

# Check if identify_pockets.py exists
POCKET_SCRIPT="$SCRIPT_DIR/identify_pockets.py"
if [[ ! -f "$POCKET_SCRIPT" ]]; then
    echo "Warning: $POCKET_SCRIPT not found in current directory" >&2
    echo "Skipping pocket identification." >&2
else
    # Run pocket identification
    python3 "$POCKET_SCRIPT" \
        -i "$CENTROID_FILE" \
        -o "$OUTPUT_DIR" \
        -t "${CLUSTER_THRESHOLD:-5.0}"

    if [[ $? -eq 0 ]]; then
        echo ""
        echo "✓ Pocket identification complete"
    else
        echo "Error: Pocket identification failed" >&2
        echo "Continuing with available results..." >&2
    fi
fi


echo ""
echo "========================================"
echo "Analysis Complete!"
echo "========================================"
echo "Output files:"
echo "  aligned_cifs/           Aligned structures"
echo "  ligand_centers.csv      Ligand centroids"
echo "  ligand_deviations.csv   Distances from average"
if [[ -f "$OUTPUT_DIR/pocket_assignments.csv" ]]; then
    echo "  pocket_assignments.csv  Pocket IDs per ligand"
    echo "  cluster_statistics.csv  Per-pocket cluster statistics"
fi
echo "========================================"
