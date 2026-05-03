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
SAVE_ALIGNED=false #Default: do not save aligned CIFs
SEED_DIR=""        # --seed: load reference.cif + pocket_anchors.json from here
SAVE_SEED_DIR=""   # --save-seed: save reference.cif + pocket_anchors.json here

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
    --save-aligned            Save aligned CIF files (default: false)
    --dry-run                 Show what files would be processed without running analysis
    --seed DIR                Load seed dir (reference.cif + pocket_anchors.json) for cross-run pocket assignment
    --save-seed DIR           After run, save reference.cif and pocket_anchors.json to this directory for use in future seeded runs
    -v, --verbose             Verbose output
    -h, --help                Show this help message

Requirements:
    - gemmi Python package
    - bc calculator
    - numpy Python package
    - scikit-learn Python package
    - (boltz) conda environment
    - align_and_extract_ligands.py (in same directory)
    - identify_pockets.py (in same directory)

Input Requirements:
    - Directory must contain .cif files

Output:
    - aligned_cifs/           Aligned .cif files (only if --save-aligned)
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
        --save-aligned)
            SAVE_ALIGNED=true
            shift
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
        --seed)
            SEED_DIR="$2"
            shift 2
            ;;
        --save-seed)
            SAVE_SEED_DIR="$2"
            shift 2
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

# Validate seed directory if provided
if [[ -n "$SEED_DIR" ]]; then
    if [[ ! -d "$SEED_DIR" ]]; then
         echo "Error: Seed directory '$SEED_DIR' not found" >&2
         exit 1
    fi
    if [[ ! -f "$SEED_DIR/reference.cif" ]]; then
        echo "Error: No reference.cif found in seed directory '$SEED_DIR'" >&2
        exit 1
    fi
    if [[ ! -f "$SEED_DIR/pocket_anchors.json" ]]; then
        echo "Error: No pocket_anchors.json found in seed directory '$SEED_DIR'" >&2
        exit 1
    fi
fi

# Validate save-seed dir is different from seed dir
if [[ -n "$SEED_DIR" && -n "$SAVE_SEED_DIR" && "$SEED_DIR" == "$SAVE_SEED_DIR" ]]; then
    echo "Error: --seed and --save-seed must point to different directories" >&2
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
#if ! command -v USalign &> /dev/null; then
#	echo "Error: US-align not found in PATH" >&2
#	echo "Please install US-align or ensure it's in your conda environment" >&2
#	exit 1
#fi

# Log function
log() {
	if [[ "$VERBOSE" == true ]]; then
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
	fi
}



# Print configuration
echo "============================================"
echo "Boltz2 Docking Analysis"
echo "============================================"
echo "Input directory:      $INPUT_DIR"
echo "Output directory:     $OUTPUT_DIR"
echo "CIF files found:      $CIF_COUNT"
echo "Cluster threashold:   $CLUSTER_THRESHOLD Å"
echo "Save aligned:         $SAVE_ALIGNED"
if [[ -n "$SEED_DIR" ]]; then
    echo "Seed directory:       $SEED_DIR"
fi
if [[ -n "$SAVE_SEED_DIR" ]]; then
    echo "Save seed to:         $SAVE_SEED_DIR"
fi
echo "============================================"
echo ""

log "Starting analysis pipeline of .cifs..."
echo ""


###############################################################################
## STEP 1-2: ALIGN STRUCTURES AND EXTRACT LIGAND CENTROIDS
###############################################################################

echo "Steps 1-2: Aligning structures and extracting ligand centroids"
echo ""

# Check if align_extract_ligands.py exists
ALIGN_SCRIPT="$SCRIPT_DIR/align_extract_ligands.py"
if [[ ! -f "$ALIGN_SCRIPT" ]]; then
    echo "Error: $ALIGN_SCRIPT not found" >&2
    echo "Please ensure align_extract_ligands.py is in the same directory as this script" >&2
    exit 1
fi

# Build commands with optional flags
ALIGN_CMD=(python3 "$ALIGN_SCRIPT" -i "$INPUT_DIR" -o "$OUTPUT_DIR")

if [[ "$SAVE_ALIGNED" == true ]]; then
    ALIGN_CMD+=(--save-aligned)
fi

if [[ "$VERBOSE" == true ]]; then
    ALIGN_CMD+=(-v)
fi

# Pass seed reference CIF to alignment script if in seeded mode
if [[ -n "$SEED_DIR" ]]; then
    ALIGN_CMD+=(--reference-cif "$SEED_DIR/reference.cif")
fi

# Run alignment and extraction
"${ALIGN_CMD[@]}"

if [[ $? -ne 0 ]]; then
    echo "Error: Alignment and extraction failed" >&2
    exit 1
fi

# Verify output file exists
CENTROID_FILE="$OUTPUT_DIR/ligand_centers.csv"
if [[ ! -f "$CENTROID_FILE" ]]; then
    echo "Error: ligand_centers.csv not created"
    exit 1
fi

echo ""


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
    # Build pocket identification command
    POCKET_CMD=(python3 "$POCKET_SCRIPT"
        -i "$CENTROID_FILE"
        -o "$OUTPUT_DIR"
        -t "${CLUSTER_THRESHOLD:-5.0}")

    # Seeded mode: load existing anchors
    if [[ -n "$SEED_DIR" ]]; then
        POCKET_CMD+=("--load-anchors" "$SEED_DIR/pocket_anchors.json")
    fi

    # Save anchors if --save-seed specified
    if [[ -n "$SAVE_SEED_DIR" ]]; then
        mkdir -p "$SAVE_SEED_DIR"
        POCKET_CMD+=("--save-anchors" "$SAVE_SEED_DIR/pocket_anchors.json")
    fi

    "${POCKET_CMD[@]}"

    if [[ $? -eq 0 ]]; then
        echo ""
        echo "✓ Pocket identification complete"
    else
        echo "Error: Pocket identification failed" >&2
        echo "Continuing with available results..." >&2
    fi
fi

# Save seed: copy reference CIF after successful run
if [[ -n "$SAVE_SEED_DIR" ]]; then
# Determine which CIF was used as reference (first alphabetically)
    if [[ -n "$SEED_DIR" ]]; then
        REF_CIF="$SEED_DIR/reference.cif"
    else
        REF_CIF="${CIF_FILES[0]}"
    fi
    cp "$REF_CIF" "$SAVE_SEED_DIR/reference.cif"
    echo "✓ Seed saved to $SAVE_SEED_DIR/"
    echo "  reference.cif   — alignment reference structure"
    echo "  pocket_anchors.json — pocket centroids and member ledger"
fi


echo ""
echo "========================================"
echo "Analysis Complete!"
echo "========================================"
echo "Output files:"
if [[ "$SAVE_ALIGNED" == true ]]; then
    echo "  aligned_cifs/           Aligned structures"
fi
echo "  ligand_centers.csv      Ligand centroids"
echo "  ligand_deviations.csv   Distances from average"
if [[ -f "$OUTPUT_DIR/pocket_assignments.csv" ]]; then
    echo "  pocket_assignments.csv  Pocket IDs per ligand"
    echo "  cluster_statistics.csv  Per-pocket cluster statistics"
fi
if [[ -n "$SAVE_SEED_DIR" ]]; then
    echo "  $SAVE_SEED_DIR/reference.cif"
    echo "  $SAVE_SEED_DIR/pocket_anchors.json"
fi
echo "========================================"
