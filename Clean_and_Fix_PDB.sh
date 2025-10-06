#!/bin/bash

#This script will use pdbfixer to fix a protein structure and then (optionally) meeko to convert it to pdbqt.  
set -e #Exit on error

#Default values
convert_flag=false
pH=7.4
chain="full"
box_size=15

#Parse arguments
while [[ $# -gt 0 ]]; do
	key="$1"
	case $key in
		-i|--input)
			input_pdb="$2"
			shift
			shift
			;;
		-o|--output)
			output_dir="$2"
			shift
			shift
			;;
		-m|--manual)
			manual_coords="$2"
			shift
			shift
			;;
		-b|--box_size)
			box_size="$2"
			shift
			shift
			;;
		--convert)
			convert_flag=true
			shift
			;;
		--chain)
			chain="$2"
			shift
			shift
			;;
		--pH)
			pH=$2
			shift
			shift
			;;
		--alphafold)
			AF_flag=true
			shift
			;;
		--find_lig)
			find_lig="$2"
			shift
			shift
			;;
		-h|--help)
			echo "Usage: $0 -i input_pdb -o output_dir -b box_size [--convert] [--chain <A|B|ect|full>] [--ph <ph_value>] [--alphafold] [--find_lig <ligand ID>]"
			exit 0
			;;
		*)
			echo "Unknown option: $1"
			echo "Usage: $0 -i input_pdb -o output_dir -b 15 [--convert] [--ph <ph_value>]"
			exit 1
			;;
	esac
done
#Check mandatory args
if [[ -z "$input_pdb" || -z "$output_dir" ]]; then
	echo "Error: input and output must be specified."
	echo "Usage: $0 -i input_pdb -o output_dir [--convert]"
	exit 1
fi

# Source the conda setup to make the `conda` command available.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate pdbfixer310

# Prepare paths
if [[ ! -d "$output_dir" ]]; then mkdir -p "$output_dir"; fi
base_name=$(basename "$input_pdb" .pdb)
	
## Manual coord generation in progress 20250928
## need to add box functionality
#Manual entry
if [[ -n "$manual_coords" ]]; then
	ligand_file="${output_dir}/VinaPrep/${base_name}_Ligand.coord" #Places in VinaPrep from start

	echo "Creating ligand coordinate file using manual coordinates"

	pos=$(echo $manual_coords)

        echo "Ligand,Xpos,Ypos,Zpos,BoxX,BoxY,BoxZ" > "$ligand_file"
        echo "Manual_Entry,${pos},${box_size},${box_size},${box_size}" >> "$ligand_file"
		        
        echo "Ligand coordinate file created: $ligand_file"
fi

# box_x, box_y, box_z
# if = 1 then do default
# elsif 3 then put all boxes valuuuuuuuuues in coord
# elsif none of those nums, If !== 1 || 3, throw errortoo few tooo maaaaaaany or       put 1 or 3
# output ligand coord to .csv format, hearder and values

#Step 1: Run structure fixing with MeekoCleanPDB
	python3 ~/Software/VerdeLeitz_scripts/MeekoCleanPDB_v9.py "$input_pdb" "$output_dir" --ph "$pH" ${AF_flag:+--alphafold}
	if [[ "$AF_flag" == true ]]; then echo "Using Alphfold model of protein PDB ID: ${base_name}"; fi 
	fixed_pdb="$output_dir/VinaPrep/${base_name}_fixed.pdb"
	raw_pdb="$output_dir/VinaPrep/${base_name}.pdb"
	if [[ -f "${base_name}.pdb" && ! -f "$raw_pdb" ]]; then
	    mv "${base_name}.pdb" "$raw_pdb"
	fi

##Step 1.5: If chain(s) are selected, keep only selected chains.	
if [[ "$chain" != "full" ]]; then
	trimmed_pdb="$output_dir/VinaPrep/${base_name}_fixed_chains.pdb"
	awk -v chains="$chain" '
		BEGIN {
			n = split(chains, chain_arr, ",")
			for (i = 1; i <= n; i++) chain_map[chain_arr[i]]
		}
		/^(ATOM  |HETATM)/ {
			c = substr($0, 22, 1)
			if (c in chain_map) print
		}
		/^TER/ || /^END/ {
			print
		}
	' "$fixed_pdb" > "$trimmed_pdb"
	mv "$trimmed_pdb" "$fixed_pdb"
	echo "✅ Retained only chains: $chain in cleaned PDB."
fi


#Step 2: Optionally remove HETATM and extract ligand centroid
if [[ -n "$find_lig" ]]; then
    lig_coord_file="$output_dir/VinaPrep/${base_name}_Ligand.coord"
    echo "Post-processing: writing cleaned PDB to $preprocessed_pdb"

# Post-Process: remove HETATM / extract ligand centroid
awk -v ligand="$find_lig" -v chain_sel="$chain" -v coord_out="$lig_coord_file" -v box="$box_size" '
    BEGIN { count=0 }
    /^HETATM/ {
        resname = substr($0, 18, 3)
        chain_id = substr($0, 22, 1)
        x = substr($0, 31, 8) + 0
        y = substr($0, 39, 8) + 0
        z = substr($0, 47, 8) + 0

        if (resname == ligand && (chain_sel == "full" || chain_id == chain_sel)) {
            lx[count] = x
            ly[count] = y
            lz[count] = z
            count++
        }
    }
    END {
        if (count > 0) {
            print "Ligand X Y Z" > coord_out
            sumx = sumy = sumz = 0
            for (i = 0; i < count; i++) {
                sumx += lx[i]; sumy += ly[i]; sumz += lz[i]
            }
            printf "%s %.2f %.2f %.2f %.2f %.2f %.2f \n", ligand, sumx/count, sumy/count, sumz/count, box, box, box > coord_out
            printf "✅ Ligand centroid written to: %s\n", coord_out
        } else {
            printf "⚠️ Ligand %s not found in HETATM records", ligand
            if (chain_sel != "full") printf " for chain %s", chain_sel
            print "."
        }
    }
' "$raw_pdb"

fi

#Step 3: Optionally convert to PDBQT
if [[ "$convert_flag" == true ]]; then
	
	#Check if mk_prepare_receptor.py is available
	if ! command -v mk_prepare_receptor.py &> /dev/null; then
		echo "Error: mk_prepare_receptor.py not found in current environment."
		echo "Make sure Meeko is installed and then environment is active"
		exit 1
	fi
	
	output_pdbqt="${fixed_pdb%.pdb}"
	echo "Converting $fixed_pdb to $output_pdbqt using mk_prepare_receptor.py"
	mk_prepare_receptor.py -i "$fixed_pdb" -o "$output_pdbqt" -p --write_json --allow_bad_res
fi


echo "Finished: Cleaned PDB ${convert_flag:+ and PDBQT }created!"
