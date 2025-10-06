#!/bin/bash

conda init
conda activate dock

#Step one: remove non-protein molecules
Receptor_raw="$1" #This should be the receptor pdb (eg 8dd0.pdb)
LigOrRes="$2" #How to identify the binding site, "L", "R", "M" or "Manual" determines the type of search that should be done to identify the binding site.
Query="$3" #If Ligand should be a three letter code. If residue include Chain ID,AA residue and Residue number in the format "D-MET-355"
BoxSize="$4" #Boxsize for the ligand file
Cut="$5" #Cut thre ceptor to desired Chains. Input can be either "full" or desired chains separated by commas (eg. B,C,D for chains B-D)
OutputDir="$6" #Optional output directory.  If blank then will be put into current directory. 

if [ "${Receptor_raw}" == "help" ]; then echo 'Welcome to PrepareVina.sh. This script will start from a pdb, convert the pdb to pdbqt using Autodock.py, and generate a Ligand.coord file that will contain docking coordinates and boxsize information to be used in Autodock Vina.  
	This scripts takes the following inputs:
	1) Receptor PDB, please include the complete path
	2) Method to identify Docking Coordinates (options are, "Manual", "L", "R", or "M"). 
	   L = Use existing Ligand positions. 
	   R = Use residue location.
	   M = Use a file containing multible residue locations.
	3) The search query.  
	   If L was selected above, provide the 3-letter ligand code in the pdb
	   If R was selected above, provide the 3-letter AA-ChainID-ResidueNumber (e.g. MET-D-355 for methionine 355 on chain D of the pdb)
	   If M was selected above, provide the path of a file that contains several residues in the same format as above (AA-ChainID-ResidueNumber)
	   If Manual was selected above, then provide the coordinates separated by some delimter (x,y,z e.g. 10,0,-15)
	4) BoxSize that will be used in the docking campaign, currently only cubes are generated, manual modification is required if rectangles are desired.
	5) Chains to be used in the structure.  If the whole pdb may be used then input "full". Otherwise denote DESIRED chains separated by commas (B,C,E for chains B C and E)
	6) Optional output directory.  If not provided a folder "VinaPrep" will be created in the current working directory.
	
	Good luck!';
	exit 1; 
fi
	   

#Checks if folder directory exists
if [ -z "$OutputDir" ]; then OutputDir=$(pwd)"/VinaPrep"; 
mkdir $OutputDir
else
	if [[ ! -d $OutputDir ]]; then mkdir $OutputDir; fi
	OutputDir=$OutputDir/VinaPrep
	mkdir $OutputDir
	echo "how bout now" $OutputDir
fi

	#This section converts the pdb to pdbqt and protonates relevant hydrogens

AD=/data/Shared/Docking_Scripts/Autodock.py
DIR=$(dirname "${Receptor_raw}")
BaseName=$(basename "${Receptor_raw}" .pdb)
Base_SM=$BaseName"_noSM.pdb"
#This section checks if the receptor should be cut or not
if [ "${Cut}" = "full" ]; then
	grep 'ATOM' $Receptor_raw > $DIR/$Base_SM
	wait
else
	echo 'cutting '$Receptor_raw' to chains '${Cut//,/ }
	grep 'ATOM' $Receptor_raw | awk -v env_var="${Cut//,/ }" 'BEGIN{split(env_var,t); for (i in t) vals[t[i]]} ($5 in vals)' > $DIR/$Base_SM
fi

python3 $AD -r $DIR/$Base_SM
wait
mv ./receptor.pdbqt $OutputDir/$BaseName".pdbqt"

#This part finds the box location based on existing ligands in the pdb
if [ "${LigOrRes}" = "" ]; then 
	echo "Ligand or Residue search not defined!"
	elif [ "${LigOrRes}" = "L" ]; then
		echo 'Searching for Ligand in '$Receptor_raw
		grep 'HETATM' $Receptor_raw | grep "${Query}" | awk '{ print $3, $4 }' | uniq > $OutputDir/tmp_Ligand.coords
		LigNum=$(cat $OutputDir/tmp_Ligand.coords | wc -l)
		#Some informational messages
			if [ "${LigNum}" = "0" ]; then 
				echo "Didn't find that ligand";
			elif [ "${LigNum}" -eq 1 ]; then
				echo "Found 1 ligand."
			elif [ "${LigNum}" -gt 1 ]; then
				echo "Found "$LigNum "ligands.  Finding All coordinates..."
			fi
		#Find coordinates for all of the ligands
		echo "Ligand code Lig ID XPos YPos ZPos Boxsize X Y Z" > $OutputDir/$BaseName"_Ligand.coord"
		for length in $(seq 1 $LigNum); 
		do 
			Lig_ID=$(sed "${length}q;d" $OutputDir/tmp_Ligand.coords | awk '{ print $2 }');
			len=$(grep "${Lig_ID}" $Receptor_raw | wc -l)
			sum=$(grep "$Lig_ID" $Receptor_raw | grep HETATM | awk '{ print $5 }' | tail -n "$len" | paste -sd+ - | bc);
			Xpos=$(echo "$sum / $len" | bc -l);

			sum=$(grep "$Lig_ID" $Receptor_raw | grep HETATM | awk '{ print $6 }' | tail -n "$len" | paste -sd+ - | bc); 
			Ypos=$(echo "$sum / $len" | bc -l);

			sum=$(grep "$Lig_ID" $Receptor_raw | grep HETATM | awk '{ print $7 }' | tail -n "$len" | paste -sd+ - | bc); 
			Zpos=$(echo "$sum / $len" | bc -l);

			LigInfo=$(sed "${length}q;d" $OutputDir/tmp_Ligand.coords)
			echo $LigInfo $Xpos $Ypos $Zpos $BoxSize $BoxSize $BoxSize >> $OutputDir/$BaseName"_Ligand.coord"
		done
		echo "Found all coordinates"
		# Find coordinates based on single Residue position with Query format of AminoAcid-ChainID-ResidueNumber (eg. MET-D-355)
	elif [ "${LigOrRes}" = "R" ]; then
		echo 'Searching for Docking coordinates using Receptor '$Query' residues'
		Xpos=$(cat $Receptor_raw | tr -s ' ' | grep "${Query//-/ }" | tail -n1 | awk '{ print $7 }')
		Ypos=$(cat $Receptor_raw | tr -s ' ' | grep "${Query//-/ }" | tail -n1 | awk '{ print $8 }')
		Zpos=$(cat $Receptor_raw | tr -s ' ' | grep "${Query//-/ }" | tail -n1 | awk '{ print $9 }')
		echo "ResInfo Xpos Ypos Zpos BoxsizeX Y Z" > $OutputDir/$BaseName"_Ligand.coord"
		echo $Query $Xpos $Ypos $Zpos $BoxSize $BoxSize $BoxSize>> $OutputDir/$BaseName"_Ligand.coord"

#Fine coordinates based on several residue positions defined in a file
	elif [ "${LigOrRes}" = "M" ]; then
		echo "Queryfile XPos Ypos Zpos Boxsize X Y Z" > $OutputDir/$BaseName"_Ligand.coord"
		echo 'Creating geometric average of coordinates provided in '$Query' file.'
		declare -a sumX=();
		declare -a interX=();
		declare -a sumY=();
		declare -a interY=();
		declare -a sumZ=();
		declare -a interZ=();
		for length in $(seq 1 $(cat ${Query} | wc -l)); 
			do	
			Xsearch=$(sed "${length}q;d" $Query);
		       echo 'searching '$Xsearch	
			sumX+=( "$(cat $Receptor_raw | tr -s ' ' | grep "${Xsearch//-/ }" | tail -n1 | awk '{ print $7 }')" );
			sumY+=( "$(cat $Receptor_raw | tr -s ' ' | grep "${Xsearch//-/ }" | tail -n1 | awk '{ print $8 }')" );
			sumZ+=( "$(cat $Receptor_raw | tr -s ' ' | grep "${Xsearch//-/ }" | tail -n1 | awk '{ print $9 }')" );
			done
			interX=(${sumX[@]%.*});
			interY=(${sumY[@]%.*});
			interZ=(${sumZ[@]%.*});
		echo ${interX[@]}	
			totX=0;
			totY=0;
			totZ=0;
			len=${#interX[@]}
			for ele in ${interX[@]}; do
				let totX=totX+ele;
			done
			let AvgX=totX/len
			for ele in ${interY[@]}; do
				let totY=totY+ele;
			done
			let AvgY=totY/len
			for ele in ${interZ[@]}; do
				let totZ=totZ+ele;
			done
			let AvgZ=totZ/len
		echo $Query $AvgX $AvgY $AvgZ $BoxSize $BoxSize $BoxSize>> $OutputDir/$BaseName"_Ligand.coord"

	#Manual entry
        elif [ "${LigOrRes}" = "Manual" ]; then	
		echo 'Searching for Docking coordinates using manual coordinates'
		#delim=$(printf '%s\n' "${Query//[[:digit:]]/}" | fold -w1 | sort -u) fucks up negatives
		pos=$(echo $Query | tr "," "\n" )
                echo "ResInfo Xpos Ypos Zpos BoxsizeX Y Z" > $OutputDir/$BaseName"_Ligand.coord"
                echo "Manual_Entry" $pos $BoxSize $BoxSize $BoxSize>> $OutputDir/$BaseName"_Ligand.coord"
fi




