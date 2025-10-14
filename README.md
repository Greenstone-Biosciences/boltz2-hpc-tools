# boltz2-hpc-tools
A collection of (maybe) useful scripts for screening using Boltz-2.

## About

Developed by Greenstone Biosciences ([https://greenstonebio.com]) 
to streamline Boltz workflows on HPC systems.

For questions or collaboration contact: Jeremy Leitz jeremyleitz@greenstonebio.com

## Requirements

- [Boltz](https://github.com/jwohlwend/boltz) - Install separately following their instructions
- Conda or Mamba
- Slurm workload manager (for HPC batch jobs)

## Installation
1. **Install Boltz** (if not already installed):
```
bash
pip install boltz[cuda] -U
# Best to follow Boltz isntallation instructions:
https://github.com/jwohlwend/boltz
```

2. Clone this repository:
```
bash
git clone https://github.com/Greenstone-Biosciences/boltz2-hpc-tools.git
cd boltz2-hpc-tools
```

3. Install additional dependencies into your Boltz environment:
```
bash
conda activate boltz # or whatever you named your Boltz environment
conda install -c conda-forge yq # YAML processor
conda install -c conda-forge jq # JSON processor
```

4. Run the setup.sh script:
```
bash
[super_cool_user]$ pwd
/home/folder_that_makes_sense/boltz2-hpc-tools
[super_cool_user]$ ./setup.sh
Path to this repository [/home/folder_that_makes_sense/boltz2-hpc-tools?]: /home/folder_that_makes_sense/boltz2-hpc-tools
Conda environment name [boltz]: my_boltz_env
Path to the affinity_template yaml [/home/folder_that_makes_sense/boltz2-hpc-tools]: /home/folder_that_makes_sense/boltz2-hpc-tools
Setup complete!
```


## Use
These programs collecivetly and affectionately called "Nutz and Boltz" consist of a couple of wrapper scripts to run Boltz-2 using SLURM array jobs.
For each program you may do -h flag to see a detailed usage. 

The general workflow is:
1) Fetch protein fasta using uniprot ID & split smiles file into n chunks for each slurm array job to process.
```
./prep_batch_boltz -p <uniprot_ID> -s <file_containing_smiles> [-n <Run_name>] [-j <Slurm_array_job_number>] [-o <path_to_output_directory>] [--rm_header] [-m <msa_file_path>]
```
The file containing SMILES strings may be comma or space separated, but the SMILES strings MUST be in the first column, followed by an Identifier in the second column.  Other columns are ignored. Optional flags: If there is a header in the SMILES file, include the flag --rm_header to ignore it. If you already have a msa file generated include the full path to that file after the -m flag. To run as only a single job, -j = 1 or omit the flag entirely. 

The output will be a .slurm submission file. This file should be modified as needed for your use case (i.e. add #SBATCH time=??? or select specific GPU/partitions to use if required by your system). Simply run using sbatch (the command itself will be echoed to the terminal).

2) Nutz.sh generates a yaml file for each line in the SMILES chunk file. Each yaml is then processed in boltz.  Note: boltz must be a command in currently in your path. 

Each SMILE string will generate a boltz_results_IDENTIFYIER folder within the Ligand_Yamls directory in Nutz. 

3) Analyze the entire Ligand_yamls directory and concatinate output into a file combine.csv in the output directory. Optionally overwrite an existing combine.csv file in the same directory using the -w flag.

```
./analyze_batch_boltz -i <Input Directory that contains 'Ligand_yamls' directory> -o <Output Directory> [-w/--overwrite]
```

Please feel free to contact us (me) with suggestions or changes. 



## Acknowledgment

**Boltz Workflow Scripts for HPC**
Developed by Jeremy Leitz @ Greenstone Biosciences
https://github.com/Greenstone-Biosciences/boltz2-hpc-tools
If you use these scripts in your work, please link to or acknowledge this repository:
```
https://github.com/Greenstone-Biosciences/boltz2-hpc-tools.git
```

If you use *Boltz-1/2* in your research, please cite:
```
@article{passaro2025boltz2,
  author = {Passaro, Saro and Corso, Gabriele and Wohlwend, Jeremy and Reveiz, Mateo and Thaler, Stephan and Somnath, Vignesh Ram and Getz, Noah and Portnoi, Tally and Roy, Julien and Stark, Hannes and Kwabi-Addo, David and Beaini, Dominique and Jaakkola, Tommi and Barzilay, Regina},
  title = {Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction},
  year = {2025},
  doi = {10.1101/2025.06.14.659707},
  journal = {bioRxiv}
}

@article{wohlwend2024boltz1,
  author = {Wohlwend, Jeremy and Corso, Gabriele and Passaro, Saro and Getz, Noah and Reveiz, Mateo and Leidal, Ken and Swiderski, Wojtek and Atkinson, Liam and Portnoi, Tally and Chinn, Itamar and Silterra, Jacob and Jaakkola, Tommi and Barzilay, Regina},
  title = {Boltz-1: Democratizing Biomolecular Interaction Modeling},
  year = {2024},
  doi = {10.1101/2024.11.19.624167},
  journal = {bioRxiv}
}
```
