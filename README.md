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

4. Create configuration file:
```
bash
cp config.template config.sh
nano config.sh #Edit relevant paths to mach your system
```

## Use
These programs collecivetly and affectionately called "Nutz and Boltz" consist of a couple of wrapper scripts to run Boltz-2 using SLURM array jobs.
For each program you may do -h flag to see a detailed usage. 

The general workflow is:


## Acknowledgment

If you use these scripts in your work, please link to or acknowledge this repository:
```
https://github.com/Greenstone-Biosciences/boltz2-hpc-tools.git
```

**Boltz Workflow Scripts for HPC**  
Developed by Jeremy Leitz @ Greenstone Biosciences  
https://github.com/Greenstone-Biosciences/boltz2-hpc-tools

If you use Boltz in your research, please cite:
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
