# boltz2-hpc-tools

Scripts for large-scale virtual screening using [Boltz-2](https://github.com/jwohlwend/boltz) on HPC systems with SLURM.

Developed by Greenstone Biosciences (https://greenstonebio.com).
For questions or collaboration: Jeremy Leitz jeremyleitz@greenstonebio.com

---

## Requirements

- [Boltz](https://github.com/jwohlwend/boltz) — install separately following their instructions
- Conda or Mamba
- SLURM workload manager
- Python 3.x with `pandas` and `rdkit` (installed into your Boltz conda environment)

---

## Installation

1. **Install Boltz** following the [official instructions](https://github.com/jwohlwend/boltz):
```bash
pip install boltz[cuda] -U
```

2. **Clone this repository:**
```bash
git clone https://github.com/Greenstone-Biosciences/boltz2-hpc-tools.git
cd boltz2-hpc-tools
```

3. **Install additional dependencies** into your Boltz environment:
```bash
conda activate boltz   # or whatever you named your environment
conda install -c conda-forge pandas rdkit
```

4. **Run setup.sh:**
```bash
[super_cool_user]$ ./setup.sh
Path to this repository [/home/folder/boltz2-hpc-tools]: /home/folder/boltz2-hpc-tools
Conda environment name [boltz]: my_boltz_env
Setup complete!
```

---

## Workflow Overview

The scripts are collectively called **Nutz and Boltz**. The general workflow is:

```
prep_batch_boltz  →  sbatch <slurm_file>  →  prep_analyze  →  sbatch <slurm_file>
```

Use `-h` on any script for detailed usage.

---

## Step 1 — Prepare and submit a screening job (`prep_batch_boltz`)

`prep_batch_boltz` fetches the protein FASTA from UniProt (or accepts a local file), splits
your SMILES file into chunks for SLURM array processing, and writes a ready-to-submit
`.slurm` file.

```bash
prep_batch_boltz -p <UniProtID_or_fasta> -s <smiles_file> [options]
```

### Required arguments

| Flag | Description |
|------|-------------|
| `-p` | UniProt ID (e.g. `P00519`) or path to a local FASTA file |
| `-s` | SMILES input file (see formats below) |

### Optional arguments

| Flag | Description |
|------|-------------|
| `-n` | Job name [default: UniProt ID or FASTA basename] |
| `-j` | Number of SLURM array chunks [default: 1] |
| `-o` | Output directory [default: current working directory] |
| `-m` | Path to a pre-computed MSA file for chain 1. If omitted, a pilot job runs automatically to generate the MSA before the main array job is written. |
| `--UniProtID2` | UniProt ID or FASTA path for a second protein chain (e.g. for a receptor/co-receptor complex) |
| `--use_msa2` | Pre-computed MSA for chain 2. If omitted and chain 2 is provided, the pilot job generates it. |
| `-b` | Max simultaneous SLURM array tasks (throttle). Requires `-j > 1`. |
| `-c` | SMILES column specifier (see below) |

### Input file formats

**CSV (`.csv`)** — recommended. Must have a header row. The SMILES column is
auto-detected from common names (`smiles`, `smi`, `canonical_smiles`, `structure`).
Use `-c <column_name>` to specify a non-standard column name.

```
smiles,name,ic50
CCO,ethanol,100
...
```

**Delimited (`.smi`, `.txt`, etc.)** — whitespace or tab-separated. SMILES is
assumed to be in column 1 and the identifier in column 2. Use `-c <column_number>`
(1-based) to specify a different column.

```
CCO ethanol
...
```

### MSA handling

`prep_batch_boltz` handles MSA generation automatically:

- **MSA provided (`-m`)**: path is baked directly into the SLURM script. No pilot job needed.
- **MSA not provided**: a single-compound pilot job is submitted immediately via `sbatch --wait`. The script blocks until it completes, locates the generated `uniref.a3m`, and bakes that path into the final SLURM script. You end up with one file and one `sbatch` command regardless.
- **Multi-chain**: each chain is evaluated independently. If chain 1 has an MSA but chain 2 does not, a single pilot job generates only the missing one.

### Examples

```bash
# Single chain, CSV input — MSA generated automatically:
prep_batch_boltz -p P00519 -s molecules.csv -n ABL1_screen -j 10 -o /data/output

# Non-standard SMILES column name:
prep_batch_boltz -p P00519 -s molecules.csv -n ABL1_screen -j 10 -c Structure

# Pre-computed MSA — pilot skipped:
prep_batch_boltz -p P00519 -s molecules.csv -n ABL1_screen -j 10 -m /data/ABL1_uniref.a3m

# Two-chain complex, no MSAs provided:
prep_batch_boltz -p P00519 --UniProtID2 Q9Y6K9 -s molecules.csv -j 10 -o /data/output

# Two-chain, one MSA provided:
prep_batch_boltz -p P00519 --UniProtID2 Q9Y6K9 -s molecules.csv -m chain1.a3m
```

### Output

```
<output_dir>/
    <UniProtID>.fasta
    Nutz/
        chunk_0000 ... chunk_NNNN   # SMILES chunks with header
        msa_pilot/                  # pilot job output (if MSA was generated)
        <name>_<stamp>.slurm        # ready-to-submit main array job
```

Submit the printed `sbatch` command to launch the screen.

> **Note:** The generated `.slurm` file can be edited before submission to add
> time limits (`#SBATCH --time=...`), memory constraints, or alternate partitions
> as needed for your system.

---

## Step 2 — Run predictions (`Nutz.sh`)

`Nutz.sh` is called automatically by the SLURM array job written by `prep_batch_boltz`.
It reads each chunk, validates SMILES with RDKit, generates Boltz-2 YAML input files,
and runs `boltz predict` on the chunk.

Compounds with unparseable SMILES are skipped and logged to
`Nutz/Nut_<chunk>/skipped_compounds.tsv` rather than crashing the job.

You should not normally need to call `Nutz.sh` directly.

---

## Step 3 — Prepare analysis (`prep_analyze`)

Once predictions are complete, run `prep_analyze` to generate an analysis SLURM script:

```bash
prep_analyze -i <input_dir> -o <output_dir> -J <job_name> -n <job_number> \
             -p <parallel_processes> -P <partition> -c <cpus> -m <memory> [-w]
```

| Flag | Description |
|------|-------------|
| `-i` | Input directory containing `Ligand_yamls/` |
| `-o` | Output directory |
| `-J` | SLURM job name |
| `-n` | Number of jobs |
| `-p` | GNU parallel `-j` value (must be ≤ `-c`) |
| `-c` | CPUs per node |
| `-m` | SLURM memory allocation |
| `-w` | Overwrite existing `combine.csv` |

---

## Step 4 — Analyze results (`analyze_batch_boltz`)

Collects all Boltz-2 outputs and concatenates them into a single `combine.csv`:

```bash
analyze_batch_boltz -i <dir_containing_Ligand_yamls> -o <output_dir> [-w]
```

---

## Inverse screening (`prep_inverse_screen`)

To dock a SMILES library against multiple protein targets in one run, use
`prep_inverse_screen`. It loops over a target list, runs a pilot MSA job per
target as needed, and generates a per-target SLURM pair with automatic
`--dependency` chaining.

```bash
prep_inverse_screen -t <targets_file_or_comma_list> -s <smiles_file> [options]
```

Targets can be provided as a file (one UniProt ID per line) or as a
comma-separated list directly on the command line:

```bash
prep_inverse_screen -t P00519,Q9Y6K9,P12345 -s molecules.csv -j 5 -o /data/output
```

---

## Updates

**2025-01** — Major overhaul of `prep_batch_boltz` (formerly `Deez.sh`) and `Nutz.sh`:
- CSV files are now parsed with `pandas`; columns are detected by header name rather than position, correctly handling files with commas in field values
- MSA generation is now fully automated via a pilot job — no need to pre-generate or manually provide MSA files
- Multi-chain support: provide multiple UniProt IDs / FASTA files and MSAs as comma-separated lists; missing MSAs are generated per-chain
- SMILES column can be specified with `-c` for non-standard column names or positions
- RDKit validation added: invalid SMILES are skipped and logged rather than crashing the job
- `--rm_header` flag removed — headers are now expected and handled automatically
- Added `prep_inverse_screen` for multi-target inverse virtual screening

**2025-12-13** — Added `prep_analyze` for aggregating Boltz-2 outputs into a combined CSV.

**2025-12-03** — `Nutz.sh` updated to use Boltz-2 batch mode (directory of YAML files).
Approximately 5x faster than previous per-ligand method. Added `-b` throttle flag.

---

## Acknowledgment

**Boltz Workflow Scripts for HPC**
Developed by Jeremy Leitz @ Greenstone Biosciences
https://github.com/Greenstone-Biosciences/boltz2-hpc-tools

If you use these scripts in your work, please acknowledge this repository:
```
https://github.com/Greenstone-Biosciences/boltz2-hpc-tools
```

If you use *Boltz-1/2* in your research, please cite:
```bibtex
@article{passaro2025boltz2,
  author  = {Passaro, Saro and Corso, Gabriele and Wohlwend, Jeremy and Reveiz, Mateo
             and Thaler, Stephan and Somnath, Vignesh Ram and Getz, Noah and Portnoi, Tally
             and Roy, Julien and Stark, Hannes and Kwabi-Addo, David and Beaini, Dominique
             and Jaakkola, Tommi and Barzilay, Regina},
  title   = {Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction},
  year    = {2025},
  doi     = {10.1101/2025.06.14.659707},
  journal = {bioRxiv}
}

@article{wohlwend2024boltz1,
  author  = {Wohlwend, Jeremy and Corso, Gabriele and Passaro, Saro and Getz, Noah
             and Reveiz, Mateo and Leidal, Ken and Swiderski, Wojtek and Atkinson, Liam
             and Portnoi, Tally and Chinn, Itamar and Silterra, Jacob and Jaakkola, Tommi
             and Barzilay, Regina},
  title   = {Boltz-1: Democratizing Biomolecular Interaction Modeling},
  year    = {2024},
  doi     = {10.1101/2024.11.19.624167},
  journal = {bioRxiv}
}
```
