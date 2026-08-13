# Pocket Analysis Pipeline — Claude Code Context

## Scope

This directory is the pocket-analysis pipeline for identifying small-molecule
binding sites on cytokine receptor subunits (IL11Ra, IL6Ra, GP130, IL11_IL11RA
complex) from Boltz2 co-folding predictions.

Free to edit, rename, restructure, split, or add files **within this directory**.
Do not modify anything in the parent repo outside `pocket_analysis/` (Jeremy
Leitz's Boltz2 prep/analysis scripts: prep_batch_boltz, analyze_batch_boltz.sh,
Nutz.sh, setup.sh, prep_analyze, Archive/) unless explicitly asked.

## Pipeline architecture

Boltz2 CIF outputs (HPC, /data/jleitz/Debarun/, other locations of boltz2 runs with .cifs)
|
pocket_pipeline.sh --save-seed <- fresh DBSCAN run, establishes pocket coords
|
pocket_anchors.json <- FROZEN: pocket centroids, radii, member ledger
reference.cif <- FROZEN: alignment coordinate frame
|
pocket_pipeline.sh --seed <- each subsequent library (ApexBio, etc.)
|
pocket_assignments.csv <- per-library compound -> pocket assignments
|
merge_pocket_runs.py <- combines all libraries into one master view
|
pocket_assignments_all.csv
pocket_summary.csv
|
query_affinity.py <- adds adj_pred scores, ranks per pocket
|
per-pocket ranked CSVs <- deliverable for Jeremy

## Scripts

- `pocket_pipeline.sh` — bash orchestrator. CIF discovery, calls the two Python
  scripts below. `--seed DIR` / `--save-seed DIR` toggle seeded vs fresh mode.
    `SCRIPT_DIR` resolves via `${BASH_SOURCE[0]}` — must stay in the same
      directory as the two Python scripts it calls.
      - `align_extract_ligands.py` — CA-atom superposition via gemmi, ligand
        centroid extraction. Accepts `--reference-cif` to lock alignment frame
          across runs (used in seeded mode).
          - `identify_pockets.py` — DBSCAN clustering (fresh mode) or bounding-sphere
            assignment against `pocket_anchors.json` (seeded mode). Star topology:
              anchors are frozen after first save, never mutated by later library runs.
              - `merge_pocket_runs.py` — stateless merge of multiple `pocket_assignments.csv`
                runs into one master CSV. Library label inferred from `cif_input_dir` in
                  each run's `run_summary.json`.
                  - `query_affinity.py` — queries Boltz2 affinity JSONs, computes
                    `adj_pred_0 = affinity_pred_value * affinity_probability_binary`, ranks
                      ascending (most negative = strongest predicted binding = rank 1).

## Conventions

- Conda env `boltz` must be active for gemmi/scikit-learn/numpy.
- Long runs (ApexBio-scale, 10k+ CIFs) go in tmux, not foreground.
- Don't run SLURM unless the job genuinely needs it — most of this pipeline
  runs fine interactively on a login/compute node.
  - Defer commits until functionality is confirmed working — don't commit
    mid-debug or mid-production-run.
    - Prefer directory names over CLI flags for library labeling where possible
      (matches how `merge_pocket_runs.py` infers labels).
      - CSV output should never include raw filesystem paths in a column that also
        gets sorted/parsed downstream — path commas break naive CSV tooling. Use
          Python's `csv` module for anything that writes or reads these files, not
            awk/sed, once a column could contain a comma.

## Known open items (as of this branch's creation)

- Fix 1: `reference_aligned.cif` appears as a phantom compound in
  `align_extract_ligands.py` output — reference structure gets its ligand
    centroid extracted like any other compound. Needs `ref_name` threaded
      through to `extract_ligand_centroids` to skip it.
      - Fix 4: strip `_aligned.cif` suffix from structure names in future runs
        (`extract_ligand_centroids` — change `f"{name}_aligned.cif"` to `name`).
          Do not apply retroactively to existing production CSVs.
          - IL11_IL11RA complex has low selfcheck rate at eps=5.0 (~31%) vs eps=7.0
            (~58.7%) — higher structural variability than monomers. eps=7.0 seed
              exists but hasn't been used for the full ApexBio assignment yet.
              - Pipeline README section still needed: pocket_pipeline.sh usage, --seed/
                --save-seed workflow, pocket_anchors.json structure, query_affinity.py
                  usage.

## Do NOT

- Do not touch files outside `pocket_analysis/` without being asked.
- Do not retroactively rename/reformat columns in already-delivered
  production CSVs under /data/cyan/ibrahim_pocket_analysis/production_*/
    or /data/Shared/ibrahim_pocket_analysis/production_*/ — those are Jeremy's
      reference copies.
      - Do not commit directly to `main` — this work happens on feature branches
        off `main`, PR'd back when validated.
