#!/usr/bin/env python3
# Written by Chris Yan
"""
Merge pocket assignments from multiple library runs into a unified cross-library view.

PIPELINE POSITION: Reporting layer only. This script does not re-run clustering,
modify the anchor file, or change any pocket assignments. It reads existing output
files produced by valid_boltz.sh and combines them.

CONTEXT: The pocket coordinate system is established once by the CHEMBL seed run
(valid_boltz.sh --save-seed) and frozen. Every subsequent library is assigned against
that frozen system (valid_boltz.sh --seed), producing a pocket_assignments.csv per
library. This script combines those per-library CSVs into a single master view and
recalculates per-pocket statistics across all libraries.

Usage:
    python3 merge_pocket_runs.py \
        --runs RUN_DIR [RUN_DIR ...] \
        --output OUTPUT_DIR \
        [--min-pocket-size INT]

    # Example: merge CHEMBL and ApexBio runs
    python3 merge_pocket_runs.py \
        --runs \
            /data/cyan/ibrahim_pocket_analysis/runs/IL11Ra_CHEMBL \
            /data/cyan/ibrahim_pocket_analysis/runs/IL11Ra_ApexBio \
        --output /data/cyan/ibrahim_pocket_analysis/merged/IL11Ra

Outputs:
    pocket_assignments_all.csv   One row per compound across all libraries
    pocket_summary.csv           Per-pocket member counts broken down by library
    merge_summary.json           Audit record of which runs were merged and when

Planned:
    --deep-stats                 Load ligand_centers.csv from each run to recalculate
                                 pocket radius and spread from all member coordinates
                                 (currently counts only)
"""

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from datetime import datetime


# ---------------------------------------------------------------------------
# Schema definitions — single source of truth for expected file structures.
# Update these if upstream scripts change their output columns.
# ---------------------------------------------------------------------------

# Required columns in pocket_assignments.csv (produced by identify_pockets.py)
ASSIGNMENTS_REQUIRED_COLS = {"Structure", "Pocket", "Cluster_Size"}

# Fields we read from run_summary.json (produced by identify_pockets.py via valid_boltz.sh)
SUMMARY_REQUIRED_FIELDS = {"output_dir"}
SUMMARY_OPTIONAL_FIELDS = {"cif_input_dir", "ligand_centers_path", "seed_dir",
                            "mode", "timestamp", "n_input", "n_assigned_existing",
                            "n_assigned_new", "params", "n_overlapping_sphere_pairs"}


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge pocket assignments from multiple library runs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--runs",
        nargs="+",
        required=True,
        metavar="RUN_DIR",
        help="One or more run output directories (each must contain "
             "pocket_assignments.csv; run_summary.json strongly recommended)"
    )
    parser.add_argument(
        "--output",
        required=True,
        metavar="OUTPUT_DIR",
        help="Output directory for merged results (created if absent)"
    )
    parser.add_argument(
        "--min-pocket-size",
        type=int,
        default=1,
        metavar="N",
        help="Exclude pockets with fewer than N total members across all libraries "
             "(default: 1 — include all pockets including singletons)"
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# File loading and validation
# ---------------------------------------------------------------------------

def load_summary(run_path):
    """Load and validate run_summary.json from a run directory.

    run_summary.json is written by identify_pockets.py (Step 4 of the pipeline)
    and records the original CIF source directory, clustering parameters, and
    assignment counts. We use it here to derive the library label and to validate
    that all runs used the same seed anchor before merging.

    Returns a summary dict with missing optional fields defaulted to None.
    Returns a minimal stub if the file is absent (backward compat with older runs).
    """
    summary_file = run_path / "run_summary.json"

    if not summary_file.exists():
        print(f"  ⚠ Warning: run_summary.json not found in {run_path} — "
              f"library label will fall back to directory name", file=sys.stderr)
        return {"output_dir": str(run_path), "cif_input_dir": None}

    with open(summary_file) as f:
        try:
            summary = json.load(f)
        except json.JSONDecodeError as e:
            print(f"Error: run_summary.json in {run_path} is malformed: {e}", file=sys.stderr)
            sys.exit(1)

    # Validate required fields exist
    missing = SUMMARY_REQUIRED_FIELDS - set(summary.keys())
    if missing:
        print(f"Error: run_summary.json in {run_path} missing required fields: {missing}",
              file=sys.stderr)
        sys.exit(1)

    # Fill optional fields with None so callers can use .get() safely
    for field in SUMMARY_OPTIONAL_FIELDS:
        summary.setdefault(field, None)

    return summary


def load_assignments(run_path, library):
    """Load and validate pocket_assignments.csv from a run directory.

    pocket_assignments.csv is produced by identify_pockets.py (Step 4).
    Each row records one compound and its assigned pocket ID — either from
    fresh DBSCAN (seed run) or bounding sphere nearest-centroid (seeded run).

    Returns a list of dicts with Library and Run_Dir fields added.
    Exits with a clear error if the file is absent or schema has changed.
    """
    assignments_file = run_path / "pocket_assignments.csv"

    if not assignments_file.exists():
        print(f"Error: pocket_assignments.csv not found in {run_path}", file=sys.stderr)
        sys.exit(1)

    rows = []
    with open(assignments_file, newline="") as f:
        reader = csv.DictReader(f)

        # Validate columns against schema — catches upstream output changes early
        actual_cols = set(reader.fieldnames or [])
        missing_cols = ASSIGNMENTS_REQUIRED_COLS - actual_cols
        if missing_cols:
            print(f"Error: pocket_assignments.csv in {run_path} missing columns: {missing_cols}",
                  file=sys.stderr)
            print(f"  Found:    {sorted(actual_cols)}", file=sys.stderr)
            print(f"  Expected: {sorted(ASSIGNMENTS_REQUIRED_COLS)}", file=sys.stderr)
            sys.exit(1)

        for row in reader:
            try:
                rows.append({
                    "Structure": row["Structure"],
                    "Pocket": int(row["Pocket"]),
                    "Library": library,
                    # Resolved absolute path — preserved for traceability
                    # even if the directory is later moved or renamed
                    "Run_Dir": str(run_path)
                })
            except (ValueError, KeyError) as e:
                print(f"Error: malformed row in {assignments_file}: {row} — {e}",
                      file=sys.stderr)
                sys.exit(1)

    return rows


def infer_library_name(summary, run_path):
    """Derive a human-readable library label from available metadata.

    Resolution order (most to least semantically stable):

    1. cif_input_dir basename — the original Boltz2 source directory, recorded
       by valid_boltz.sh via --cif-input-dir. Most stable because it names the
       actual library directory set by Jeremy's upstream pipeline.
       e.g. /data/.../IL11Ra_monomer/ApexBio  →  'ApexBio'

    2. output_dir basename — the run output directory. Stable by naming convention
       but depends on how the user named the output dir.
       e.g. /data/.../IL11Ra_ApexBio_r1  →  'IL11Ra_ApexBio_r1'

    3. run_path basename — the directory passed to --runs. Last resort.

    Note: older runs recorded the ligand_centers.csv path in cif_input_dir rather
    than the CIF directory (predates the --cif-input-dir fix). We detect this by
    checking for a .csv suffix and skip to the next fallback.
    """
    cif_dir = summary.get("cif_input_dir")
    if cif_dir and Path(cif_dir).suffix != ".csv":
        return Path(cif_dir).name

    output_dir = summary.get("output_dir")
    if output_dir:
        return Path(output_dir).name

    return run_path.name


def load_run(run_dir):
    """Load all data from a single run directory.

    Returns:
        rows:    list of assignment dicts (Structure, Pocket, Library, Run_Dir)
        summary: parsed run_summary.json dict
        library: inferred library name string
    """
    run_path = Path(run_dir).resolve()

    if not run_path.exists():
        print(f"Error: run directory not found: {run_dir}", file=sys.stderr)
        sys.exit(1)

    summary = load_summary(run_path)
    library = infer_library_name(summary, run_path)
    rows = load_assignments(run_path, library)

    print(f"  {len(rows):>6} compounds   library='{library}'")
    print(f"           source: {summary.get('cif_input_dir') or 'unknown (pre-traceability run)'}")

    return rows, summary, library


# ---------------------------------------------------------------------------
# Validation across runs
# ---------------------------------------------------------------------------

def validate_seed_consistency(summaries, run_dirs):
    """Verify all seeded runs used the same pocket_anchors.json.

    If two runs were seeded against different anchor files their pocket IDs
    refer to different coordinate systems — merging them would be meaningless.
    Fresh (seed generation) runs have seed_dir=None and are exempt.

    Exits with a clear error listing the conflicting seed paths if inconsistency
    is detected. This must run before any output is written.
    """
    seed_dirs = {
        run_dir: summary["seed_dir"]
        for run_dir, summary in zip(run_dirs, summaries)
        if summary.get("seed_dir") is not None
    }

    unique_seeds = set(seed_dirs.values())
    if len(unique_seeds) > 1:
        print("\n⚠ ERROR: Runs were seeded against different anchor files. "
              "Pocket IDs are not comparable — aborting merge.", file=sys.stderr)
        for run_dir, seed_dir in seed_dirs.items():
            print(f"  {run_dir}\n    seed: {seed_dir}", file=sys.stderr)
        sys.exit(1)


def check_duplicate_structures(all_rows):
    """Warn if the same Structure name appears in more than one library.

    Does not exit — duplicates may be legitimate (same compound in both CHEMBL
    and ApexBio) but should be surfaced so they can be inspected if the merged
    counts look unexpectedly high.
    """
    seen = defaultdict(list)
    for row in all_rows:
        seen[row["Structure"]].append(row["Library"])

    dupes = {s: libs for s, libs in seen.items() if len(libs) > 1}
    if dupes:
        print(f"\n  ⚠ {len(dupes)} structure(s) appear in multiple libraries "
              f"(may be expected for compounds in multiple screens):")
        for struct, libs in list(dupes.items())[:5]:
            print(f"    {struct}: {libs}")
        if len(dupes) > 5:
            print(f"    ... and {len(dupes) - 5} more")
    else:
        print("  ✓ No duplicate structures across libraries")


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate_by_pocket(all_rows):
    """Group compound assignments by pocket ID and library.

    Core aggregation step. Produces a nested count:
        pockets[pocket_id][library] = n_compounds

    Pocket IDs are stable across runs because all seeded runs reference the
    same frozen CHEMBL anchor file (validated by validate_seed_consistency).

    Returns:
        pockets:       defaultdict {pocket_id: {library: count}}
        all_libraries: set of all library names present
    """
    pockets = defaultdict(lambda: defaultdict(int))
    all_libraries = set()

    for row in all_rows:
        pockets[row["Pocket"]][row["Library"]] += 1
        all_libraries.add(row["Library"])

    return pockets, all_libraries


# ---------------------------------------------------------------------------
# Deep stats stub
# Coordinate-based radius and spread recalculation — requires ligand_centers.csv.
# When implemented, reads per-compound centroid coordinates from each run and
# recalculates true geometric spread and max radius of each pocket across all
# libraries combined (as Jeremy requested for size updates).
# The seed centroid coordinates are NEVER updated — only descriptive stats change.
# ---------------------------------------------------------------------------

def calculate_coordinate_stats(pockets, run_dirs):
    """[PLANNED -- deep stats] Recalculate pocket radius and spread from coordinates.

    Currently a no-op returning None (counts-only mode is active).

    When implemented:
      1. Load ligand_centers.csv from each run_dir
         (columns: Structure, X, Y, Z — produced by align_extract_ligands.py)
      2. Join compound centroids to pocket assignments by Structure name
      3. Load seed pocket_anchors.json to get frozen centroid coordinates
      4. For each pocket, compute across all libraries combined:
         - spread: mean distance of all member centroids from seed centroid
         - radius: max distance of any member from seed centroid
      5. Return dict: {pocket_id: {'spread': float, 'radius': float}}

    To activate: replace the call in main() with this function's implementation
    and add --deep-stats flag to parse_args().
    """
    return None


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_assignments_all(all_rows, output_dir):
    """Write the master compound-level file combining all libraries.

    Concatenation of every run's pocket_assignments.csv with Library and
    Run_Dir columns added. One row per compound across all libraries.

    This is the file Jeremy can use to pull filenames for a given pocket:
        awk -F',' '$2==1 {print $1, $3}' pocket_assignments_all.csv
    """
    output_path = output_dir / "pocket_assignments_all.csv"
    fieldnames = ["Structure", "Pocket", "Library", "Run_Dir"]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"  ✓ pocket_assignments_all.csv  ({len(all_rows)} compounds)")


def build_pocket_summary_rows(pockets, all_libraries, min_pocket_size, coord_stats=None):
    """Build sorted, filtered per-pocket summary rows.

    Sorted by Total_Members descending (dominant pocket first).
    Filtered to exclude pockets below min_pocket_size.

    coord_stats: if provided by --deep-stats, adds Spread_A and Radius_A columns.
    Currently always None.

    Returns (rows list, sorted library names list) — separated from the writer
    so the table can be printed to stdout and written to CSV from the same data.
    """
    sorted_libs = sorted(all_libraries)
    rows = []

    for pid, lib_counts in pockets.items():
        total = sum(lib_counts.values())
        if total < min_pocket_size:
            continue

        row = {"Pocket": pid, "Total_Members": total}
        for lib in sorted_libs:
            row[lib] = lib_counts.get(lib, 0)

        if coord_stats is not None:
            stats = coord_stats.get(pid, {})
            row["Spread_A"] = round(stats.get("spread", float("nan")), 3)
            row["Radius_A"] = round(stats.get("radius", float("nan")), 3)

        rows.append(row)

    rows.sort(key=lambda r: r["Total_Members"], reverse=True)
    return rows, sorted_libs


def write_pocket_summary(pockets, all_libraries, output_dir, min_pocket_size, coord_stats=None):
    """Write per-pocket summary CSV and print a quick table to stdout.

    Primary deliverable for Jeremy — shows how many compounds from each library
    land in each pocket. Pocket IDs are stable across runs (frozen at CHEMBL seed).

    Columns: Pocket, Total_Members, [Library1], [Library2], ..., [LibraryN]
    Optional (--deep-stats): Spread_A, Radius_A
    """
    rows, sorted_libs = build_pocket_summary_rows(
        pockets, all_libraries, min_pocket_size, coord_stats
    )

    fieldnames = ["Pocket", "Total_Members"] + sorted_libs
    if coord_stats is not None:
        fieldnames += ["Spread_A", "Radius_A"]

    output_path = output_dir / "pocket_summary.csv"
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"  ✓ pocket_summary.csv  ({len(rows)} pockets, min_size≥{min_pocket_size})")

    # Stdout preview table
    col_w = 15
    header = f"  {'Pocket':<10} {'Total':<10}" + "".join(f"{lib:<{col_w}}" for lib in sorted_libs)
    divider = f"  {'-' * (len(header) - 2)}"
    print(f"\n  Pocket Summary (top 20):")
    print(divider)
    print(header)
    print(divider)
    for row in rows[:20]:
        line = f"  {row['Pocket']:<10} {row['Total_Members']:<10}"
        line += "".join(f"{row[lib]:<{col_w}}" for lib in sorted_libs)
        print(line)
    if len(rows) > 20:
        print(f"  ... and {len(rows) - 20} more pockets (see pocket_summary.csv)")


def write_merge_summary(args, run_dirs, summaries, libraries, all_rows, pockets, output_dir):
    """Write a programmatic audit record of this merge operation.

    Captures which runs were merged, their resolved library labels, original CIF
    source directories, clustering parameters, and compound counts at merge time.
    Fully reproducible and auditable even if run directories are later moved,
    because all relevant paths and metadata are snapshotted here.
    """
    library_counts = defaultdict(int)
    for row in all_rows:
        library_counts[row["Library"]] += 1

    merge_summary = {
        "timestamp": datetime.now().isoformat(),
        "n_runs_merged": len(run_dirs),
        "n_compounds_total": len(all_rows),
        "n_pockets_total": len(pockets),
        "min_pocket_size_filter": args.min_pocket_size,
        "library_counts": dict(library_counts),
        "runs": [
            {
                # Snapshot of resolved path at merge time
                "run_dir_at_merge": str(Path(d).resolve()),
                "library": lib,
                # Original Boltz2 CIF source — most stable identifier
                "cif_input_dir": s.get("cif_input_dir"),
                "ligand_centers_path": s.get("ligand_centers_path"),
                "seed_dir": s.get("seed_dir"),
                "mode": s.get("mode"),
                "run_timestamp": s.get("timestamp"),
                "n_input": s.get("n_input"),
                "n_assigned_existing": s.get("n_assigned_existing"),
                "n_assigned_new": s.get("n_assigned_new"),
                "params": s.get("params"),
            }
            for d, s, lib in zip(run_dirs, summaries, libraries)
        ]
    }

    output_path = output_dir / "merge_summary.json"
    with open(output_path, "w") as f:
        json.dump(merge_summary, f, indent=2)
    print(f"  ✓ merge_summary.json")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Pocket Assignment Merge")
    print("=" * 60)
    print(f"Runs:    {len(args.runs)}")
    print(f"Output:  {args.output}")
    print()

    # Load all runs
    print("Loading runs:")
    all_rows, all_summaries, all_libraries = [], [], []
    for run_dir in args.runs:
        rows, summary, library = load_run(run_dir)
        all_rows.extend(rows)
        all_summaries.append(summary)
        all_libraries.append(library)

    print(f"\n  Total compounds: {len(all_rows)}")
    print(f"  Libraries:       {', '.join(sorted(set(all_libraries)))}")

    # Validate before writing anything
    print("\nValidating:")
    validate_seed_consistency(all_summaries, args.runs)
    print("  ✓ Seed consistency — all seeded runs reference the same anchor")
    check_duplicate_structures(all_rows)

    # Aggregate counts by pocket and library
    pockets, unique_libraries = aggregate_by_pocket(all_rows)

    # Coordinate stats — None until --deep-stats is implemented
    coord_stats = calculate_coordinate_stats(pockets, args.runs)

    # Write outputs
    print("\nWriting outputs:")
    write_assignments_all(all_rows, output_dir)
    write_pocket_summary(pockets, unique_libraries, output_dir,
                         args.min_pocket_size, coord_stats)
    write_merge_summary(args, args.runs, all_summaries, all_libraries,
                        all_rows, pockets, output_dir)

    print()
    print("=" * 60)
    print("✓ Merge complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
