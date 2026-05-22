#!/usr/bin/env python3
# Written by Chris Yan
"""
Query Boltz2 affinity prediction values for a list of compounds, optionally
filtered and ranked by pocket assignment.

PIPELINE POSITION: Downstream analysis layer. Reads Boltz2 affinity JSON
outputs and pocket_assignments_all.csv produced by merge_pocket_runs.py.
Does not modify any pipeline files.

BOLTZ2 AFFINITY OUTPUT (per compound):
    predictions/{ID}_affinity/affinity_{ID}_affinity.json
    Fields: affinity_pred_value, affinity_probability_binary (x3 models)

ADJ_PRED (Jeremy's adjusted prediction score):
    adj_pred = affinity_pred_value * affinity_probability_binary
    Computed for each of the three Boltz2 models (0, 1, 2).
    Mean adj_pred across all three models is also reported.
    Sort by most negative mean_adj_pred for ranking.

INPUT MODES:
    --ids           Space-separated compound ID stems on the command line
    --ids-file      One compound ID stem per line in a text file
    --pocket-csv    pocket_assignments_all.csv from merge_pocket_runs.py
                    Use with --pocket to filter to a specific pocket,
                    and --receptor to filter to a specific receptor

OUTPUT:
    CSV keyed by compound stem ID — joins cleanly to pocket_assignments_all.csv
    on the Structure column (strip _aligned.cif suffix to get the stem).

Usage examples:

    # Single compound across all predictions in a receptor dir
    python3 query_affinity.py \\
        --input /data/Shared/Ibrahim/Debarun/IL11Ra_monomer/ApexBio \\
        --ids ApexBio_333 Drugbank_11072 \\
        --output results.csv

    # All compounds in pocket 1 ranked by adj_pred (IL11Ra)
    python3 query_affinity.py \\
        --input /data/Shared/Ibrahim/Debarun/IL11Ra_monomer/ApexBio \\
        --pocket-csv /data/Shared/ibrahim_pocket_analysis/production_20260512/merged/IL11Ra/pocket_assignments_all.csv \\
        --pocket 1 \\
        --receptor ApexBio \\
        --output pocket1_IL11Ra_ranked.csv

    # From a file of IDs
    python3 query_affinity.py \\
        --input /data/Shared/Ibrahim/Debarun/IL11Ra_monomer/ApexBio \\
        --ids-file my_compounds.txt \\
        --output results.csv
"""

import argparse
import csv
import json
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Schema — field names in Boltz2 affinity JSON.
# Update here if upstream Boltz2 output format changes.
# ---------------------------------------------------------------------------

# Three model replicates produced by Boltz2
AFFINITY_MODELS = [
    ("affinity_pred_value",  "affinity_probability_binary"),   # model 0 (primary)
    ("affinity_pred_value1", "affinity_probability_binary1"),  # model 1
    ("affinity_pred_value2", "affinity_probability_binary2"),  # model 2
]

# Output CSV columns
OUTPUT_FIELDS = [
    "ID",
    "Receptor_Dir",
    "Pocket",
    "Rank_In_Pocket",
    "affinity_pred_value",
    "affinity_probability_binary",
    "adj_pred_0",
    "affinity_pred_value1",
    "affinity_probability_binary1",
    "adj_pred_1",
    "affinity_pred_value2",
    "affinity_probability_binary2",
    "adj_pred_2",
    "mean_adj_pred",
    "affinity_json_path",
]


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Query Boltz2 affinity values and compute adj_pred scores",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Input: Boltz2 base directory
    parser.add_argument(
        "--input",
        required=True,
        metavar="BOLTZ2_DIR",
        help="Boltz2 output base directory to search for affinity JSONs "
             "(e.g. /data/Shared/Ibrahim/Debarun/IL11Ra_monomer/ApexBio)"
    )

    # Compound ID input — mutually exclusive modes
    id_group = parser.add_mutually_exclusive_group()
    id_group.add_argument(
        "--ids",
        nargs="+",
        metavar="ID",
        help="One or more compound stem IDs to query "
             "(e.g. ApexBio_333 Drugbank_11072 FDA_65)"
    )
    id_group.add_argument(
        "--ids-file",
        metavar="FILE",
        help="Text file with one compound stem ID per line"
    )
    id_group.add_argument(
        "--pocket-csv",
        metavar="FILE",
        help="pocket_assignments_all.csv from merge_pocket_runs.py — "
             "use with --pocket and optionally --receptor to filter compounds"
    )

    # Pocket / receptor filters (used with --pocket-csv)
    parser.add_argument(
        "--pocket",
        type=int,
        default=None,
        metavar="N",
        help="Filter pocket_assignments_all.csv to this pocket ID"
    )
    parser.add_argument(
        "--receptor",
        default=None,
        metavar="LIBRARY",
        help="Filter pocket_assignments_all.csv to this library/receptor label "
             "(e.g. ApexBio, CHEMBL)"
    )

    parser.add_argument(
        "--output",
        required=True,
        metavar="OUTPUT_CSV",
        help="Output CSV file path"
    )
    parser.add_argument(
        "--top",
        type=int,
        default=None,
        metavar="N",
        help="Only output top N compounds by mean_adj_pred (most negative first)"
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# ID collection
# ---------------------------------------------------------------------------

def collect_ids_from_args(args):
    """Collect compound stem IDs from whichever input mode was specified.

    Returns a list of (stem_id, pocket_label) tuples.
    pocket_label is None when not filtering by pocket.
    """
    if args.ids:
        return [(id_.strip(), None) for id_ in args.ids]

    if args.ids_file:
        ids_path = Path(args.ids_file)
        if not ids_path.exists():
            print(f"Error: --ids-file not found: {args.ids_file}", file=sys.stderr)
            sys.exit(1)
        with open(ids_path) as f:
            return [(line.strip(), None) for line in f if line.strip()]

    if args.pocket_csv:
        return collect_ids_from_pocket_csv(args.pocket_csv, args.pocket, args.receptor)

    print("Error: one of --ids, --ids-file, or --pocket-csv is required", file=sys.stderr)
    sys.exit(1)


def collect_ids_from_pocket_csv(pocket_csv, pocket_filter, receptor_filter):
    """Read compound IDs from pocket_assignments_all.csv.

    Strips _aligned.cif suffix from Structure column to recover the stem ID
    that matches Boltz2 output filenames.

    pocket_filter: if set, only include compounds assigned to this pocket ID
    receptor_filter: if set, only include compounds from this library
    """
    csv_path = Path(pocket_csv)
    if not csv_path.exists():
        print(f"Error: --pocket-csv not found: {pocket_csv}", file=sys.stderr)
        sys.exit(1)

    results = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)

        required = {"Structure", "Pocket", "Library"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            print(f"Error: pocket_assignments_all.csv missing columns: {missing}",
                  file=sys.stderr)
            sys.exit(1)

        for row in reader:
            pocket_id = int(row["Pocket"])
            library = row["Library"]

            if pocket_filter is not None and pocket_id != pocket_filter:
                continue
            if receptor_filter is not None and library != receptor_filter:
                continue

            # Strip _aligned.cif suffix — this is the stem used in Boltz2 filenames
            stem = row["Structure"].replace("_aligned.cif", "")
            results.append((stem, pocket_id))

    if not results:
        filters = []
        if pocket_filter is not None:
            filters.append(f"pocket={pocket_filter}")
        if receptor_filter:
            filters.append(f"receptor={receptor_filter}")
        print(f"Warning: no compounds found in {pocket_csv}"
              + (f" with filters: {', '.join(filters)}" if filters else ""),
              file=sys.stderr)

    return results


# ---------------------------------------------------------------------------
# Affinity JSON discovery and parsing
# ---------------------------------------------------------------------------

def find_affinity_json(base_dir, stem_id):
    """Find the affinity JSON for a compound stem ID in the Boltz2 output tree.

    Boltz2 output structure:
        {base_dir}/**/predictions/{stem_id}_affinity/affinity_{stem_id}_affinity.json

    Returns Path to the JSON file, or None if not found.
    Uses glob to handle arbitrary nesting (Nut_0000, Nut_0001, etc.).
    """
    pattern = f"**/predictions/{stem_id}_affinity/affinity_{stem_id}_affinity.json"
    matches = list(Path(base_dir).glob(pattern))

    if not matches:
        return None
    if len(matches) > 1:
        # Duplicate compound across batches — take first, warn
        print(f"  ⚠ Warning: {stem_id} found in {len(matches)} batch dirs — "
              f"using first: {matches[0]}", file=sys.stderr)
    return matches[0]


def parse_affinity_json(json_path):
    """Parse a Boltz2 affinity JSON and compute adj_pred for all three models.

    adj_pred (Jeremy's adjusted prediction score):
        adj_pred = affinity_pred_value * affinity_probability_binary

    Computed for each of the three Boltz2 model replicates.
    mean_adj_pred is the arithmetic mean across all three.

    Returns a flat dict of all values, or None if the file is malformed.
    """
    try:
        with open(json_path) as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        print(f"  Error reading {json_path}: {e}", file=sys.stderr)
        return None

    result = {}
    adj_preds = []

    for i, (pred_key, prob_key) in enumerate(AFFINITY_MODELS):
        pred_val = data.get(pred_key)
        prob_val = data.get(prob_key)

        suffix = "" if i == 0 else str(i)
        result[f"affinity_pred_value{suffix}"]        = pred_val
        result[f"affinity_probability_binary{suffix}"] = prob_val

        # adj_pred — None if either component is missing
        if pred_val is not None and prob_val is not None:
            adj = pred_val * prob_val
            result[f"adj_pred_{i}"] = round(adj, 6)
            adj_preds.append(adj)
        else:
            result[f"adj_pred_{i}"] = None

    result["mean_adj_pred"] = round(sum(adj_preds) / len(adj_preds), 6) \
        if adj_preds else None

    return result


# ---------------------------------------------------------------------------
# Core query
# ---------------------------------------------------------------------------

def query_compounds(base_dir, id_pocket_pairs, receptor_dir_label):
    """Query affinity values for a list of (stem_id, pocket_label) pairs.

    Returns a list of result dicts, one per compound.
    Compounds with no affinity JSON found are included with null values
    so they're visible in the output (not silently dropped).
    """
    results = []
    n_found = 0
    n_missing = 0

    for stem_id, pocket_label in id_pocket_pairs:
        json_path = find_affinity_json(base_dir, stem_id)

        if json_path is None:
            print(f"  ⚠ Not found: {stem_id}", file=sys.stderr)
            n_missing += 1
            # Include in output with nulls so the compound is visible
            results.append({
                "ID": stem_id,
                "Receptor_Dir": receptor_dir_label,
                "Pocket": pocket_label,
                "Rank_In_Pocket": None,
                "affinity_pred_value": None,
                "affinity_probability_binary": None,
                "adj_pred_0": None,
                "affinity_pred_value1": None,
                "affinity_probability_binary1": None,
                "adj_pred_1": None,
                "affinity_pred_value2": None,
                "affinity_probability_binary2": None,
                "adj_pred_2": None,
                "mean_adj_pred": None,
                "affinity_json_path": None,
            })
            continue

        n_found += 1
        affinity = parse_affinity_json(json_path)
        if affinity is None:
            continue

        results.append({
            "ID": stem_id,
            "Receptor_Dir": receptor_dir_label,
            "Pocket": pocket_label,
            "Rank_In_Pocket": None,          # filled in after sorting
            "affinity_pred_value":            affinity.get("affinity_pred_value"),
            "affinity_probability_binary":    affinity.get("affinity_probability_binary"),
            "adj_pred_0":                     affinity.get("adj_pred_0"),
            "affinity_pred_value1":           affinity.get("affinity_pred_value1"),
            "affinity_probability_binary1":   affinity.get("affinity_probability_binary1"),
            "adj_pred_1":                     affinity.get("adj_pred_1"),
            "affinity_pred_value2":           affinity.get("affinity_pred_value2"),
            "affinity_probability_binary2":   affinity.get("affinity_probability_binary2"),
            "adj_pred_2":                     affinity.get("adj_pred_2"),
            "mean_adj_pred":                  affinity.get("mean_adj_pred"),
            "affinity_json_path":             str(json_path),
        })

    print(f"  Found: {n_found}  Missing: {n_missing}")
    return results


def assign_ranks(results):
    """Assign Rank_In_Pocket sorted by mean_adj_pred ascending (most negative first).

    Compounds with None mean_adj_pred are ranked last.
    Ranking is per-pocket if pocket labels are present, otherwise global.
    """
    # Group by pocket
    by_pocket = {}
    for row in results:
        key = row.get("Pocket")
        by_pocket.setdefault(key, []).append(row)

    for pocket_rows in by_pocket.values():
        # Sort: valid scores first (ascending = most negative first), then None
        pocket_rows.sort(key=lambda r: (
            r["mean_adj_pred"] is None,
            r["mean_adj_pred"] if r["mean_adj_pred"] is not None else float("inf")
        ))
        for rank, row in enumerate(pocket_rows, start=1):
            row["Rank_In_Pocket"] = rank

    return results


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_output(results, output_path, top_n=None):
    """Write results CSV and print a summary table to stdout."""
    if top_n is not None:
        results = results[:top_n]

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    print(f"\n✓ Saved {len(results)} rows to {output_path}")

    # Stdout summary table
    col = [20, 8, 6, 10, 10, 10]
    header = (f"  {'ID':<{col[0]}} {'Pocket':<{col[1]}} {'Rank':<{col[2]}} "
              f"{'adj_pred_0':<{col[3]}} {'adj_pred_1':<{col[4]}} {'mean_adj':<{col[5]}}")
    divider = f"  {'-' * (sum(col) + 5)}"

    print(f"\n  Results (sorted by mean_adj_pred, most negative first):")
    print(divider)
    print(header)
    print(divider)

    for row in results[:30]:
        adj0     = f"{row['adj_pred_0']:.4f}"     if row['adj_pred_0']     is not None else "N/A"
        adj1     = f"{row['adj_pred_1']:.4f}"     if row['adj_pred_1']     is not None else "N/A"
        mean_adj = f"{row['mean_adj_pred']:.4f}"  if row['mean_adj_pred']  is not None else "N/A"
        pocket   = str(row['Pocket'])  if row['Pocket']          is not None else "—"
        rank     = str(row['Rank_In_Pocket']) if row['Rank_In_Pocket'] is not None else "—"

        print(f"  {row['ID']:<{col[0]}} {pocket:<{col[1]}} {rank:<{col[2]}} "
              f"{adj0:<{col[3]}} {adj1:<{col[4]}} {mean_adj:<{col[5]}}")

    if len(results) > 30:
        print(f"  ... and {len(results) - 30} more rows (see {output_path})")
    print(divider)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    base_dir = Path(args.input)
    if not base_dir.exists():
        print(f"Error: --input directory not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    receptor_label = base_dir.name  # e.g. 'ApexBio' from the path

    print("=" * 60)
    print("Boltz2 Affinity Query")
    print("=" * 60)
    print(f"Input dir:  {args.input}")
    print(f"Output:     {args.output}")
    if args.pocket is not None:
        print(f"Pocket:     {args.pocket}")
    if args.receptor is not None:
        print(f"Receptor:   {args.receptor}")
    print()

    # Collect compound IDs
    print("Collecting compound IDs...")
    id_pocket_pairs = collect_ids_from_args(args)
    print(f"  {len(id_pocket_pairs)} compounds to query")
    print()

    # Query affinity values
    print("Querying affinity JSONs:")
    results = query_compounds(base_dir, id_pocket_pairs, receptor_label)

    # Sort and rank
    results = assign_ranks(results)

    # Write output
    write_output(results, args.output, args.top)

    print()
    print("=" * 60)
    print("✓ Complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
