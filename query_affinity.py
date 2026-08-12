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
    Sorted by most negative mean_adj_pred — stronger predicted binding first.

INPUT MODES (mutually exclusive):
    --ids           One or more compound stem IDs on the command line
    --ids-file      Text file with one compound stem ID per line
    --pocket-csv    pocket_assignments_all.csv from merge_pocket_runs.py
                    Use with --pocket to filter to one or more pockets

OUTPUT:
    CSV keyed by compound stem ID. Joins to pocket_assignments_all.csv
    on the Structure column (strip _aligned.cif suffix to get the stem).
    Ranked within each pocket by mean_adj_pred (most negative first).

Usage examples:

    # Query specific compounds across multiple library directories
    python3 query_affinity.py \\
        --input \\
            /data/Shared/Ibrahim/Debarun/IL11Ra_monomer/CHEMBL \\
            /data/Shared/Ibrahim/Debarun/IL11Ra_monomer/ApexBio \\
        --ids ApexBio_333 Drugbank_11072 FDA_65 \\
        --output genistein_IL11Ra.csv

    # Rank all compounds in pockets 1 and 2 by adj_pred
    python3 query_affinity.py \\
        --input \\
            /data/Shared/Ibrahim/Debarun/IL11Ra_monomer/CHEMBL \\
            /data/Shared/Ibrahim/Debarun/IL11Ra_monomer/ApexBio \\
        --pocket-csv merged/IL11Ra/pocket_assignments_all.csv \\
        --pocket 1 2 \\
        --output pocket1_2_IL11Ra_ranked.csv

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

# Output CSV columns — ID joins to Structure in pocket_assignments_all.csv
OUTPUT_FIELDS = [
    "Structure",
    "Pocket",
    "Library",
    "Rank_In_Pocket",
    "affinity_pred_value",
    "affinity_probability_binary",
    "adj_pred_0",
]


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Query Boltz2 affinity values and compute adj_pred scores",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--input",
        nargs="+",
        required=True,
        metavar="BOLTZ2_DIR",
        help="One or more Boltz2 output directories to search for affinity JSONs. "
             "Each directory is searched in order — first match wins. "
             "Pass multiple to cover both CHEMBL and ApexBio in one query. "
             "e.g. --input .../IL11Ra_monomer/CHEMBL .../IL11Ra_monomer/ApexBio"
    )

    # Compound ID input — mutually exclusive modes
    id_group = parser.add_mutually_exclusive_group(required=True)
    id_group.add_argument(
        "--ids",
        nargs="+",
        metavar="ID",
        help="One or more compound stem IDs "
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
        help="pocket_assignments_all.csv from merge_pocket_runs.py. "
             "Use with --pocket to filter to specific pockets."
    )

    parser.add_argument(
        "--pocket",
        nargs="+",
        type=int,
        default=None,
        metavar="N",
        help="One or more pocket IDs to filter when using --pocket-csv. "
             "Compounds are ranked within each pocket separately. "
             "e.g. --pocket 1 2 25"
    )

    parser.add_argument(
        "--output",
        required=True,
        metavar="OUTPUT_CSV",
        help="Output CSV file path"
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# ID collection
# ---------------------------------------------------------------------------

def collect_ids_from_args(args):
    """Collect compound stem IDs from whichever input mode was specified.

    Returns a list of (stem_id, pocket_label) tuples.
    pocket_label is None when not using --pocket-csv.
    """
    if args.ids:
        return [(id_.strip(), None, None) for id_ in args.ids]

    if args.ids_file:
        ids_path = Path(args.ids_file)
        if not ids_path.exists():
            print(f"Error: --ids-file not found: {args.ids_file}", file=sys.stderr)
            sys.exit(1)
        with open(ids_path) as f:
            return [(line.strip(), None, None) for line in f if line.strip()]

    if args.pocket_csv:
        return collect_ids_from_pocket_csv(args.pocket_csv, args.pocket)

    # Should never reach here due to mutually exclusive group being required
    print("Error: one of --ids, --ids-file, or --pocket-csv is required",
          file=sys.stderr)
    sys.exit(1)


def collect_ids_from_pocket_csv(pocket_csv, pocket_filter):
    """Read compound IDs from pocket_assignments_all.csv.

    Strips _aligned.cif suffix from Structure column to recover the stem ID
    that matches Boltz2 output filenames.

    pocket_filter: list of pocket IDs to include, or None for all pockets.
    """
    csv_path = Path(pocket_csv)
    if not csv_path.exists():
        print(f"Error: --pocket-csv not found: {pocket_csv}", file=sys.stderr)
        sys.exit(1)

    results = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)

        required = {"Structure", "Pocket"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            print(f"Error: pocket_assignments_all.csv missing columns: {missing}",
                  file=sys.stderr)
            sys.exit(1)

        for row in reader:
            pocket_id = int(row["Pocket"])

            if pocket_filter is not None and pocket_id not in pocket_filter:
                continue

            # Strip _aligned.cif suffix to get stem matching Boltz2 filenames
            stem = row["Structure"].replace("_affinity_model_0_aligned.cif", "")
            library = row.get("Library", "")
            results.append((stem, pocket_id, library))

    if not results:
        filter_str = f" with --pocket {pocket_filter}" if pocket_filter else ""
        print(f"Warning: no compounds found in {pocket_csv}{filter_str}",
              file=sys.stderr)

    return results


# ---------------------------------------------------------------------------
# Affinity JSON discovery and parsing
# ---------------------------------------------------------------------------

def build_affinity_index(base_dirs):
    """Pre-index all affinity JSONs across all base directories.

    Scans each base_dir once using glob, building a dict:
        stem_id -> Path to affinity JSON

    This avoids per-compound glob scans which are extremely slow at scale.
    One glob per directory instead of one glob per compound.

    Stem extraction from filename: affinity_{stem}_affinity.json
    e.g. affinity_ApexBio_333_affinity.json -> ApexBio_333

    Warns on duplicate stems (same compound in multiple batches).
    """
    index = {}
    total = 0

    for base_dir in base_dirs:
        pattern = "**/predictions/*_affinity/affinity_*_affinity.json"
        matches = list(Path(base_dir).glob(pattern))
        total += len(matches)

        for match in matches:
            # filename: affinity_{stem}_affinity.json
            fname = match.stem                  # affinity_ApexBio_333_affinity
            stem = fname[len("affinity_"):]     # ApexBio_333_affinity
            if stem.endswith("_affinity"):
                stem = stem[:-len("_affinity")]        # ApexBio_333

            if stem in index:
                print(f"  ⚠ Duplicate stem '{stem}' — keeping first", file=sys.stderr)
            else:
                index[stem] = match

    print(f"  Indexed {total} affinity JSONs → {len(index)} unique stems")
    return index


def find_affinity_json(index, stem_id):
    """ Look up a compound stem ID in the pre-built affinity index.

    0(1) dict lookup instead of per-compound glob scan.
    Returns (Path, source_dir_name) or (None, None) if not found.
    """
    path = index.get(stem_id)
    if path is None:
        return None, None
    return path, Path(path).parts[-6]

def parse_affinity_json(json_path):
    """Parse Boltz2 affinity JSON and compute adj_pred for all three models.

    adj_pred (Jeremy's adjusted prediction score):
        adj_pred = affinity_pred_value * affinity_probability_binary

    This combines the raw predicted affinity with Boltz2's confidence that
    the prediction is in the active range. More negative = stronger predicted
    binding. Computed for each of the three Boltz2 model replicates.
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
        result[f"affinity_pred_value{suffix}"]         = pred_val
        result[f"affinity_probability_binary{suffix}"] = prob_val

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

def query_compounds(affinity_index, id_pocket_pairs, libraries_label):
    """Query affinity values for a list of (stem_id, pocket_label) pairs.

    Searches across all provided base_dirs for each compound.
    Compounds not found are included with null values — never silently dropped.

    Returns a list of result dicts, one per compound.
    """
    results = []
    n_found = 0
    n_missing = 0

    null_row_template = {
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
    }

    for stem_id, pocket_label, library in id_pocket_pairs:
        json_path, source_dir = find_affinity_json(affinity_index, stem_id)

        if json_path is None:
            print(f"  ⚠ Not found: {stem_id}", file=sys.stderr)
            n_missing += 1
            results.append({
                "Structure": stem_id,
                "Pocket": pocket_label,
                "Library": library or "",
                "Rank_In_Pocket": None,
                **null_row_template,
            })
            continue

        n_found += 1
        affinity = parse_affinity_json(json_path)
        if affinity is None:
            results.append({
                "Structure": stem_id,
                "Library": library or "",
                "Pocket": pocket_label,
                "Rank_In_Pocket": None,
                **null_row_template,
            })
            continue

        results.append({
            "Structure": stem_id,
            "Library": library or "",
            "Pocket": pocket_label,
            "Rank_In_Pocket": None,
            "affinity_pred_value":           affinity.get("affinity_pred_value"),
            "affinity_probability_binary":   affinity.get("affinity_probability_binary"),
            "adj_pred_0":                    affinity.get("adj_pred_0"),
            "affinity_pred_value1":          affinity.get("affinity_pred_value1"),
            "affinity_probability_binary1":  affinity.get("affinity_probability_binary1"),
            "adj_pred_1":                    affinity.get("adj_pred_1"),
            "affinity_pred_value2":          affinity.get("affinity_pred_value2"),
            "affinity_probability_binary2":  affinity.get("affinity_probability_binary2"),
            "adj_pred_2":                    affinity.get("adj_pred_2"),
            "mean_adj_pred":                 affinity.get("mean_adj_pred"),
            "affinity_json_path":            str(json_path),
        })

    print(f"  Found: {n_found}  Missing: {n_missing}")
    return results


def assign_ranks(results):
    """Assign Rank_In_Pocket sorted by adj_pred_0 ascending (most negative first).

    Ranking is per pocket when pocket labels are present, otherwise global.
    Compounds with None adj_pred_0 are ranked last within their pocket.
    """
    by_pocket = {}
    for row in results:
        key = row.get("Pocket")
        by_pocket.setdefault(key, []).append(row)

    for pocket_rows in by_pocket.values():
        pocket_rows.sort(key=lambda r: (
            r["adj_pred_0"] is None,
            r["adj_pred_0"] if r["adj_pred_0"] is not None else float("inf")
        ))
        for rank, row in enumerate(pocket_rows, start=1):
            row["Rank_In_Pocket"] = rank

    # Sort List in place and return
    results.sort(key=lambda r: (
        r["adj_pred_0"] is None,
        r["adj_pred_0"] if r["adj_pred_0"] is not None else float("inf")
        ))

    return results


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def write_output(results, output_path):
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    print(f"\n✓ Saved {len(results)} rows to {output_path}")

    # results is already sorted by adj_pred_0 from assign_ranks
    ranked = [r for r in results if r["adj_pred_0"] is not None]
    unranked = [r for r in results if r["adj_pred_0"] is None]

    col = [30, 8, 12, 6, 10]
    hdr = f"  {'Structure':<{col[0]}} {'Pocket':<{col[1]}} {'Library':<{col[2]}} {'Rank':<{col[3]}} {'adj_pred_0':<{col[4]}}"
    div = f"  {'-' * (sum(col) + 4)}"

    def fmt_row(row):
        adj0   = f"{row['adj_pred_0']:.4f}" if row['adj_pred_0'] is not None else "N/A"
        pocket = str(row['Pocket']) if row['Pocket'] is not None else "—"
        rank   = str(row['Rank_In_Pocket']) if row['Rank_In_Pocket'] is not None else "—"
        lib    = str(row['Library'] or "")
        return (f"  {row['Structure']:<{col[0]}} {pocket:<{col[1]}} "
                f"{lib:<{col[2]}} {rank:<{col[3]}} {adj0:<{col[4]}}")

    print(f"\n  Results (sorted by adj_pred_0, most negative first):")
    print(f"  {len(ranked)} with scores, {len(unranked)} without (N/A — no affinity JSON)")
    print(div); print(hdr); print(div)

    show = 15
    if len(ranked) <= show * 2:
        for r in ranked:
            print(fmt_row(r))
    else:
        for r in ranked[:show]:
            print(fmt_row(r))
        print(f"  ... {len(ranked) - show*2} more ...")
        for r in ranked[-show:]:
            print(fmt_row(r))

    print(div)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # Validate all input directories exist
    for d in args.input:
        if not Path(d).exists():
            print(f"Error: --input directory not found: {d}", file=sys.stderr)
            sys.exit(1)

    libraries_label = "+".join(Path(d).name for d in args.input)

    print("=" * 60)
    print("Boltz2 Affinity Query")
    print("=" * 60)
    print(f"Searching:  {libraries_label}")
    print(f"Output:     {args.output}")
    if args.pocket:
        print(f"Pockets:    {args.pocket}")
    print()

    print("Indexing affinity JSONs...")
    affinity_index = build_affinity_index(args.input)
    print()

    # Collect compound IDs
    print("Collecting compound IDs...")
    id_pocket_pairs = collect_ids_from_args(args)
    print(f"  {len(id_pocket_pairs)} compounds to query")
    print()

    # Query affinity values across all input dirs
    print("Querying affinity JSONs:")
    results = query_compounds(affinity_index, id_pocket_pairs, libraries_label)

    # Sort and rank within each pocket
    results = assign_ranks(results)

    # Write output
    write_output(results, args.output)

    print()
    print("=" * 60)
    print("✓ Complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
