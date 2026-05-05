#!/usr/bin/env python3
#Written by Chris Yan
"""
Identify binding pockets via DBSCAN clustering of ligand centroids.

Usage:
    python3 identify_pockets.py -i ligand_centers.csv -o output_dir [-t 5.0]

"""


import argparse
import sys
import numpy as np
import csv
from pathlib import Path
import json
from datetime import datetime

def parse_args():
    """Parse command line arguments of ligand centroids and coordinates."""
    parser = argparse.ArgumentParser(
        description="Identify binding pockets by clustering ligand coordinates"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input file containing ligand centroids (e.g. ligand_centers.csv)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory for results"
    )
    parser.add_argument(
        "-t", "--threshold",
        type=float,
        default=5.0,
        help="Distance threshold for clustering in Angstroms (default: 5.0)"
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=1,
        help="Minimum samples per cluster (default: 1)"
    )
    parser.add_argument(
        "--save-anchors",
        default=None,
        help="Path to save pocket anchor JSON file for use in future seeded runs"
    )
    parser.add_argument(
        "--load-anchors",
        default=None,
        help="Path to pocket_anchors.json from a previous seed run for cross-run pocket assignment"
     )
    return parser.parse_args()

def read_centroids(filepath):
    """Read ligand centroids from file. File contain header, .cif of protein with ligand then x, y, z coordinates of ligand centroids.

    Args:
        filepath: Path to ligand_centers.csv

    Returns:
        ligand_names: List of structure names
        coords: numpy array of (x, y, z) coordinates

    """
    ligand_names = []
    coords = []

    try:
        with open(filepath, 'r') as f:
            reader = csv.reader(f)
            next(reader, None)
            for row in reader:
                if not row:
                    continue
                parts = row

                if len(parts) < 4:
                    continue

                ligand_names.append(parts[0])
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

        if not ligand_names:
            print("Error: No valid centroid data found", file=sys.stderr)
            sys.exit(1)

        return ligand_names, np.array(coords)

    except FileNotFoundError:
        print(f"Error: Input file '{filepath}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading centroids: {e}", file=sys.stderr)
        sys.exit(1)

def cluster_ligands(coords, eps, min_samples):
    """Cluster ligand positions using DBSCAN.

    Args:
        coords: numpy array of (x, y, z) coordinates
        eps: Maximum distance between neighbors (Angstroms)
        min_samples: Minimum samples per cluster

    Returns:
        labels: Cluster labels (0-indexed, -1 for noise)

    """
    try:
        from sklearn.cluster import DBSCAN
    except ImportError:
        print("Error: scikit-learn not installed", file=sys.stderr)
        print("Install with: pip install scikit-learn or appropriate Python environment", file=sys.stderr)
        sys.exit(1)

    clustering = DBSCAN(eps=eps, min_samples=min_samples, metric='euclidean')
    labels = clustering.fit_predict(coords)

    return labels


def reassign_noise(labels):
    """Reassign noise points labeled as (-1) to unique pocket IDs. They will be made the highest number.

    Args:
        labels: Original cluster labels

    Returns:
        new_labels: Labels with noise points reassigned
    """

    new_labels = labels.copy()
    max_label = max(labels) if len(labels) > 0 else -1

    for i, label in enumerate(new_labels):
        if label == -1:
            max_label += 1
            new_labels[i] = max_label

    return new_labels

def calculate_cluster_stats(coords, labels):
    """Calculate statistics for each cluster.

    Args:
        coords: numpy array of coordinates
        labels: Cluster labels

    Returns:
        stats: Dictionary of cluster statistics
    """
    stats = {}
    unique_labels = np.unique(labels)

    for label in unique_labels:
        mask = labels == label
        cluster_coords = coords[mask]

        # Calculate centroid (center of a cluster here)
        centroid = np.mean(cluster_coords, axis=0)

        # Calculate spread (mean distance of cluster from centroid)
        distances = np.linalg.norm(cluster_coords - centroid, axis=1)
        spread = np.mean(distances)

        stats[label] = {
            'size': np.sum(mask),
            'centroid': centroid,
            'spread': spread
        }

    return stats

def write_pocket_assignments(ligand_names, labels, output_dir, already_indexed=False):
    """Write pocket assignments to file (of receptor-ligand combos).

    Args:
        ligand_names: List of structure names
        labels: Cluster labels (0-indexed)
        output_dir: Output directory

    """
    output_path = Path(output_dir) / "pocket_assignments.csv"

    # Calculate cluster sizes
    unique_labels = np.unique(labels)
    cluster_sizes = {label: np.sum(labels==label) for label in unique_labels}

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Structure', 'Pocket', 'Cluster_Size'])
        for name, label in zip(ligand_names, labels):
            cluster_size = np.sum(labels == label)
            # Handle indexing
            pocket_num = label if already_indexed else label + 1
            writer.writerow([name, pocket_num, cluster_size]) # Change pocket IDs to 1-indexed and just an integer

def write_cluster_stats(stats, output_dir, already_indexed=False):
    """Write cluster statistics to file.

    Args:
        stats: Dictionary of cluster statistics
        output_dir: Output directory
    """
    output_path = Path(output_dir) / "cluster_statistics.csv"

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Pocket', 'N_ligands', 'Center_X', 'Center_Y', 'Center_Z', 'Spread_Angstrom'])
        for label in sorted(stats.keys()):
            pocket_id = label if already_indexed else label + 1
            size = stats[label]['size']
            cx, cy, cz = stats[label]['centroid']
            spread = stats[label]['spread']
            writer.writerow([pocket_id, size, f"{cx:.3f}", f"{cy:.3f}", f"{cz:.3f}", f"{spread:.3f}"])

def save_pocket_anchors(ligand_names, coords, labels, stats, eps, min_samples, output_path):
    """Save pocket anchor data to JSON for cross-run pocket assignment.
    
    Args:
        ligand_names: List of structure names (CHEMBL IDs)
        coords: numpy array of (x, y, z) centroids
        labels: cluster label per compound (0-indexed)
        stats: cluster stats dict from calculate_cluster_stats()
        eps: DBSCAN eps value used
        output_path: path to write pocket_anchors.json
    """
    pockets = {}
    for label in sorted(stats.keys()):
#       We skip any anchors that are < than the min_samples provided in seed
        if stats[label]['size'] < min_samples:
            continue
        pocket_id = label + 1
        mask = label == labels
        centroid = stats[label]['centroid'].tolist()

        # Radius: 100th percentile distance from centroid. Adding 2.5Å buffer
        cluster_coords = coords[mask]
        distances = np.linalg.norm(cluster_coords - stats[label]['centroid'], axis=1)
        radius = float(np.max(distances)) + 2.5 if len(distances) > 1 else float(eps) + 2.5
        # Minimum radius is eps so single-compound pockets still catch nearby new ligands
        radius = max(radius, float(eps))

        members = {}
        for name, is_member, coord in zip(ligand_names, mask, coords):
            if is_member:
                # Strip _aligned.cif suffix if present, keep ID
                clean_name = name.replace('_aligned.cif', '')
                members[clean_name] = coord.tolist()

        pockets[str(pocket_id)] = {
            'centroid': centroid,
            'radius': radius,
            'n_members': int(stats[label]['size']),
            'spread': float(stats[label]['spread']),
            'members': members
        }

    anchor_data = {
        'metadata': {
            'created': datetime.now().isoformat(),
            'eps': float(eps),
            'n_pockets': len(pockets),
            'n_compounds': len(ligand_names)
        },
        'pockets': pockets
    }

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, 'w') as f:
        json.dump(anchor_data, f, indent=2)

    # Check for overlapping spheres/pockets
    check_sphere_overlaps(pockets)

    print(f"\n✓ Pocket anchors saved to {output_path}")
    print(f"  {len(pockets)} pockets, {len(ligand_names)} compounds")
    

def load_pocket_anchors(anchor_path):
    """Load pocket anchors from previous run.

    Args:
        anchor_path: Path to anchors.json

    Returns:
        anchor_data: Full anchor dictionary
        pocket_cenroids: dict of pocket_id (int) -> np.array centroid
        pocket_radii: dict of pocket_id (int) -> float radius
    """
    with open(anchor_path, 'r') as f:
        anchor_data = json.load(f)

    pocket_centroids = {}
    pocket_radii = {}
    for pid_str, pdata in anchor_data['pockets'].items():
        pid = int(pid_str)
        pocket_centroids[pid] = np.array(pdata['centroid'])
        pocket_radii[pid] = float(pdata['radius'])

    print(f"✓ Loaded {len(pocket_centroids)} pocket anchors from {anchor_path}")
    print(f"  Original run: {anchor_data['metadata']['n_compounds']} compounds, "
          f"eps={anchor_data['metadata']['eps']} Å")
    return anchor_data, pocket_centroids, pocket_radii


def assign_to_anchors(coords, ligand_names, pocket_centroids, pocket_radii, eps, min_samples):
    """Assign new ligands to existing pockets via bounding sphere, with DBSCAN on residuals.

    Args:
        coords: numpy array of (x, y, z) centroids for new compounds
        ligand_names: list of structure names
        pocket_centroids: dict of pocket_id -> np.array centroid
        pocket_radii: dict of pocket_id -> float radius
        eps: DBSCAN eps
        min_samples: DBSCAN min_samples for residual clustering

    Returns:
        labels: array of pocket assignments (1-indexed, matching existing IDso or new)
        matched: boolean array, True if assigned to existing pocket
    """
    labels = np.full(len(coords), -1, dtype=int)
    matched = np.zeros(len(coords), dtype=bool)

    pocket_ids = sorted(pocket_centroids.keys())

    for i, coord in enumerate(coords):
        best_pid = None
        best_dist = np.inf

        for pid in pocket_ids:
            dist = np.linalg.norm(coord - pocket_centroids[pid])
            if dist <= pocket_radii[pid] and dist < best_dist:
                best_dist = dist
                best_pid = pid

        if best_pid is not None:
            labels[i] = best_pid
            matched[i] = True

#   DBSCAN on unmatched results
    unmatched_idx = np.where(~matched)[0]
    if len(unmatched_idx) > 0:
        max_existing = max(pocket_ids)
        unmatched_coords = coords[unmatched_idx]

        if len(unmatched_coords) == 1:
            # Single residual get next ID
            labels[unmatched_idx[0]] = max_existing + 1
        else:
            residual_labels = cluster_ligands(unmatched_coords, eps, min_samples)
            residual_labels = reassign_noise(residual_labels)
#           Offset new labels above existing pocket IDs
            for i, idx in enumerate(unmatched_idx):
                labels[idx] = max_existing + residual_labels[i] + 1

    n_matched = int(np.sum(matched))
    n_new = len(unmatched_idx)
    print(f"\nAssignment summary:")
    print(f"  {n_matched} compounds assigned to existing pockets")
    print(f"  {n_new} compounds assigned to new pockets")

    return labels, matched

def check_sphere_overlaps(pockets):
    """Check if any pocket bounding spheres overlap and warn.

    Args:
        pockets: dict of pocket_id -> {centroid, radius, ...}
    """
    ids = sorted(pockets.keys())
    overlaps = []
    for i, pid_a in enumerate(ids):
        for pid_b in ids[i+1:]:
            ca = np.array(pockets[pid_a]['centroid'])
            cb = np.array(pockets[pid_b]['centroid'])
            dist = np.linalg.norm(ca - cb)
            combined_radii = pockets[pid_a]['radius'] + pockets[pid_b]['radius']
            if dist < combined_radii:
                overlaps.append((pid_a, pid_b, dist, combined_radii))

    if overlaps:
        print(f"\n⚠ Warning: {len(overlaps)} overlapping pocket sphere pair(s):")
        for pid_a, pid_b, dist, combined in overlaps:
            print(f"  Pocket {pid_a} and Pocket {pid_b}: "
                f"centroid dist={dist:.2f}Å, combined radii={combined:.2f}Å")
    else:
        print("\n✓ No overlapping pocket spheres")

    return overlaps



# Write a run summary into a json
def write_run_summary(args, mode, n_input, n_assigned_existing, n_assigned_new, n_existing_pockets, n_new_pockets, seed_dir=None):
    """Write run summary JSON for programmatic parsing. Logs input and output dirs, seed dirs, what parameters were used in runs and the pockets numbers from summary."""
    summary = {
        'timestamp': datetime.now().isoformat(),
        'mode': mode,
        'input_dir': str(args.input),
        'output_dir': str(args.output),
        'seed_dir': str(seed_dir) if seed_dir else None,
        'params': {
              'eps': args.threshold,
              'min_samples': args.min_samples,
        },
        'n_input': n_input,
        'n_assigned_existing': n_assigned_existing,
        'n_assigned_new': n_assigned_new,
        'n_existing_pockets': n_existing_pockets,
        'n_new_pockets': n_new_pockets,
        'n_total_pockets': n_existing_pockets + n_new_pockets
    }
    output_path = Path(output_dir) / 'run_summary.json'
    with open(output_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"✓ Run summary saved to run_summary.json")


def main():
    """Main function"""

    # Parse arguments
    args = parse_args()

    # Reads centroids from file
    print(f"Reading centroids from: {args.input}")
    ligand_names, coords = read_centroids(args.input)
    print(f"Found {len(ligand_names)} ligands")

    # Debug lines
#    print(f"DEBUG, ligand names: {ligand_names}")
#    print(f"DEBUG, coords: {coords}")

    # Clusters or assign to anchors depending on mode
    if args.load_anchors:
        print(f"\nSeeded mode: assigning to existing pockets from {args.load_anchors}")
        anchor_data, pocket_centroids, pocket_radii = load_pocket_anchors(args.load_anchors)
        labels, matched = assign_to_anchors(
            coords, ligand_names, pocket_centroids, pocket_radii,
            args.threshold, args.min_samples
        )

        # Convert to 0-indexed for stats functions, then back
        # Labels are already 1-indexed from anchor IDs - adjust stats
        n_clusters = len(np.unique(labels))
        print(f"  {n_clusters} total pockets ({len(pocket_centroids)} existing + new)")
    else:
        print(f"\nClustering with DBSCAN (eps={args.threshold} Å, min_samples={args.min_samples})")
        labels = cluster_ligands(coords, args.threshold, args.min_samples)
        labels = reassign_noise(labels)
        n_clusters = len(np.unique(labels))
        print(f"Identified {n_clusters} binding pocket(s)")


    # Calculate stats of the clusters
    stats = calculate_cluster_stats(coords, labels)
   #  print(f"DEBUG, cluster stats: {stats}")

    # Write pocket assingment outputs and clusters stats to files
    print(f"\nWriting results to: {args.output}/")
    write_pocket_assignments(ligand_names, labels, args.output, already_indexed=bool(args.load_anchors))
    print(f" - pocket_assignments.csv")
    write_cluster_stats(stats, args.output, already_indexed=bool(args.load_anchors))
    print(f" - cluster_statistics.csv")

    # Saving pocket anchors
    if args.save_anchors:
        save_pocket_anchors(ligand_names, coords, labels, stats, args.threshold, args.min_samples,  args.save_anchors)



    # Print a summary
    print("\nPocket Summary:")
    print(f"{'Pocket':<10} {'N Ligands':<12} {'Spread (Å)':<12} {'Quality'}")
    print("-" * 50)
    for label in sorted(stats.keys()):
        pocket_id = label if args.load_anchors else label + 1
        size = stats[label]['size']
        spread = stats[label]['spread']

    # Spread qualifiers, can adjust
        if size ==1:
            quality = "○ Singleton (no clustering)"
        elif spread < 2.0:
            quality = "✓ Tight"
        elif spread < 5.0:
            quality = "~ Moderate"
        else:
            quality = "✗ Loose"

        print(f"{pocket_id:<10} {size:<12} {spread:<12.3f} {quality}")

    # Write a run summary to JSON in output dir
    write_run_summary(
        args, 
        mode='seeded' if args.load_anchors else 'fresh',
        n_input=len(ligand_names),
        n_assigned_existing=int(np.sum(matched)) if args.load_anchors else len(ligand_names),
        n_assigned_new=int(np.sum(~matched)) if args.load_anchors else 0,
        n_existing_pockets=len(pocket_centroids) if args.load_anchors else 0,
        n_new_pockets=n_clusters - (len(pocket_centroids) if args.load_anchors else 0),
        seed_dir=args.load_anchors
    )


if __name__ == "__main__":
    main()
