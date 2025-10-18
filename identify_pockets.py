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

def write_pocket_assignments(ligand_names, labels, output_dir):
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
            writer.writerow([name, f"pocket{label + 1}", cluster_size]) # Change pocket IDs to 1-indexed

def write_cluster_stats(stats, output_dir):
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
            pocket_id = f"pocket{label + 1}"
            size = stats[label]['size']
            cx, cy, cz = stats[label]['centroid']
            spread = stats[label]['spread']
            writer.writerow([pocket_id, size, f"{cx:.3f}", f"{cy:.3f}", f"{cz:.3f}", f"{spread:.3f}"])


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

    # Clusters ligand centroids
    print(f"\nClustering with DBSCAN (eps={args.threshold} Å, min_samples={args.min_samples})")
    labels = cluster_ligands(coords, args.threshold, args.min_samples)
#    print(f"DEBUG, labels after clustering: {labels}")

    # Reassign noise points from clustering
    labels = reassign_noise(labels)
    # print(f"DEBUG, labels after reassigning noise: {labels}")

    n_clusters = len(np.unique(labels))
    print(f"Identified {n_clusters} binding pocket(s)")

    # Calculate stats of the clusters
    stats = calculate_cluster_stats(coords, labels)
   #  print(f"DEBUG, cluster stats: {stats}")

    # Write pocket assingment outputs and clusters stats to files
    print(f"\nWriting results to: {args.output}/")
    write_pocket_assignments(ligand_names, labels, args.output)
    print(f" - pocket_assignments.csv")
    write_cluster_stats(stats, args.output)
    print(f" - cluster_statistics.csv")

    # Print a summary
    print("\nPocket Summary:")
    print(f"{'Pocket':<10} {'N Ligands':<12} {'Spread (Å)':<12} {'Quality'}")
    print("-" * 50)
    for label in sorted(stats.keys()):
        pocket_id = f"pocket{label + 1}"
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



if __name__ == "__main__":
    main()
