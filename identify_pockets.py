#!/usr/bin/env python3

"""
Identify binding pockets via DBSCAN clustering of ligand centroids.

Usage:
    python3 identify_pockets.py -i ligand_centers.txt -o output_dir [-t 5.0]

"""


import argparse
import sys
import numpy as np
from pathlib import Path

def parse_args():
    """Parse command line arguments of ligand centroids and coordinates."""
    parser = argparse.ArgumentParser(
        description="Identify binding pockets by clustering ligand coordinates"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input file containing ligand centroids (e.g. ligand_centers.txt)"
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
        filepath: Path to ligand_centers.txt

    Returns:
        ligand_names: List of structure names
        coords: numpy array of (x, y, z) coordinates

    """
    ligand_names = []
    coords = []

    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split()
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
    print(f"DEBUG, labels after reassigning noise: {labels}")

    n_clusters = len(np.unique(labels))
    print(f"Identified {n_clusters} binding pocket(s)")



if __name__ == "__main__":
    main()
