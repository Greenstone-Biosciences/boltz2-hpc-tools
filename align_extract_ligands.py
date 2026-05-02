#!/usr/bin/env python3
"""
Align Boltz2 protein structures and extract ligand centroids.

Can be run standalone or imported as a module.

Usage:
    python align_extract_ligands.py -i INPUT_DIR -o OUTPUT_DIR [--save-aligned] [-v]
"""

import argparse
import sys
import csv
import numpy as np
from pathlib import Path


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Align structures and extract ligand centroids"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input directory containing CIF files"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory for results"
    )
    parser.add_argument(
        "--save-aligned",
        action='store_true',
        help="Save aligned CIF structures to disk"
    )
    parser.add_argument(
        "-v", "--verbose",
        action='store_true',
        help="Verbose output"
    )
    parser.add_argument(
        "--reference-cif",
        default=None,
        help="Path to a specific CIF file to use as alignment reference (overrides auto-selection)"
    )
    return parser.parse_args()


def find_cif_files(input_dir, verbose=False):
    """Find CIF files using Boltz2 structure detection.
    
    Args:
        input_dir: Input directory path
        verbose: Print detailed information
        
    Returns:
        List of CIF file paths
    """
    input_path = Path(input_dir)
    
    # Check if Boltz2 structure exists (look for Ligand_yamls directory)
    ligand_yamls_dirs = list(input_path.glob("**/Ligand_yamls"))
    
    if ligand_yamls_dirs:
        # Boltz2 structure: search within Ligand_yamls for predictions pattern
        if verbose:
            print("  Detected Boltz2 structure, using predictions pattern")
        cif_files = []
        for lyd in ligand_yamls_dirs:
            cif_files.extend(lyd.parent.glob("**/predictions/*/*.cif"))
        cif_files = sorted([str(f) for f in cif_files])
    else:
        # Flat structure: find any CIF files
        if verbose:
            print("  No Ligand_yamls detected, searching for all CIF files")
        cif_files = sorted([str(f) for f in input_path.glob("**/*.cif")])
    
    if not cif_files:
        print("Error: No .cif files found", file=sys.stderr)
        sys.exit(1)
    
    return cif_files


def align_structures(cif_files, reference_cif=None, verbose=False):
    """Align all structures to the first as reference using Cα superposition.
    
    Args:
        cif_files: List of CIF file paths
        verbose: Print alignment details
        
    Returns:
        Dictionary mapping structure names to aligned gemmi Structure objects
    """
    try:
        import gemmi
    except ImportError:
        print("Error: gemmi not installed. Install with: pip install gemmi", file=sys.stderr)
        sys.exit(1)
    
    aligned_structures = {}
    
    # Use explicit reference if provided, otherwise first file alphabetically
    ref_file = reference_cif if reference_cif else cif_files[0]
    ref_basename = Path(ref_file).name
    ref_name = Path(ref_file).stem
    
    if verbose:
        print(f"Reference structure: {ref_basename}")
    
    # Load reference structure
    ref_st = gemmi.read_structure(ref_file)
    ref_model = ref_st[0]
    
    # Get reference polymer (protein chain for alignment)
    ref_polymer = None
    for chain in ref_model:
        poly = chain.get_polymer()
        if poly:
            ref_polymer = poly
            if verbose:
                print(f"  Found polymer in chain {chain.name}: {len(poly)} residues")
            break
    
    if not ref_polymer:
        print("Error: No polymer found in reference structure", file=sys.stderr)
        sys.exit(1)
    
    # Store reference structure
    aligned_structures[ref_name] = ref_st
    
    # Print alignment table header
    print(f"\n{'Structure':<50} {'RMSD (Å)'}")
    print("-" * 65)
    print(f"{ref_basename:<50} {'0.000'}")
    
    # Align all other structures
    for mobile_file in cif_files[1:]:
        mobile_basename = Path(mobile_file).name
        mobile_name = Path(mobile_file).stem
        
        try:
            # Load mobile structure
            mobile_st = gemmi.read_structure(mobile_file)
            mobile_model = mobile_st[0]
            
            # Get mobile polymer
            mobile_polymer = None
            for chain in mobile_model:
                poly = chain.get_polymer()
                if poly:
                    mobile_polymer = poly
                    break
            
            if not mobile_polymer:
                print(f"{mobile_basename:<50} No polymer found", file=sys.stderr)
                continue
            
            # Calculate superposition on CA atoms
            ptype = gemmi.PolymerType.PeptideL
            sup = gemmi.calculate_superposition(
                ref_polymer,
                mobile_polymer,
                ptype,
                gemmi.SupSelect.CaP
            )
            
            # Apply transformation to entire model (protein + ligand)
            for chain in mobile_model:
                for residue in chain:
                    for atom in residue:
                        atom.pos = gemmi.Position(sup.transform.apply(atom.pos))
            
            print(f"{mobile_basename:<50} {sup.rmsd:.3f}")
            
            # Store aligned structure in memory
            aligned_structures[mobile_name] = mobile_st
            
        except Exception as e:
            print(f"{mobile_basename:<50} Error: {e}", file=sys.stderr)
            continue
    
    print(f"\n✓ Aligned {len(aligned_structures)} structures")
    
    return aligned_structures


def extract_ligand_centroids(aligned_structures, ligand_resname="LIG1", verbose=False):
    """Extract ligand centroids from aligned structures.
    
    Args:
        aligned_structures: Dictionary of aligned gemmi Structure objects
        ligand_resname: Residue name of ligand to extract
        verbose: Print extraction details
        
    Returns:
        List of dictionaries with structure name and centroid coordinates
    """
    ligand_data = []
    
    print("\nExtracting ligand centroids from aligned structures...")
    print(f"{'Structure':<50} {'Centroid (x, y, z)':<35} {'# Atoms'}")
    print("-" * 100)
    
    for name, structure in aligned_structures.items():
        coords = []
        
        # Extract ligand atom coordinates
        for chain in structure[0]:
            for residue in chain:
                if residue.name == ligand_resname:
                    for atom in residue:
                        coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
        
        if coords:
            # Calculate centroid
            coords_array = np.array(coords)
            centroid = coords_array.mean(axis=0)
            
            ligand_data.append({
                'structure': f"{name}_aligned.cif",
                'x': centroid[0],
                'y': centroid[1],
                'z': centroid[2],
                'n_atoms': len(coords)
            })
            
            print(f"{name}_aligned.cif"[:50].ljust(50) + 
                  f"({centroid[0]:8.3f}, {centroid[1]:8.3f}, {centroid[2]:8.3f})".ljust(35) + 
                  f"{len(coords):>10}")
        else:
            if verbose:
                print(f"{name}: No {ligand_resname} ligand found", file=sys.stderr)
    
    print(f"\n✓ Extracted {len(ligand_data)} ligand centroids")
    
    return ligand_data


def write_aligned_structures(aligned_structures, output_dir, verbose=False):
    """Write aligned structures to CIF files.
    
    Args:
        aligned_structures: Dictionary of aligned gemmi Structure objects
        output_dir: Output directory path
        verbose: Print writing details
    """
    output_path = Path(output_dir) / "aligned_cifs"
    output_path.mkdir(parents=True, exist_ok=True)
    
    if verbose:
        print(f"\nWriting aligned structures to {output_path}/")
    
    for name, structure in aligned_structures.items():
        output_file = output_path / f"{name}_aligned.cif"
        structure.make_mmcif_document().write_file(str(output_file))
        if verbose:
            print(f"  Wrote {output_file.name}")
    
    print(f"✓ Saved {len(aligned_structures)} aligned structures")


def write_centroid_csv(ligand_data, output_dir):
    """Write ligand centroid data to CSV file.
    
    Args:
        ligand_data: List of dictionaries with centroid information
        output_dir: Output directory path
    """
    output_file = Path(output_dir) / "ligand_centers.csv"
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Structure', 'X', 'Y', 'Z', 'N_atoms'])
        
        for entry in ligand_data:
            writer.writerow([
                entry['structure'],
                f"{entry['x']:.3f}",
                f"{entry['y']:.3f}",
                f"{entry['z']:.3f}",
                entry['n_atoms']
            ])
    
    print(f"✓ Ligand centers saved to ligand_centers.csv")


def main():
    """Main function."""
    args = parse_args()
    
    print("=" * 70)
    print("Structure Alignment and Ligand Extraction")
    print("=" * 70)
    print(f"Input directory:   {args.input}")
    print(f"Output directory:  {args.output}")
    print(f"Save aligned:      {args.save_aligned}")
    print("=" * 70)
    print()
    
    # Create output directory
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find CIF files
    print("Step 1: Finding CIF files...")
    cif_files = find_cif_files(args.input, args.verbose)
    print(f"  Found {len(cif_files)} CIF files")
    
    # Align structures
    print("\nStep 2: Aligning structures...")
    aligned_structures = align_structures(cif_files, args.reference_cif, args.verbose)
    
    # Extract ligand centroids
    print("\nStep 3: Extracting ligand centroids...")
    ligand_data = extract_ligand_centroids(aligned_structures, verbose=args.verbose)
    
    # Write ligand centers CSV (always)
    print()
    write_centroid_csv(ligand_data, args.output)
    
    # Conditionally write aligned structures
    if args.save_aligned:
        print()
        write_aligned_structures(aligned_structures, args.output, args.verbose)
    else:
        print("\n(Aligned structures not saved to disk)")
    
    print()
    print("=" * 70)
    print("✓ Complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
