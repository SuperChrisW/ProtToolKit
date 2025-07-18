#!/usr/bin/env python3
"""
Example usage of the BackboneExtractor class

This script demonstrates how to use the backbone coordinate extraction
functionality programmatically.
"""

import sys
from pathlib import Path

# Add the current directory to Python path to import the extractor
sys.path.append(str(Path(__file__).parent))

from extract_backbone_coordinates import BackboneExtractor


def example_basic_extraction(pdb_file):
    """Example of basic backbone coordinate extraction."""
    print("=== Basic Extraction Example ===")
    
    extractor = BackboneExtractor()
    
    # Load structure
    structure = extractor.load_structure(pdb_file)
    
    # Extract all backbone coordinates
    coordinates = extractor.extract_backbone_coordinates(structure)
    
    # Print summary
    extractor.print_summary(coordinates)
    
    # Save as CSV
    extractor.save_coordinates(coordinates, "output_basic.csv", "csv")
    print("✓ Saved basic extraction to output_basic.csv")
    
    return coordinates


def example_selective_extraction(pdb_file):
    """Example of selective extraction (specific chains and atoms)."""
    print("\n=== Selective Extraction Example ===")
    
    extractor = BackboneExtractor()
    
    # Load structure
    structure = extractor.load_structure(pdb_file)
    
    # Extract only CA atoms from specific chains
    coordinates = extractor.extract_backbone_coordinates(
        structure, 
        chains=['A', 'B'],  # Only chains A and B
        atoms=['CA']        # Only CA (alpha carbon) atoms
    )
    
    # Print summary
    extractor.print_summary(coordinates)
    
    # Save as JSON
    extractor.save_coordinates(coordinates, "output_selective.json", "json")
    print("✓ Saved selective extraction to output_selective.json")
    
    return coordinates


def example_dataframe_analysis(coordinates):
    """Example of analyzing extracted coordinates using pandas."""
    print("\n=== DataFrame Analysis Example ===")
    
    extractor = BackboneExtractor()
    
    # Convert to DataFrame for analysis
    df = extractor.coordinates_to_dataframe(coordinates)
    
    print(f"DataFrame shape: {df.shape}")
    print("\nFirst few rows:")
    print(df.head())
    
    print("\nBasic statistics:")
    print(df[['x', 'y', 'z']].describe())
    
    print("\nAtom type distribution:")
    print(df['atom_name'].value_counts())
    
    print("\nChain distribution:")
    print(df['chain_id'].value_counts())
    
    # Calculate center of mass for each chain
    print("\nCenter of mass by chain:")
    for chain_id in df['chain_id'].unique():
        chain_df = df[df['chain_id'] == chain_id]
        com_x = chain_df['x'].mean()
        com_y = chain_df['y'].mean()
        com_z = chain_df['z'].mean()
        print(f"  Chain {chain_id}: ({com_x:.2f}, {com_y:.2f}, {com_z:.2f})")


def example_missing_atoms_handling(pdb_file):
    """Example showing how missing atoms are handled."""
    print("\n=== Missing Atoms Handling Example ===")
    
    extractor = BackboneExtractor()
    structure = extractor.load_structure(pdb_file)
    
    # Try to extract some atoms that might not exist in all residues
    coordinates = extractor.extract_backbone_coordinates(
        structure, 
        atoms=['N', 'CA', 'C', 'O', 'CB']  # CB might be missing in glycine
    )
    
    # Count missing atoms
    missing_count = 0
    total_count = 0
    
    for chain_data in coordinates.values():
        for residue in chain_data:
            for atom_data in residue['atoms'].values():
                total_count += 1
                if atom_data['coordinates'][0] is None:
                    missing_count += 1
    
    print(f"Missing atoms: {missing_count}/{total_count} ({missing_count/total_count*100:.1f}%)")


def main():
    """Main example function."""
    # You can replace this with an actual PDB file path
    example_pdb = "example.pdb"
    
    print("Backbone Coordinate Extraction Examples")
    print("="*50)
    
    # Check if example file exists
    if not Path(example_pdb).exists():
        print(f"Example PDB file '{example_pdb}' not found.")
        print("Please provide a valid PDB or mmCIF file path.")
        print("\nTo use this script:")
        print("1. Download a PDB file (e.g., from https://www.rcsb.org/)")
        print("2. Update the 'example_pdb' variable with the file path")
        print("3. Run this script again")
        return
    
    try:
        # Run examples
        coords_basic = example_basic_extraction(example_pdb)
        coords_selective = example_selective_extraction(example_pdb)
        example_dataframe_analysis(coords_basic)
        example_missing_atoms_handling(example_pdb)
        
        print("\n" + "="*50)
        print("✓ All examples completed successfully!")
        
    except Exception as e:
        print(f"Error running examples: {e}")
        print("Make sure you have a valid PDB/mmCIF file and required dependencies installed.")


if __name__ == "__main__":
    main() 