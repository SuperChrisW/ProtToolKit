#!/usr/bin/env python3
"""
Backbone Coordinate Extractor for PDB and mmCIF Files

This script extracts backbone coordinates (N, CA, C, O atoms) from protein structures
in both PDB and mmCIF formats, organized by chain.

Dependencies:
    - biopython: pip install biopython
    - numpy: pip install numpy
    - pandas: pip install pandas

Usage:
    python extract_backbone_coordinates.py input_file.pdb --output output.csv
    python extract_backbone_coordinates.py input_file.cif --format json
    python extract_backbone_coordinates.py input_file.pdb --chains A,B --atoms N,CA,C
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import json
import csv

try:
    from Bio.PDB import PDBParser, MMCIFParser, Structure
    from Bio.PDB.Entity import Entity
    import numpy as np
    import pandas as pd
except ImportError as e:
    print(f"Error: Missing required dependency - {e}")
    print("Please install required packages:")
    print("pip install biopython numpy pandas")
    sys.exit(1)


class BackboneExtractor:
    """Class to extract backbone coordinates from protein structures."""
    
    def __init__(self):
        self.pdb_parser = PDBParser(QUIET=True)
        self.mmcif_parser = MMCIFParser(QUIET=True)
        
        # Standard backbone atoms for proteins
        self.backbone_atoms = ['N', 'CA', 'C', 'O']
        
    def load_structure(self, file_path: str) -> Structure:
        """
        Load protein structure from PDB or mmCIF file.
        
        Args:
            file_path: Path to the structure file
            
        Returns:
            Bio.PDB.Structure object
            
        Raises:
            ValueError: If file format is not supported
            FileNotFoundError: If file doesn't exist
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        
        file_extension = Path(file_path).suffix.lower()
        structure_id = Path(file_path).stem
        
        if file_extension == '.pdb':
            return self.pdb_parser.get_structure(structure_id, file_path)
        elif file_extension in ['.cif', '.mmcif']:
            return self.mmcif_parser.get_structure(structure_id, file_path)
        else:
            raise ValueError(f"Unsupported file format: {file_extension}")
    
    def extract_backbone_coordinates(
        self, 
        structure: Structure, 
        chains: Optional[List[str]] = None,
        atoms: Optional[List[str]] = None
    ) -> Dict:
        """
        Extract backbone coordinates from the structure.
        
        Args:
            structure: Bio.PDB.Structure object
            chains: List of chain IDs to extract (default: all chains)
            atoms: List of atom names to extract (default: ['N', 'CA', 'C', 'O'])
            
        Returns:
            Dictionary with chain coordinates organized by chain ID
        """
        if atoms is None:
            atoms = self.backbone_atoms
        
        coordinates = {}
        
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                
                # Skip chain if not in specified chains list
                if chains and chain_id not in chains:
                    continue
                
                chain_coords = []
                
                for residue in chain:
                    # Skip hetero residues (water, ligands, etc.)
                    if residue.get_id()[0] != ' ':
                        continue
                    
                    residue_data = {
                        'residue_number': residue.get_id()[1],
                        'residue_name': residue.get_resname(),
                        'insertion_code': residue.get_id()[2],
                        'atoms': {}
                    }
                    
                    for atom_name in atoms:
                        if atom_name in residue:
                            atom = residue[atom_name]
                            residue_data['atoms'][atom_name] = {
                                'coordinates': atom.get_coord().tolist(),
                                'bfactor': atom.get_bfactor(),
                                'occupancy': atom.get_occupancy()
                            }
                        else:
                            # Handle missing atoms
                            residue_data['atoms'][atom_name] = {
                                'coordinates': [None, None, None],
                                'bfactor': None,
                                'occupancy': None
                            }
                    
                    chain_coords.append(residue_data)
                
                if chain_coords:  # Only add non-empty chains
                    coordinates[chain_id] = chain_coords
        
        return coordinates
    
    def coordinates_to_dataframe(self, coordinates: Dict) -> pd.DataFrame:
        """
        Convert coordinates dictionary to pandas DataFrame.
        
        Args:
            coordinates: Dictionary of chain coordinates
            
        Returns:
            pandas DataFrame with flattened coordinate data
        """
        rows = []
        
        for chain_id, chain_data in coordinates.items():
            for residue in chain_data:
                for atom_name, atom_data in residue['atoms'].items():
                    if atom_data['coordinates'][0] is not None:  # Skip missing atoms
                        row = {
                            'chain_id': chain_id,
                            'residue_number': residue['residue_number'],
                            'residue_name': residue['residue_name'],
                            'insertion_code': residue['insertion_code'],
                            'atom_name': atom_name,
                            'x': atom_data['coordinates'][0],
                            'y': atom_data['coordinates'][1],
                            'z': atom_data['coordinates'][2],
                            'bfactor': atom_data['bfactor'],
                            'occupancy': atom_data['occupancy']
                        }
                        rows.append(row)
        
        return pd.DataFrame(rows)
    
    def save_coordinates(
        self, 
        coordinates: Dict, 
        output_path: str, 
        format_type: str = 'csv'
    ) -> None:
        """
        Save coordinates to file in specified format.
        
        Args:
            coordinates: Dictionary of chain coordinates
            output_path: Output file path
            format_type: Output format ('csv', 'json', 'pdb')
        """
        if format_type.lower() == 'csv':
            df = self.coordinates_to_dataframe(coordinates)
            df.to_csv(output_path, index=False)
            
        elif format_type.lower() == 'json':
            with open(output_path, 'w') as f:
                json.dump(coordinates, f, indent=2)
                
        elif format_type.lower() == 'pdb':
            self._save_as_pdb(coordinates, output_path)
            
        else:
            raise ValueError(f"Unsupported output format: {format_type}")
    
    def _save_as_pdb(self, coordinates: Dict, output_path: str) -> None:
        """Save coordinates in PDB format (backbone atoms only)."""
        with open(output_path, 'w') as f:
            atom_id = 1
            
            for chain_id, chain_data in coordinates.items():
                for residue in chain_data:
                    res_num = residue['residue_number']
                    res_name = residue['residue_name']
                    insertion = residue['insertion_code']
                    
                    for atom_name, atom_data in residue['atoms'].items():
                        coords = atom_data['coordinates']
                        if coords[0] is not None:  # Skip missing atoms
                            bfactor = atom_data['bfactor'] or 0.0
                            occupancy = atom_data['occupancy'] or 1.0
                            
                            f.write(
                                f"ATOM  {atom_id:5d} {atom_name:^4s} {res_name:3s} "
                                f"{chain_id:1s}{res_num:4d}{insertion:1s}   "
                                f"{coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}"
                                f"{occupancy:6.2f}{bfactor:6.2f}          "
                                f"{atom_name[0]:>2s}\n"
                            )
                            atom_id += 1
            
            f.write("END\n")
    
    def print_summary(self, coordinates: Dict) -> None:
        """Print summary statistics of extracted coordinates."""
        total_residues = sum(len(chain_data) for chain_data in coordinates.values())
        total_atoms = 0
        
        for chain_data in coordinates.values():
            for residue in chain_data:
                for atom_data in residue['atoms'].values():
                    if atom_data['coordinates'][0] is not None:
                        total_atoms += 1
        
        print(f"\nExtraction Summary:")
        print(f"├── Chains: {len(coordinates)}")
        print(f"├── Residues: {total_residues}")
        print(f"├── Atoms: {total_atoms}")
        print(f"└── Chain IDs: {', '.join(coordinates.keys())}")
        
        # Per-chain statistics
        for chain_id, chain_data in coordinates.items():
            chain_atoms = sum(
                1 for residue in chain_data 
                for atom_data in residue['atoms'].values()
                if atom_data['coordinates'][0] is not None
            )
            print(f"    Chain {chain_id}: {len(chain_data)} residues, {chain_atoms} atoms")


def main():
    """Main function to handle command line interface."""
    parser = argparse.ArgumentParser(
        description="Extract backbone coordinates from PDB and mmCIF files"
    )
    
    parser.add_argument(
        'input_file',
        help='Input PDB or mmCIF file path'
    )
    
    parser.add_argument(
        '--output', '-o',
        help='Output file path (default: input_basename_backbone.csv)',
        default=None
    )
    
    parser.add_argument(
        '--format', '-f',
        choices=['csv', 'json', 'pdb'],
        default='csv',
        help='Output format (default: csv)'
    )
    
    parser.add_argument(
        '--chains', '-c',
        help='Comma-separated list of chain IDs to extract (default: all chains)',
        default=None
    )
    
    parser.add_argument(
        '--atoms', '-a',
        help='Comma-separated list of atom names to extract (default: N,CA,C,O)',
        default='N,CA,C,O'
    )
    
    parser.add_argument(
        '--summary', '-s',
        action='store_true',
        help='Print summary statistics'
    )
    
    args = parser.parse_args()
    
    # Parse chains and atoms
    chains = args.chains.split(',') if args.chains else None
    atoms = args.atoms.split(',') if args.atoms else None
    
    # Generate output filename if not provided
    if args.output is None:
        input_path = Path(args.input_file)
        output_name = f"{input_path.stem}_backbone.{args.format}"
        args.output = str(input_path.parent / output_name)
    
    try:
        # Initialize extractor and process file
        extractor = BackboneExtractor()
        
        print(f"Loading structure from: {args.input_file}")
        structure = extractor.load_structure(args.input_file)
        
        print("Extracting backbone coordinates...")
        coordinates = extractor.extract_backbone_coordinates(
            structure, chains=chains, atoms=atoms
        )
        
        if not coordinates:
            print("Warning: No coordinates extracted. Check chain IDs and file format.")
            return
        
        print(f"Saving coordinates to: {args.output}")
        extractor.save_coordinates(coordinates, args.output, args.format)
        
        if args.summary:
            extractor.print_summary(coordinates)
        
        print("✓ Extraction completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 