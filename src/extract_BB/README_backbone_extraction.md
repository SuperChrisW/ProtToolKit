# Backbone Coordinate Extraction Tools

This repository contains Python scripts for extracting backbone coordinates from protein structure files in both PDB and mmCIF formats.

## Features

- **Dual Format Support**: Works with both PDB (.pdb) and mmCIF (.cif, .mmcif) files
- **Flexible Atom Selection**: Extract standard backbone atoms (N, CA, C, O) or specify custom atoms
- **Chain Filtering**: Extract coordinates from specific chains or all chains
- **Multiple Output Formats**: Save results as CSV, JSON, or PDB format
- **Missing Atom Handling**: Gracefully handles missing atoms in structures
- **Detailed Statistics**: Provides summary information about extracted coordinates
- **Programmatic API**: Use as a Python library or command-line tool

## Installation

1. Install required dependencies:
```bash
pip install -r requirements.txt
```

Or install manually:
```bash
pip install biopython numpy pandas
```

2. Make the script executable (optional):
```bash
chmod +x extract_backbone_coordinates.py
```

## Usage

### Command Line Interface

#### Basic Usage
```bash
# Extract all backbone atoms from all chains
python extract_backbone_coordinates.py structure.pdb

# Extract from mmCIF file
python extract_backbone_coordinates.py structure.cif
```

#### Advanced Options
```bash
# Specify output file and format
python extract_backbone_coordinates.py structure.pdb --output results.csv --format csv

# Extract specific chains
python extract_backbone_coordinates.py structure.pdb --chains A,B,C

# Extract specific atoms
python extract_backbone_coordinates.py structure.pdb --atoms N,CA,C

# Extract CA atoms from chain A only, save as JSON with summary
python extract_backbone_coordinates.py structure.pdb --chains A --atoms CA --format json --summary
```

#### Command Line Options
- `input_file`: Input PDB or mmCIF file (required)
- `--output, -o`: Output file path (default: `input_basename_backbone.csv`)
- `--format, -f`: Output format - csv, json, or pdb (default: csv)
- `--chains, -c`: Comma-separated chain IDs to extract (default: all chains)
- `--atoms, -a`: Comma-separated atom names to extract (default: N,CA,C,O)
- `--summary, -s`: Print summary statistics

### Programmatic Usage

```python
from extract_backbone_coordinates import BackboneExtractor

# Initialize extractor
extractor = BackboneExtractor()

# Load structure
structure = extractor.load_structure("protein.pdb")

# Extract backbone coordinates
coordinates = extractor.extract_backbone_coordinates(structure)

# Convert to pandas DataFrame for analysis
df = extractor.coordinates_to_dataframe(coordinates)

# Save in different formats
extractor.save_coordinates(coordinates, "output.csv", "csv")
extractor.save_coordinates(coordinates, "output.json", "json")
extractor.save_coordinates(coordinates, "backbone_only.pdb", "pdb")

# Print summary
extractor.print_summary(coordinates)
```

## Output Formats

### CSV Format
Contains flattened coordinate data with columns:
- `chain_id`: Chain identifier
- `residue_number`: Residue sequence number
- `residue_name`: Three-letter residue code
- `insertion_code`: Insertion code (if any)
- `atom_name`: Atom name (N, CA, C, O, etc.)
- `x`, `y`, `z`: Cartesian coordinates
- `bfactor`: B-factor (temperature factor)
- `occupancy`: Atomic occupancy

### JSON Format
Hierarchical structure organized by chain → residue → atom:
```json
{
  "A": [
    {
      "residue_number": 1,
      "residue_name": "MET",
      "insertion_code": " ",
      "atoms": {
        "N": {
          "coordinates": [x, y, z],
          "bfactor": 20.0,
          "occupancy": 1.0
        },
        "CA": { ... },
        ...
      }
    },
    ...
  ],
  ...
}
```

### PDB Format
Standard PDB format containing only the extracted backbone atoms.

## Examples

See `example_backbone_extraction.py` for comprehensive usage examples including:
- Basic coordinate extraction
- Selective extraction (specific chains/atoms)
- DataFrame analysis and statistics
- Handling missing atoms
- Center of mass calculations

## File Structure

```
script/
├── extract_backbone_coordinates.py    # Main extraction script
├── example_backbone_extraction.py     # Usage examples
├── requirements.txt                   # Python dependencies
└── README_backbone_extraction.md      # This documentation
```

## Technical Details

### Backbone Atoms
The script extracts standard protein backbone atoms by default:
- **N**: Backbone nitrogen
- **CA**: Alpha carbon
- **C**: Carbonyl carbon
- **O**: Carbonyl oxygen

### Chain Handling
- Extracts coordinates from all models (takes first model if multiple)
- Skips hetero residues (water, ligands, etc.) by default
- Handles insertion codes and missing atoms gracefully

### Missing Atoms
Missing atoms are represented with `null`/`None` coordinates and can be:
- Filtered out in CSV output
- Preserved as null values in JSON output
- Skipped in PDB output

## Error Handling

The script includes comprehensive error handling for:
- Invalid file formats
- Missing files
- Corrupted structures
- Missing dependencies
- Empty coordinate sets

## Performance

The script is optimized for:
- Large protein structures (tested with structures >10,000 residues)
- Batch processing multiple files
- Memory-efficient coordinate storage
- Fast parsing using BioPython

## Dependencies

- **BioPython**: PDB/mmCIF parsing and structure handling
- **NumPy**: Numerical operations and coordinate arrays
- **Pandas**: Data manipulation and CSV output

## License

This script is provided as-is for research and educational purposes.

## Contributing

Feel free to submit issues or pull requests for improvements. 