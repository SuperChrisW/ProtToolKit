import os
import argparse
from Bio.PDB import PDBParser, PDBIO, Select
from tqdm import tqdm

class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

def extract_chain_from_pdb(input_folder, output_folder, chain_id):
    """
    Extracts the specified chain from all PDB files in input_folder and saves them in output_folder
    with filenames ending in _chainX.pdb, where X is the chain_id.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    for filename in tqdm(os.listdir(input_folder), desc="Processing PDB files"):
        if not filename.endswith('.pdb'):
            continue
        input_path = os.path.join(input_folder, filename)
        output_filename = filename[:-4] + f'_chain{chain_id}.pdb'
        output_path = os.path.join(output_folder, output_filename)
        try:
            structure = parser.get_structure('structure', input_path)
            io.set_structure(structure)
            io.save(output_path, select=ChainSelect(chain_id))
        except Exception as e:
            print(f"Failed to process {filename}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract a specific chain from all PDB files in a folder using Biopython.")
    parser.add_argument("--input", type=str, help="Path to the input folder containing PDB files.")
    parser.add_argument("--output", type=str, help="Path to the output folder to save extracted chain PDB files.")
    parser.add_argument("--chain", type=str, default="B", help="Chain ID to extract (default: B).")
    args = parser.parse_args()

    extract_chain_from_pdb(args.input, args.output, args.chain)