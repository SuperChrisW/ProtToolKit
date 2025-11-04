import pandas as pd

def csv_to_fasta(csv_file, fasta_file, design_col='Design', seq_col='Sequence'):
    """
    Reads a CSV file and writes protein sequences in FASTA format.

    Args:
        csv_file (str): Path to the input CSV file.
        fasta_file (str): Path to the output FASTA file.
        design_col (str): Column name for design/identifier.
        seq_col (str): Column name for sequence.
    """
    df = pd.read_csv(csv_file)
    with open(fasta_file, 'w') as f:
        for idx, row in df.iterrows():
            design = str(row[design_col])
            seq = str(row[seq_col])
            f.write(f">{design}\n{seq}\n")

# Example usage:
# csv_to_fasta('input.csv', 'output.fasta')
if __name__ == "__main__":
    csv_file = '/Users/liyao/Documents/Xtalpi_Intern/AFdesign_proj/IL23-IL23R_binder_design/IL23_results/additional_designs_filtered.csv'
    df = pd.read_csv(csv_file)
    df['Design'] = df['src_folder'].apply(lambda x: x.replace('/','_'))+ '_' + df['Design']
    df.to_csv(csv_file, index=False)

    csv_to_fasta(csv_file, 
                csv_file.replace('.csv', '.fasta'),
                design_col='Design',
                seq_col='Sequence')