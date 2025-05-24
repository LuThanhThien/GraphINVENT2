import os
import pandas as pd
from rdkit import Chem
from argparse import ArgumentParser

def parse_args_preprocess():
    """
    Parse command line arguments for the script.
    """
    parser = ArgumentParser(description="Filter SMILES strings based on atom types and charges.")
    parser.add_argument("--input_file", type=str, required=True, help="Path to the input CSV file.")
    parser.add_argument("--output_file", type=str, required=True, help="Path to the output .smi file.")
    parser.add_argument("--smi_column_name", type=str, default='SMILES', help="Column name for SMILES strings in the input CSV.")
    parser.add_argument("--overwrite", type=bool, default=True, help="Overwrite existing output file without prompt.")

    return parser.parse_args()

def filter_smi(can_smiles):
    """
    Converts a canonical SMILES string to a non-canonical form.
    """
    # 1. Parse it into an RDKit Mol
    mol = Chem.MolFromSmiles(can_smiles)

    # 2. Write it back out with canonicalization turned off
    plain_smiles = Chem.MolToSmiles(mol, canonical=False)

    return plain_smiles

def filter_smiles(args):
    """
    Filters SMILES strings from an input CSV file and writes valid ones to an output CSV file.

    Args:
        input_file (str): Path to the input CSV file containing SMILES strings.
        output_file (str): Path to the output CSV file for valid SMILES strings.
    """
    # Get input and output file paths from arguments
    input_file = args.input_file
    output_file = args.output_file
    smi_column_name = args.smi_column_name # smile column name in the input CSV
    overwrite = args.overwrite  # whether to disable warning options
    
    # Check if output file already exists
    if not overwrite and os.path.exists(output_file):
        # Ask user for confirmation to overwrite
        confirm = input(f"Output file '{output_file}' already exists. Overwrite? (y/n): ")
        if confirm.lower() != 'y':
            print("Operation cancelled. Output file not overwritten.")
            return
        
    # Read the input CSV file
    df = pd.read_csv(input_file)
    
    # Check if the specified SMILES column exists
    if smi_column_name not in df.columns:
        raise ValueError(f"Column '{smi_column_name}' not found in the input CSV file.")
    
    # Filter out rows with empty or null SMILES strings, and apply the filter_smi function
    valid_smiles = df[df[smi_column_name].notnull() & df[smi_column_name].str.strip().ne('')]
    valid_smiles = valid_smiles.drop_duplicates(subset=[smi_column_name])
    valid_smiles['SMILES'] = valid_smiles[smi_column_name].apply(filter_smi)
    valid_smiles['Name'] = valid_smiles.index.astype(str)
    valid_smiles = valid_smiles[['SMILES', 'Name']]

    # Write the valid SMILES strings to the output CSV file
    with open(output_file, 'w') as f:
        f.write("SMILES Name\n")
        for idx, row in valid_smiles.iterrows():
            f.write(f"{row['SMILES']} {row['Name']}\n")

    print(f"Filtered SMILES strings written to {output_file}")
    
if __name__ == "__main__":
    args = parse_args_preprocess()
    filter_smiles(args)

    