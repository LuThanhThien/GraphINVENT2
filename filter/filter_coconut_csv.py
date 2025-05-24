import pandas as pd
from rdkit import Chem

def to_smiles_original(can_smiles):
    """
    Converts a canonical SMILES string to a non-canonical form.
    """
    # 1. Parse it into an RDKit Mol
    mol = Chem.MolFromSmiles(can_smiles)

    # 2. Write it back out with canonicalization turned off
    plain_smiles = Chem.MolToSmiles(mol, canonical=False)

    return plain_smiles
    
def filter_smiles(input_file, output_file):
    """
    Filters SMILES strings from an input CSV file and writes valid ones to an output CSV file.

    Args:
        input_file (str): Path to the input CSV file containing SMILES strings.
        output_file (str): Path to the output CSV file for valid SMILES strings.
    """
    # Read the input CSV file
    print(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file)
    valid_smiles = df[df['canonical_smiles'].notnull() & df['canonical_smiles'].str.strip().ne('')]
    print(f"Number of valid SMILES: {len(valid_smiles)}")

    # Convert canonical SMILES to non-canonical form
    print("Converting canonical SMILES to non-canonical form...")
    valid_smiles['SMILES'] = valid_smiles['canonical_smiles'].apply(to_smiles_original)
    
    # Get only smiles column, add new name column to the right as index
    valid_smiles['Name'] = valid_smiles.index.astype(str)
    valid_smiles = valid_smiles[['SMILES', 'Name']]
    
    # Write the valid SMILES strings to the output CSV file
    with open(output_file, 'w') as f:
        f.write("SMILES Name\n")
        for index, row in valid_smiles.iterrows():
            f.write(f"{row['SMILES']} {row['Name']}\n")
    
if __name__ == "__main__":
    input_file = "/root/projects/GraphINVENT2/data/pre-training/coconut/raw/coconut_complete-10-2024.csv"
    output_file = "/root/projects/GraphINVENT2/data/pre-training/coconut/raw/coconut_complete-10-2024.smi"
    
    filter_smiles(input_file, output_file)
    print(f"Filtered SMILES strings written to {output_file}")

    