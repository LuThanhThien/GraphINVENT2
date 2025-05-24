import pandas as pd
from rdkit import Chem


included_atoms = {'C', 'N', 'O', 'S'}
included_charges = {-1, 0, 1}
max_n_nodes = 50
min_n_nodes = 15
max_n_molecules = 20000

def do_filter_smiles(smiles):
    """
    Parse a SMILES string, check atom types, and return a non-canonical SMILES.
    Returns:
      - non-canonical SMILES (str) if valid and only allowed atoms present
      - None otherwise
    """
    global included_atoms, included_charges, max_n_nodes, min_n_nodes
    
    # try to build molecule
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return None

    if mol is None:
        return None

    # check atom symbols
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in included_atoms:
            return None

    # check formal charge
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() not in included_charges:
            return None
        
    # check number of nodes
    if mol.GetNumAtoms() > max_n_nodes:
        return None
    
    # check min number of nodes
    if mol.GetNumAtoms() < min_n_nodes:
        return None
    
    # export as non-canonical SMILES
    return Chem.MolToSmiles(mol, canonical=False)


def filter_to_smi(input_file, output_file):
    """
    Filters SMILES strings from an input CSV file and writes valid, non-canonical ones
    to an output .smi file (space-delimited SMILES + Name).
    
    Args:
        input_file (str): Path to the input CSV (expects columns "SMILES" and "Name").
        output_file (str): Path to the output .smi file.
    """
    global max_n_molecules
    
    print(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file, sep=' ')
    print(f" Total rows read: {len(df)}")

    # apply filter + conversion
    print("Filtering and converting SMILES...")
    df['filtered_smiles'] = df['SMILES'].apply(do_filter_smiles)

    # drop any rows where do_filter_smiles returned None
    df = df[df['filtered_smiles'].notna()]
    print(f" Rows after filtering: {len(df)}")

    # write out to .smi
    print(f"Writing results to: {output_file}")
    with open(output_file, 'w') as f:
        f.write("SMILES Name\n")
        for idx, row in df.iterrows():
            f.write(f"{row['filtered_smiles']} {row['Name']}\n")
            if idx >= max_n_molecules - 1:
                break
            
    print(f" Total rows written: {idx + 1}")
    print("Done.")

if __name__ == "__main__":
    # Example usage
    input_file = "/root/projects/GraphINVENT2/data/pre-training/minicoconut/flavonoid_new.smi"
    output_file = "/root/projects/GraphINVENT2/data/pre-training/minicoconut/flavonoid_new_filtered.smi"

    filter_to_smi(input_file, output_file)
    print(f"Filtered SMILES strings written to {output_file}")