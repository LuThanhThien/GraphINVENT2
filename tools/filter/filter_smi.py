import pandas as pd
from rdkit import Chem
from argparse import ArgumentParser
from config import dict2str, str2dict

def parse_args():
    """
    Parse command line arguments for the script.
    """
    parser = ArgumentParser(description="Filter SMILES strings based on specified criteria.")
    parser.add_argument("input_file", type=str, help="Path to the input CSV file containing SMILES strings.")
    parser.add_argument("output_file", type=str, help="Path to the output .smi file for valid SMILES strings.")
    parser.add_argument("--config", type=int, default=20000, help="Maximum number of molecules to write to output file.")
    args = parser.parse_args()
    config_dict = str2dict(args.config)
    for key, value in config_dict.items():
        setattr(args, key, value)
    return args

def do_filter_smiles(smiles, config):
    """
    Parse a SMILES string, check atom types, and return a non-canonical SMILES.
    Returns:
      - non-canonical SMILES (str) if valid and only allowed atoms present
      - None otherwise
    """
    # try to build molecule
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return None

    if mol is None:
        return None

    # check atom symbols
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in config.included_atoms:
            return None

    # check formal charge
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() not in config.included_charges:
            return None
        
    # check number of nodes
    if mol.GetNumAtoms() > config.max_n_nodes:
        return None
    
    # check min number of nodes
    if mol.GetNumAtoms() < config.min_n_nodes:
        return None

    # export as non-canonical SMILES
    out = Chem.MolToSmiles(mol, canonical=False)
    
    return out

def filter_to_smi(config):
    """
    Filters SMILES strings from an input CSV file and writes valid, non-canonical ones
    to an output .smi file (space-delimited SMILES + Name).
    
    Args:
        input_file (str): Path to the input CSV (expects columns "SMILES" and "Name").
        output_file (str): Path to the output .smi file.
    """
    input_file = config.input_file
    output_file = config.output_file
    print(f"Configuration: {dict2str(config.__dict__)}")
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    
    print(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file, sep=' ')
    print(f" Total rows read: {len(df)}")

    # apply filter + conversion
    print("Filtering and converting SMILES...")
    df['filtered_smiles'] = df['SMILES'].apply(lambda x: do_filter_smiles(x, config))

    # drop any rows where do_filter_smiles returned None
    df = df[df['filtered_smiles'].notna()]
    print(f" Rows after filtering: {len(df)}")

    # write out to .smi
    print(f"Writing results to: {output_file}")
    with open(output_file, 'w') as f:
        f.write("SMILES Name\n")
        for idx, row in df.iterrows():
            f.write(f"{row['filtered_smiles']} {row['Name']}\n")
            if idx >= config.max_n_molecules - 1:
                break

    print(f" Total rows written: {idx + 1}")
    print("Done.")

if __name__ == "__main__":
    args = parse_args()
    filter_to_smi(args)
    