import os
import pandas as pd
from rdkit import Chem
from config import dict2str
import ast
from argparse import ArgumentParser

def parse_args_preprocess():
    parser = ArgumentParser(description="Filter SMILES strings based on atom types and charges.")
    parser.add_argument("--input_file", type=str, required=True)
    parser.add_argument("--output_file", type=str, required=True)
    parser.add_argument("--included_atoms", type=str, default="{}",
                        help="Set of allowed atoms as a string, e.g. \"{'C','N','O'}\"")
    parser.add_argument("--included_charges", type=str, default="{}",
                        help="Set of allowed charges as a string, e.g. \"{-1,0,1}\"")
    parser.add_argument("--derivative", type=str, default="{}",
                        help="Set of allowed derivatives as a string, e.g. benzene smi \"{'c1ccccc1'}\"")
    parser.add_argument("--max_n_nodes", type=int, default=None)
    parser.add_argument("--min_n_nodes", type=int, default=None)
    parser.add_argument("--max_n_molecules", type=int, default=None)
    parser.add_argument("--overwrite", type=bool, default=True, help="Overwrite existing output file without prompt.")

    args = parser.parse_args()

    # Convert the string sets into actual Python sets/lists
    args.included_atoms = list(ast.literal_eval(args.included_atoms))
    args.included_charges = list(ast.literal_eval(args.included_charges))
    if args.derivative is not None:
        args.derivative = list(ast.literal_eval(args.derivative))
    else:
        args.derivative = None

    return args

def filter_smi(smiles, config):
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
    if config.included_atoms:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() not in config.included_atoms:
                return None

    # check formal charge
    if config.included_charges:
        for atom in mol.GetAtoms():
            if atom.GetFormalCharge() not in config.included_charges:
                return None
        
    # check number of nodes
    if config.max_n_nodes:
        # check max number of nodes
        if mol.GetNumAtoms() > config.max_n_nodes:
            return None
    
    # check min number of nodes
    if config.min_n_nodes:
        # check min number of nodes
        if mol.GetNumAtoms() < config.min_n_nodes:
            return None

    # check if derivative of a known molecule
    if config.derivative:
        is_derivative = False  # assume it's not a derivative
        for derivative in config.derivative:
            dev_mol = Chem.MolFromSmiles(derivative)
            if dev_mol is None:
                continue
            if mol.HasSubstructMatch(dev_mol):
                is_derivative = True
                break  # no need to check further

        if not is_derivative:
            return None  # exclude if not derivative of any

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
    overwrite = args.overwrite  # whether to disable warning options
    
    # Check if output file already exists
    # Print configuration and input/output file paths
    print(f"Configuration: {config.__dict__}")
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    
    if not overwrite and os.path.exists(output_file):
        # Ask user for confirmation to overwrite
        confirm = input(f"Output file '{output_file}' already exists. Overwrite? (y/n): ")
        if confirm.lower() != 'y':
            print("Operation cancelled. Output file not overwritten.")
            return

    # read file and check if required columns exist
    print(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file, sep=' ')
    if 'SMILES' not in df.columns or 'Name' not in df.columns:
        raise ValueError("Input CSV must contain 'SMILES' and 'Name' columns.")
    
    # filter out rows with empty or null SMILES strings
    df['filtered_smiles'] = df['SMILES'].apply(lambda x: filter_smi(x, config))
    df = df[df['filtered_smiles'].notna()]
    df = df.drop_duplicates(subset=['filtered_smiles'])
    df['Name'] = df.index.astype(str)
    df = df[['filtered_smiles', 'Name']]
    print(f" Rows after filtering: {len(df)}")

    # write out to .smi
    print(f"Writing results to: {output_file}")
    with open(output_file, 'w') as f:
        f.write("SMILES Name\n")
        for idx, row in df.iterrows():
            f.write(f"{row['filtered_smiles']} {row['Name']}\n")
            if config.max_n_molecules and idx >= config.max_n_molecules - 1:
                break

    print("Done.")

if __name__ == "__main__":
    args = parse_args_preprocess()
    filter_to_smi(args)
    