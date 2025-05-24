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
    parser.add_argument("--smi_column_name", type=str, default='SMILES', help="Column name for SMILES strings in the input CSV.")
    parser.add_argument("--included_atoms", type=str, default="{}",
                        help="Set of allowed atoms as a string, e.g. \"{'C','N','O'}\"")
    parser.add_argument("--included_charges", type=str, default="{}",
                        help="Set of allowed charges as a string, e.g. \"{-1,0,1}\"")
    parser.add_argument("--derivative", type=str, default="{}",
                        help="Set of allowed derivatives as a string, e.g. benzene smi \"{'c1ccccc1'}\"")
    parser.add_argument("--max_n_nodes", type=int, default=0)
    parser.add_argument("--min_n_nodes", type=int, default=0)
    parser.add_argument("--max_n_molecules", type=int, default=0)
    parser.add_argument("--overwrite", type=bool, default=True, help="Overwrite existing output file without prompt.")
    parser.add_argument("--split", type=bool, default=False, help="Split to train, valid and test sets.")

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
    if config.max_n_nodes > 0:
        # check max number of nodes
        if mol.GetNumAtoms() > config.max_n_nodes:
            return None
    
    # check min number of nodes
    if config.min_n_nodes > 0:
        # check min number of nodes
        if mol.GetNumAtoms() < config.min_n_nodes:
            return None

    # check if derivative of a known molecule
    if config.derivative:
        is_derivative = False  # assume it's not a derivative
        for derivative in config.derivative:
            dev_mol = Chem.MolFromSmiles(derivative)
            if dev_mol is None:
                print(f"Warning: Invalid derivative SMILES '{derivative}'")
                continue
            if mol.HasSubstructMatch(dev_mol):
                is_derivative = True
                break  # no need to check further

        if not is_derivative:
            return None  # exclude if not derivative of any

    # export as non-canonical SMILES
    out = Chem.MolToSmiles(mol, canonical=False)
    
    return out

def read_file(input_file, smi_column_name):
    """
    Reads a CSV file and returns a DataFrame.
    
    Args:
        input_file (str): Path to the input CSV file.
        smi_column_name (str): Column name for SMILES strings in the input CSV.
    
    Returns:
        pd.DataFrame: DataFrame containing the data from the CSV file.
    """
    
    if input_file.endswith('.csv'):
        df = pd.read_csv(input_file)
        # Check if the specified SMILES column exists
        if smi_column_name not in df.columns:
            raise ValueError(f"Column '{smi_column_name}' not found in the input CSV file.")
        # Filter out rows with empty or null SMILES strings
        df = df[df[smi_column_name].notnull() & df[smi_column_name].str.strip().ne('')]
        df['SMILES'] = df[smi_column_name]
    elif input_file.endswith('.smi'):
        # Read the file as a space-delimited file
        df = pd.read_csv(input_file, sep=' ')
        
    else:
        raise ValueError("Unsupported file format. Please provide a .csv or .smi file.")
    
    df = df[df['SMILES'].notna()]
    df = df.drop_duplicates(subset=['SMILES'])
    return df


def filter_to_smi(config):
    """
    Filters SMILES strings from an input CSV file and writes valid, non-canonical ones
    to an output .smi file (space-delimited SMILES + Name).
    
    Args:
        input_file (str): Path to the input CSV (expects columns "SMILES" and "Name").
        output_file (str): Path to the output .smi file.
    """
    input_file          = config.input_file
    output_file         = config.output_file
    overwrite           = args.overwrite  # whether to disable warning options
    smi_column_name     = config.smi_column_name  # smile column name in the input CSV
    
    print(f"Configuration: {config.__dict__}")
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    
    # Check if output file is .smi
    if not output_file.endswith('.smi'):
        raise ValueError("Output file must be a .smi file.")
    
    # Check if output file already exists
    if not overwrite and os.path.exists(output_file):
        # Ask user for confirmation to overwrite
        confirm = input(f"Output file '{output_file}' already exists. Overwrite? (y/n): ")
        if confirm.lower() != 'y':
            print("Operation cancelled. Output file not overwritten.")
            return

    # read file and check if required columns exist
    print(f"Reading input file: {input_file}")
    df = read_file(input_file, smi_column_name)
    if 'SMILES' not in df.columns or 'Name' not in df.columns:
        raise ValueError("Input CSV must contain 'SMILES' and 'Name' columns.")
    print(f"Number of molecules before filtering: {len(df)}")
    
    # filter out rows with empty or null SMILES strings
    df['filtered_smiles'] = df['SMILES'].apply(lambda x: filter_smi(x, config))
    print(f"Number of molecules after filtering: {len(df)}")

    # write out to .smi
    print(f"Writing results to: {output_file}")
    idx = 0
    with open(output_file, 'w') as f:
        f.write("SMILES Name\n")
        for _, row in df.iterrows():
            f.write(f"{row['filtered_smiles']} {idx}\n")
            idx += 1
            if config.max_n_molecules > 0 and idx >= config.max_n_molecules:
                break

    print("Done.")

if __name__ == "__main__":
    args = parse_args_preprocess()
    filter_to_smi(args)
    