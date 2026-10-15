import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

# Helper to clean affinity values (strings into numbers)
def clean_affinity(value):
    # Check if it's already empty
    if pd.isna(value):
        return np.nan, None
    
    # Convert to string and remove whitespace
    value = str(value).strip()
    
    # If it's empty after stripping, return NaN
    if value == "":
        return np.nan, None
    
    # Default sign is equals
    sign = "="
    
    # Check if there is a < or > at the start
    if value.startswith(">"):
        sign = ">"
        value = value[1:] # remove the sign
    elif value.startswith("<"):
        sign = "<"
        value = value[1:] # remove the sign
        
    value = value.strip()
    
    # Get rid of commas (thousands separators)
    # Note: BindingDB uses "." for decimals
    value = value.replace(",", "")
    
    try:
        num = float(value)
    except ValueError:
        return np.nan, None
    
    # Ignore negative or zero values
    if num <= 0:
        return np.nan, None
        
    return num, sign

# Get some basic properties from a SMILES string
def compute_descriptors(smiles):
    # Make a dict of NaNs to reuse if anything goes wrong
    empty = {
        "mol_weight": np.nan,
        "logp": np.nan,
        "tpsa": np.nan,
        "h_donors": np.nan,
        "h_acceptors": np.nan,
    }
    
    # Check if the input is empty or missing
    if pd.isna(smiles):
        return empty
    
    # Convert to string and remove extra spaces
    smiles_str = str(smiles).strip()
    
    # If it's just an empty string, return NaN too
    if smiles_str == "":
        return empty
    
    # Try to make the molecule
    # RDKit gives back None if the SMILES is broken
    molecule = Chem.MolFromSmiles(smiles_str)
    
    if molecule is None:
        return empty
    
    # If we got here, the molecule is valid, so calculate everything
    return {
        "mol_weight": Descriptors.MolWt(molecule),
        "logp": Descriptors.MolLogP(molecule),
        "tpsa": rdMolDescriptors.CalcTPSA(molecule),
        "h_donors": Lipinski.NumHDonors(molecule),
        "h_acceptors": Lipinski.NumHAcceptors(molecule),
    }

# Turn a raw BindingDB table into a clean feature table for analysis
# Adds: affinity_sign, pKi, mol_weight, logp, tpsa, h_donors, h_acceptors
def build_feature_table(df, exact_only=False, verbose=True):
    
    # Check that we have the columns we need
    for col in ["smiles", "affinity_nM"]:
        if col not in df.columns:
            raise ValueError("Missing required column: " + col)
    
    # Make a copy so we don't mess up the original dataframe
    table = df.copy()
    n_start = len(table)
    
    # Clean the affinity values
    # clean_affinity returns a tuple (number, sign), so we split it
    cleaned = table["affinity_nM"].apply(clean_affinity)
    table["affinity_nM"] = cleaned.apply(lambda x: x[0])
    table["affinity_sign"] = cleaned.apply(lambda x: x[1])
    
    # Convert nM to pKi
    table["pKi"] = 9 - np.log10(table["affinity_nM"])
    
    # Compute descriptors for each SMILES and add as columns
    desc_list = table["smiles"].apply(compute_descriptors).tolist()
    descriptors = pd.DataFrame(desc_list, index=table.index)
    table = pd.concat([table, descriptors], axis=1)
    
    # Find rows with problems
    bad_affinity = table["affinity_nM"].isna()
    bad_smiles = table["mol_weight"].isna()
    
    # Keep only rows where both affinity and smiles are ok
    keep = (bad_affinity == False) & (bad_smiles == False)
    n_unusable = n_start - keep.sum()
    
    # If exact_only, also drop the > and < rows
    n_censored_dropped = 0
    if exact_only:
        exact = (table["affinity_sign"] == "=")
        n_censored_dropped = (keep & (exact == False)).sum()
        keep = keep & exact
    
    result = table[keep].reset_index(drop=True)
    
    # Print a summary if verbose is True
    if verbose:
        n_censored_kept = result["affinity_sign"].isin([">", "<"]).sum()
        print("--- build_feature_table summary ---")
        print("Rows in:", n_start)
        print("Bad/missing affinity:", bad_affinity.sum())
        print("Invalid/missing SMILES:", bad_smiles.sum())
        print("Dropped (unusable):", n_unusable)
        if exact_only:
            print("Dropped (> or < values):", n_censored_dropped)
        print("Rows out:", len(result))
        print("  of which > or < values:", n_censored_kept)
    
    return result