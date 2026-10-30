import numpy as np
import pandas as pd

from .affinity import clean_affinity
from .descriptors import compute_descriptors

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