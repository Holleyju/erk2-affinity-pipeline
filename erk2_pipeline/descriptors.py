import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

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
