import numpy as np
import pandas as pd

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