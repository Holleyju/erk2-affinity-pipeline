"""Clean BindingDB affinity values for analysis"""

import numpy as np
import pandas as pd

def clean_affinity(value):
    """Convert affinity values to numbers and keep any inequality signs.
    """
    if pd.isna(value):
        return np.nan, None

    text = str(value).strip()
    if not text:
        return np.nan, None

    qualifier = "="
    if text.startswith((">", "<")):
        qualifier = text[0]
        text = text[1:].strip()

    # Assumes "," is a thousands separator (BindingDB uses "." for decimals).
    text = text.replace (",","")

    try:
        number = float(text)
    except ValueError:
        return np.nan, None

    if number <= 0:
        return np.nan, None

    return number, qualifier