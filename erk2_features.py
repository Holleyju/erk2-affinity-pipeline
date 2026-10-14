"""Feature engineering helpers for the ERK2 binding-affinity pipeline."""

import numpy as np
import pandas as pd


def clean_affinity(value):
    """Convert a raw BindingDB affinity entry to a float in nM.

    Handles strings like "4.2", ">10000", "<0.5", " 12 " and blanks.
    Returns (number, qualifier), where qualifier is "=", ">" or "<".
    Returns (np.nan, None) if the value cannot be interpreted.
    """
    if pd.isna(value):
        return np.nan, None

    text = str(value).strip()
    if text == "":
        return np.nan, None

    qualifier = "="
    if text[0] in (">", "<"):
        qualifier = text[0]
        text = text[1:].strip()

    try:
        number = float(text)
    except ValueError:
        return np.nan, None

    if number <= 0:
        return np.nan, None

    return number, qualifier