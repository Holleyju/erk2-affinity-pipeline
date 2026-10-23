# erk2-affinity-pipeline

Turns raw BindingDB-style binding data into a clean table with cleaned
affinities, pKi, and RDKit descriptors.

Part 1 of a bigger project I'm working on. The eventual goal is to see
if molecular properties (weight, logP, TPSA, H-bond counts) can predict
how strongly something binds to ERK2.

**Caveat:** `data/test_molecules.csv` is a toy dataset. The molecules are
real, but I made up the affinity numbers to test the pipeline. Nothing
here is real ERK2 data yet.

## Files

| File | What it is |
|---|---|
| `erk2_pipeline/functions.py` | The three pipeline functions |
| `erk2_pipeline/__init__.py` | Marks `erk2_pipeline` as a Python package |
| `erk2_pipeline/config.py` | Empty for now, a place for settings later |
| `setup.py` | Lets you install the package with `pip install -e .` |
| `environment.yml` | Conda env: Python 3.11, pandas, numpy, RDKit, Jupyter, pytest |
| `data/test_molecules.csv` | 10 fake-ish rows with messy values on purpose |
| `notebook/test_notebook.ipynb` | Notebook where I poke at each function |

## Setup

```bash
git clone https://github.com/Holleyju/erk2-affinity-pipeline.git
cd erk2-affinity-pipeline
conda env create -f environment.yml
conda activate erk2-pipeline
pip install -e .
```

`pip install -e .` installs `erk2_pipeline` as a package, so you can import it
from any folder. The `-e` means changes to the code take effect without reinstalling.

Then open `notebook/test_notebook.ipynb`, pick the `erk2-pipeline` kernel, and Run All.

## The functions

**`clean_affinity(value)`** takes one raw text value, gives back `(number, sign)`.
Deals with spaces, commas in numbers like `"1,500"`, and the `>` / `<` qualifiers.
Blank / garbage / zero / negative returns `(nan, None)`.

**`compute_descriptors(smiles)`** is an RDKit wrapper that returns mol weight,
logP, TPSA, and H-bond donor/acceptor counts. Bad SMILES returns all NaN.

**`build_feature_table(df, exact_only=False, verbose=True)`** runs both of
the above over a whole dataframe.
- Needs `smiles` and `affinity_nM` columns.
- Adds `affinity_sign`, `pKi`, and the 5 descriptor columns.
- Drops rows that can't be used, prints a summary of what got dropped.
- `exact_only=True` also drops `>` / `<` rows.

pKi = 9 - log10(Ki in nM). Higher means tighter binding.

## Quick example

```python
import pandas as pd
from erk2_pipeline.functions import build_feature_table

raw = pd.read_csv("data/test_molecules.csv", dtype={"affinity_nM": str})
features = build_feature_table(raw)
```

Run this from the repo folder. From inside `notebook/`, use
`"../data/test_molecules.csv"` instead.

On the toy data this should print 10 rows in, 6 rows out (4 with `exact_only=True`).

## Stuff to watch out for

- **Commas = thousands separators.** So `"4,2"` would be read as 42, not 4.2.
  Fine for BindingDB (which uses periods) but would break on European data.
- **`>` and `<` values are kept by default.** Their pKi is a bound, not a real
  value, so filter on `affinity_sign` before doing any modeling.
- **Assumes nM already.** Real BindingDB exports have separate Ki / Kd / IC50
  columns, so I'll need to add a column-picking step before this works on real data.
- **RDKit prints red errors** for bad SMILES. That's normal, the row still gets dropped.
- Only tested on the toy data so far.