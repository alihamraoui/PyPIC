# PyPIC : Protein interaction calculator 

**PyPIC** is a tool designed to calculate non-covalent interactions in a `.pdb` file.
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/alihamraoui/PyPIC/HEAD) 

## Dependancies
- Python 3.9.6
- Biopython 1.77
- reduce 3.16
- tabulate 0.8.9
## install Dependancies

```c
git clone https://github.com/alihamraoui/PyPIC.git
conda env create -f environment.yml
```

## Example

```c
python PyPIC.py -pdb test/data/pdb/1c8ca.pdb  -mmh -smh -ssh -tsv 
```

The above example computes the hydrogene bond in `1c8ca.pdb`, displays the results, and writes them in the file `main_main_chain_1c8ca.tsv`

## List of interactions available

- hydrophobic interactions
- disulphide bridges
- ionic interactions
- aromatic-aromatic interactions
- aromatic-sulphur interactions
- cation-pi interactions
- hydrogen bonds 

