# PyPIC : Protein interaction calculator 

**PyPIC** is a tool designed to calculate non-covalent interactions in a `.pdb` file.

## Dependancies
- Python 3
- Biopython

## install Dependancies

```c
git clone https://github.com/maxibor/protinter.git
pip install -r requirements.txt
```
## How to install PyPIC

```c
cd PyPIC
chmod +x PyPIC
```

## Example

```c
PyPIC --PDB test/data/3i40.pdb --hydrophobic -csv   
```

The above example computes the hydrophobic interactions in `3i40.pdb`, displays the results, and writes them in the file `result_hydrophobic.csv`

## List of interactions available

- hydrophobic interactions
- disulphide bridges
- ionic interactions
- aromatic-aromatic interactions
- aromatic-sulphur interactions
- cation-pi interactions
- hydrogen bonds (use results with caution)

#### Binder 
## envirenement.yml
