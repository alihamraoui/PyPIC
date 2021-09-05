
from Bio.PDB.PDBParser import PDBParser
import numpy as np

def parse_PDB(file):
    """
    This function parse pdb file and return Biopython pdb object
    input : pdb file.
    output : Bio.pdb object contain pdb structure.
    """
    parser=PDBParser()
    structure = parser.get_structure('tmp',file)
    model = structure[0]
    return model 

def calc_dis (atom1,atom2):
    """
    This function calculate distance between two atoms
    input : atoms object of class Bio.pdb.
    output : distance between atoms (Angchtrom).
    """
    no_inter_atom = ["CB", "CD", "CD1", "CD2", "CE", "CE1", "CH2",
                    "CE2", "CE3", "CG", "CG1", "CG2","CZ", "CZ2", "CZ3", "SD"]

    dist = np.linalg.norm(atom1.get_coord() - atom2.get_coord())
    if (atom1.get_name() in no_inter_atom and
        atom2.get_name() in no_inter_atom ) :
        if dist<5 :
            return dist


def hydrophobic(model):
    """
    This function calculate hydrophobic interactions
    input : model = Bio.pdb object contain pdb structure.
    output : list of residues interacted and distance between them.
    """
    amino = {"hydrophobic": [ "ALA", "VAL", "LEU",  "ILE", "MET", "PHE",  "TRP",  "PRO", "TYR"]}

    res_hydrophobic = []
    for chain in model.get_list():
        for residue in chain.get_list():
            if residue.get_resname() in amino['hydrophobic']:
                res_hydrophobic.append(residue)

    intra = True
    inter_hydr = []
    for i, res1 in enumerate(res_hydrophobic):
        for res2 in res_hydrophobic[i + 1:]:
            for atom1 in res1:
                for atom2 in res2:
                    dist = calc_dis (atom1,atom2)
                    if dist : 
                        if intra :
                            if res1.parent.id == res2.parent.id :
                              inter_hydr.append([res1.resname, res1.id[1], res1.parent.id, res2.resname, res2.id[1] , res2.parent.id, dist])
                        else : 
                            inter_hydr.append([res1.resname, res1.id[1], res1.parent.id, res2.resname, res2.id[1] , res2.parent.id, dist])
                        break
                else : 
                    continue
                break
    return inter_hydr