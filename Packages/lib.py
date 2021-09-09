
from Bio.PDB.PDBParser import PDBParser
import numpy as np
from math import degrees
from tabulate import tabulate


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


def output(table,headers,title):
    tabulate(table, headers, floatfmt=(".2f"), tablefmt="fancy_grid")
    print(tabulate([],headers=title, tablefmt="pipe"))
    print(tabulate(table, headers, floatfmt=(".2f"), tablefmt="orgtbl"))


def calc_dis (atom1,atom2, d, filter):
    """
    This function calculate distance between two atoms
    input : atoms object of class Bio.pdb.
    output : distance between atoms (Angchtrom).
    """
    if filter:
        no_inter_atom = ["CB", "CD", "CD1", "CD2", "CE", "CE1", "CH2",
                    "CE2", "CE3", "CG", "CG1", "CG2","CZ", "CZ2", "CZ3", "SD"]
        if (atom1.get_name() in no_inter_atom and
            atom2.get_name() in no_inter_atom ) :
            dist = np.linalg.norm(atom1.get_coord() - atom2.get_coord())
            if dist< d : return dist
    else :
        dist = np.linalg.norm(atom1.get_coord() - atom2.get_coord())
        if dist< d : return dist



def get_angle(atom1, H, atom2):
    vecA = atom1.coord - H.coord
    vecB = atom2.coord - H.coord
    angle = np.arccos(np.dot(vecA,vecB) / (np.linalg.norm(vecA) * np.linalg.norm(vecB)))
    return degrees(angle)


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
                    dist = calc_dis (atom1,atom2,5, filter=True)
                    if dist : 
                        if intra :
                            if res1.parent.id == res2.parent.id:
                                inter_hydr.append([res1.resname, res1.id[1], res1.parent.id, res2.resname, res2.id[1] , res2.parent.id, dist])
                        else : 
                            inter_hydr.append([res1.resname, res1.id[1], res1.parent.id, res2.resname, res2.id[1] , res2.parent.id, dist])
                        break
                else : 
                    continue
                break
    return inter_hydr


def main_main_hydrogene(model) :
    """
    This function calculate hydrogene bonds
    input : model = Bio.pdb object contain pdb structure.
    output : list of donor and accepto, distance between them and angle between donor-H-acceptor.
    """
    mm_hydrogenes = []
    for chain in model.get_list():
        res_list = chain.get_list()
        for res1 in res_list:
            for res2 in res_list:
                if (res1.parent.id == res2.parent.id 
                    and res1.resname != 'HOH' and res2.resname != 'HOH'):
                    for donor in res1:
                        if donor.name == 'N':
                            for accept in res2 :
                                if accept.name == 'O':
                                    dist = calc_dis (donor, accept, 3.5, filter = False)
                                    if abs(res1.id[1] - res2.id[1]) >= 2 :
                                        if dist :
                                            try :
                                                donor_H = donor.parent['H']
                                                Dh = donor_H - accept
                                                ang = get_angle(donor, donor_H, accept)
                                                mm_hydrogenes.append([res1.resname, res1.id[1], donor.id ,res1.parent.id, res2.resname, res2.id[1] ,\
                                                                accept.id , res2.parent.id, dist, Dh, ang],)
                                            except : pass
    return mm_hydrogenes