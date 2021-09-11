
from Bio.PDB.PDBParser import PDBParser
import numpy as np
import sys
from math import degrees
from tabulate import tabulate

#main = ["N", "O", "OXT"]
side_all = ["ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ", "OD1", "OD2", "OE1",
            "OE2", "OG", "OG1", "OH", "SD", "SG"]
side_don = {"ARG": ["NE", "NH1", "NH2"], "ASN": ["ND2"], "CYS": ["SG"],
            "GLN": ["NE2"], "HIS": ["ND1", "NE2"], "LYS": ["NZ"],
            "SER": ["OG"], "THR": ["OG1"], "TRP": ["NE1"], "TYR": ["OH"]}
side_acc = {"ASN": ["OD1"], "ASP": ["OD1", "OD2"], "GLN": ["OE1"],
            "GLU": ["OE1", "OE2"], "HIS": ["ND1", "NE2"], "MET": ["SD"], "SER": ["OG"],
            "THR": ["OG1"], "TYR": ["OH"]}


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
    '''
    '''
    tabulate(table, headers, floatfmt=(".2f"), tablefmt="fancy_grid")
    print(tabulate([],headers=title, tablefmt="pipe"))
    print(tabulate(table, headers, floatfmt=(".2f"), tablefmt="orgtbl"))

def tsv_fun (table,headers,name): 
    try : 
        tsv_table = tabulate(table, headers, floatfmt=(".2f"), tablefmt="tsv")
        text_file=open("test/results/{}_table.tsv".format(name),"w")
        text_file.write(tsv_table)
        text_file.close()
    except : 
        sys.stderr.write("Error in writing tsv file  \n")


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
            if dist < d : return dist
    else :
        dist = np.linalg.norm(atom1.get_coord() - atom2.get_coord())
        if dist < d : return dist


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
                    dist = calc_dis (atom1,atom2, 5, filter=True)
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


def main_main_hydrogene(res1, res2) :
    """
    This function calculate hydrogene bonds
    input : model = Bio.pdb object contain pdb structure.
    output : list of donor and accepto, distance between them and angle between donor-H-acceptor.
    """
    mm_hydrogenes = []
    if (res1.parent.id == res2.parent.id 
        and res1.resname != 'HOH' and res2.resname != 'HOH'):
        for donor in res1:
            if donor.name == 'N':
                for accept in res2 :
                    if accept.name in ['O','OXT']:
                        dist = calc_dis (donor, accept, 3.5, filter = False)
                        if abs(res1.id[1] - res2.id[1]) >= 2 :
                            if dist :
                                if donor.parent.has_id('H'):
                                    donor_H = donor.parent['H']
                                    Dh = donor_H - accept
                                    ang = get_angle(donor, donor_H, accept)
                                    mm_hydrogenes.append([res1.resname, res1.id[1], donor.id ,res1.parent.id, res2.resname, res2.id[1] ,\
                                                    accept.id , res2.parent.id, dist, Dh, ang],)
                                else : pass
    return mm_hydrogenes


def get_hydrogene(donor):
    '''
    '''
    neib_Hydrogene=[]
    residue = donor.parent
    if donor.name == 'N':
        if residue.has_id('H'):
            neib_Hydrogene.append(residue['H'])
            return neib_Hydrogene
        else : 
            return [donor]
    else :
        atoms = donor.parent.get_list()
        if len(donor.name) > 2:
            idx = len(donor.name)-1
        else : 
            idx = len(donor.name)
        Hname = 'H'+ donor.name[1:idx]
        for a in atoms :
            if len(a.name) > 3:
                Hname = 'H'+ donor.name[1:]
                if a.name [:(idx+1)] == Hname:
                    neib_Hydrogene.append(a)
            else :
                if a.name [:idx] == Hname:
                    neib_Hydrogene.append(a)
        if len(neib_Hydrogene) != 0 : 
            return neib_Hydrogene
        else : 
            return [residue['H']]


def check_msHydrogene (donor, accept, dis):
    '''
    '''
    bond = []
    res1 = donor.parent
    res2 = accept.parent
    donor_H = get_hydrogene(donor)
    for H in donor_H : 
        Dh = H - accept
        ang = get_angle(donor, H, accept)
        bond.append([res1.resname, res1.id[1], donor.id ,res1.parent.id, res2.resname, res2.id[1] ,\
                            accept.id , res2.parent.id, dis, Dh, ang]) 
    return bond


def side_main_hydrogene (res1, res2):
    '''
    '''
    sm_hydrogenes=[]
    if res1.parent.id == res2.parent.id :
        for atom1 in res1:
            for atom2 in res2 :
                donor = atom1
                accept = atom2
                if ((donor.name == 'N' and res2.resname in side_acc and accept.name in side_acc[str(res2.resname)]) 
                or (accept.name == 'O' and donor.name in side_all)):
                    dist = calc_dis (donor, accept, 3.5, filter= False)
                    if abs(res1.id[1] - res2.id[1]) >= 2 :
                        if dist :
                            smBond=check_msHydrogene (donor,accept, dist)
                            for bond in smBond :
                                sm_hydrogenes.append(bond)
    return sm_hydrogenes


def side_side_hydrogene (res1,res2):
    '''
    '''
    ss_hydrogenes= []
    if res1.parent.id == res2.parent.id :
        for atom1 in res1:
            for atom2 in res2 :
                donor = atom1
                accept = atom2
                if (res1.resname in side_don and donor.name in side_don[str(res1.resname)] \
                    and res2.resname in side_acc and accept.name in side_acc[str(res2.resname)]):
                    dist = calc_dis (donor, accept, 3.5, filter =False)
                    if abs(res1.id[1] - res2.id[1]) >= 2 :
                        if dist :
                            smBond=check_msHydrogene (donor,accept, dist)
                            for bond in smBond :
                                ss_hydrogenes.append(bond)
    return ss_hydrogenes