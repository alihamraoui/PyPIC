from Bio.PDB.PDBParser import PDBParser
import numpy as np
from tabulate import tabulate

def prepare_parser():
    '''This function parses and return arguments passed in'''
    desc = "Calculating INTRAPROTEIN INTERACTIONS."

    exmp = "Example: PyPIC -PDB 3i40.pdb --hydrophobics "

    parser = argparse.ArgumentParser(prog='PyPIC', description=desc, epilog=exmp)

    parser.add_argument('--PDB', metavar='File.pdb', dest='pdbfile', required=True,
                        help='An PDB file.')

    parser.add_argument('--hydrophobics', metavar='H', dest='hydrophobic', type=str,
                        help='Checking for Hydrophobic Interaction')

    parser.add_argument('--disulphide', metavar='Disulphide-Bridges', dest='disulph', type=str,
                        help='Checking for Disulphide Bridges')

    parser.add_argument('--option', metavar='option', dest='intra', type=str, default="both", choices=['intra', 'inter', 'both'],
                        help='Which type of interactions? (intra, inter, both). e.g : For "inter" chain (e.g. hydrophobics interaction) only interaction between two different chains will be listed. choose reverse.')

    return parser