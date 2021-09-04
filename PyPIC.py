from Packages.lib import *
from tabulate import tabulate
import argparse
import sys

def parse_args() :
    '''
    This function parses and return arguments passed in
    '''

    desc = "Calculating INTRAPROTEIN INTERACTIONS."

    exmp = "Example: python PyPIC -PDB 3i40.pdb --hydrophobics "

    parser = argparse.ArgumentParser(prog='PyPIC', description=desc, epilog=exmp)

    parser.add_argument('--PDB', metavar='File.pdb', dest='pdbfile', required=True,
                        help='An PDB file.')

    parser.add_argument('--hydrophobics', metavar='H', dest='hydrophobic', type=str,
                        help='Checking for Hydrophobic Interaction')

    parser.add_argument('--disulphide', metavar='Disulphide-Bridges', dest='disulph', type=str,
                        help='Checking for Disulphide Bridges')

    parser.add_argument('--option', metavar='option', dest='intra', type=str, default="both", choices=['intra', 'inter', 'both'],
                        help='Which type of interactions? (intra, inter, both). e.g : For "inter" chain (e.g. hydrophobics interaction) only interaction between two different chains will be listed. choose reverse.')
    return parser.parse_args()
args = parse_args()

def main ():
    """Start PyPIC......parse options......"""
    if args.pdbfile[-4:] == '.pdb':
        try : 
            model = parse_PDB(args.pdbfile)
        except : 
            sys.stderr.write("Error attempting to load : file does not exist :  %s\n" % args.pdbfile )
            sys.exit(1)
    else:
        sys.stderr.write("Error : structure file extension not recognized, it needs to be .pdb \n")
        sys.exit(1)

    inter_hydr = hydrophobic(model)
    print('Hydrophobic Interactions within 5 Angstroms')
    print(tabulate(inter_hydr, headers=['Residue', 'Position', 'Chain','Residue', 'Position', 'Chain', 'Distance (A)'], tablefmt='orgtbl'))

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt! \n")
        sys.exit(0)