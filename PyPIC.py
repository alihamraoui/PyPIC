from numpy.lib.shape_base import tile
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
            sys.stderr.write("Error attempting to load  %s\n" % args.pdbfile )
            sys.exit(1)
    else:
        sys.stderr.write("Error : structure file extension not recognized, it needs to be .pdb \n")
        sys.exit(1)

    table_hydrophobic = hydrophobic(model)
    headers=['Residue', 'Position', 'Chain','Residue', 'Position', 'Chain', 'Distance (A)']
    title = ['\n ' + '\033[1m' + 'Hydrophobic Interactions within 5 Angstroms' + '\033[0m']
    table = tabulate(table_hydrophobic, headers, floatfmt=(".2f"), tablefmt="fancy_grid")


    """ Output """

    print(tabulate([],headers=title, tablefmt="pipe"))
    print(table)

    """ """

    """ Repporting """

    #reporting(table)

    """ """


    tsv = True
    if tsv : 
        try : 
            tsv_table = tabulate(table_hydrophobic, headers, floatfmt=(".2f"), tablefmt="tsv")
            text_file=open("test/results/table.tsv","w")
            text_file.write(tsv_table)
            text_file.close()
        except : 
            sys.stderr.write("Error in writing tsv file")

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt! \n")
        sys.exit(0)