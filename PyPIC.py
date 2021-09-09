from Packages.lib import *
import argparse
import sys

def parse_args() :
    '''
    This function parses and return arguments passed in
    '''

    desc = "PYPIC : CALCULATING INTRAPROTEIN INTERACTIONS."

    exmp = "Example: python PyPIC -PDB 3i40.pdb -hydrophobics "

    parser = argparse.ArgumentParser(prog='PyPIC', description=desc, epilog=exmp)

    parser.add_argument('-pdb', metavar='File.pdb', dest='pdbfile', required=True,
                        help='An PDB file.')

    parser.add_argument('-hydrophob', dest='hydrophob', action='store_true',
                        help='Checking for Hydrophobic Interaction')

    parser.add_argument('-disulphide', metavar='', dest='disulph', type=str,
                        help='Checking for Disulphide Bridges')

    parser.add_argument('-mmh', dest='mmhydrogenes', action='store_true',
                        help='Checking for main main Hydrogen bond')

    parser.add_argument('-tsv', metavar='', dest='tsv', type=str,
                        help='Export the result tables in tsv format')

    parser.add_argument('-option', metavar='option', dest='intra', type=str, default="both", choices=['intra', 'inter', 'both'],
                        help='Which type of interactions? (intra, inter, both). e.g : For "inter" chain (e.g. hydrophobics interaction) only interaction between two different chains will be listed. choose reverse.')
    
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    return parser.parse_args()


def main (): 
    """
    Start PyPIC......parse options......
    """
    args = parse_args()
    if args.pdbfile[-4:] == '.pdb':
        try : 
            model = parse_PDB(args.pdbfile)
        except : 
            sys.stderr.write("Error attempting to load  %s\n" % args.pdbfile )
            sys.exit(1)
    else:
        sys.stderr.write("Error : structure file extension not recognized, it needs to be .pdb \n")
        sys.exit(1)

    if args.hydrophob :
        table= hydrophobic(model)
        headers=['Residue', 'Position', 'Chain','Residue', 'Position', 'Chain', 'Distance (A)']
        title = ['\n ' + '\033[1m' + 'Hydrophobic Interactions within 5 Angstroms' + '\033[0m']
        output(table,headers,title)

    if args.mmhydrogenes:
        table= main_main_hydrogene( model)
        headers=['Residue', 'Position', 'donor', 'Chain','Residue', 'Position','acceptor', 'Chain', 'd(don-acc)', 'd(Hdon-acc)', 'agnle(don-H-acc)']
        title = ['\n ' + '\033[1m' + 'Main-Main chain Hydrogene bonds' + '\033[0m']
        output(table,headers,title)

    if args.tsv: 
        try : 
            tsv_table = tabulate(table, headers, floatfmt=(".2f"), tablefmt="tsv")
            text_file=open("test/results/table.tsv","w")
            text_file.write(tsv_table)
            text_file.close()
        except : 
            sys.stderr.write("Error in writing tsv file  \n")

    """ Repporting """
    #reporting(table)
    """ """


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt! \n")
        sys.exit(0)