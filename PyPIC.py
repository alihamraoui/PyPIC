from Packages.lib import *
import argparse
import sys, os

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
                        help='Checking for main chain-main chain Hydrogen bond')

    parser.add_argument('-smh', dest='smhydrogenes', action='store_true',
                        help='Checking for side chain-main chain Hydrogen bond')

    parser.add_argument('-ssh', dest='sshydrogenes', action='store_true',
                        help='Checking for side chain-side chain Hydrogen bond')

    parser.add_argument('-tsv', dest='tsv', action='store_true',
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
            try : 
                chain = model.get_list()
                test = chain[0]
                test[2]['H']
            except :
                print("Adding Hydrogens to PDB ...")
                os.system ('reduce -version')
                cmd ='reduce ' + args.pdbfile + ' -Quiet > '+ args.pdbfile+'h'
                pdbh=os.system (cmd)
                pdbh = args.pdbfile+'h'
                model = parse_PDB(pdbh)
        except : 
            sys.stderr.write("Error attempting to load  %s\n" % args.pdbfile )
            sys.exit(1)
    else:
        sys.stderr.write("Error : structure file extension not recognized, it needs to be .pdb \n")
        sys.exit(1)
    try : 
        chain = model.get_list()
        test = chain[0]
        test[2]['H']
    except :
        cmd ='reduce ' + ' my3i40.pdb'
        os.system (cmd)



    table_mmh=[]
    table_smh=[]
    table_ssh=[]
    for chain in model.get_list():
        for res1 in chain.get_list():
            for res2 in chain.get_list():

                if args.mmhydrogenes:
                    mmh = main_main_hydrogene(res1, res2)
                    for bond in mmh :
                        table_mmh.append(bond)

                if args.smhydrogenes:
                    smh = side_main_hydrogene(res1, res2)
                    for bond in smh :
                        table_smh.append(bond)

                if args.sshydrogenes:
                    ssh = side_side_hydrogene(res1, res2)
                    for bond in ssh :
                        table_ssh.append(bond)
    if args.hydrophob :
        table_hydrophobic=hydrophobic(model)
        headers_hydrophobic=['Residue', 'Position', 'Chain','Residue', 'Position', 'Chain', 'Distance (A)']
        title_hydrophobic = ['\n ' + '\033[1m' + 'Hydrophobic Interactions within 5 Angstroms' + '\033[0m']
        output(table_hydrophobic, headers_hydrophobic, title_hydrophobic)
        if args.tsv: tsv_fun (table_hydrophobic, headers_hydrophobic, 'Hydrophobic')

    if args.mmhydrogenes:
        headers_mmh=['Residue', 'Position', 'donor', 'Chain','Residue', 'Position','acceptor', 'Chain', 'd(don-acc)', 'd(Hdon-acc)', 'agnle(don-H-acc)']
        title_mmh = ['\n ' + '\033[1m' + 'Main Chain-Main Chain Hydrogene bonds' + '\033[0m']
        output(table_mmh, headers_mmh, title_mmh)
        if args.tsv:
            outfile = 'Main_MainChaine_' + str(args.pdbfile)[-8:-4] 
            tsv_fun (table_mmh, headers_mmh, outfile)

    if args.smhydrogenes:
        headers_smh=['Residue', 'Position', 'donor', 'Chain','Residue', 'Position','acceptor', 'Chain', 'd(don-acc)', 'd(Hdon-acc)', 'agnle(don-H-acc)']
        title_smh= ['\n ' + '\033[1m' + 'Side Chain -Main Chain Hydrogene bonds' + '\033[0m']
        output(table_smh, headers_smh, title_smh)
        if args.tsv: 
            outfile = 'Side_MainChaine_' + str(args.pdbfile)[-8:-4] 
            tsv_fun (table_smh, headers_smh, outfile)

    if args.sshydrogenes:
        headers_ssh=['Residue', 'Position', 'donor', 'Chain','Residue', 'Position','acceptor', 'Chain', 'd(don-acc)', 'd(Hdon-acc)', 'agnle(don-H-acc)']
        title_ssh = ['\n ' + '\033[1m' + 'Side Chain-Side Chain Hydrogene bonds' + '\033[0m']
        output(table_ssh, headers_ssh, title_ssh)
        if args.tsv: 
            outfile = 'Side_SideChaine_' + str(args.pdbfile)[-8:-4] 
            tsv_fun (table_ssh, headers_ssh, outfile)
        

    """ Repporting """
    #reporting(table)
    """ """


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt! \n")
        sys.exit(0)