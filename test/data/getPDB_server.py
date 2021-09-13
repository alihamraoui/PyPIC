import os
import sys
from urllib import request

def download_pdb(pdbfn):
  pdbfn = pdbfn[:-2]
  pdbfn = pdbfn + ".pdb"
  downloadurl="https://files.rcsb.org/download/"
  url = downloadurl + pdbfn
  #outfnm = os.path.join(datadir, pdbfn)
  try:
      request.urlretrieve(url,  pdbfn)
  except Exception as err:
      print(str(err), file=sys.stderr)

with open ('PDBs_metafolds.txt','r') as infile:
  for struc in infile:
      download_pdb(struc)