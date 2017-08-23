import json
# Generate images of failing atoms
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from rdkit.Chem.Draw import SimilarityMaps


## Import JSON and create mol / atom indices
in_data = json.load(open("ref.json"))["data"]

## Open molecules
mols = []
for i,x in enumerate(Chem.SDMolSupplier("parm/zinc.sdf")):
    if i > 100:
        break
    mols.append(x)
out_d = {}
## Create Molecule images (atom type hilighted) with failure types as column header
for val in in_data:
    print val
    key = val[2]+">>"+val[3]
    if key not in out_d:
       out_d[key]=[]
    mol_ind = val[0]
    mol = mols[mol_ind]
    AllChem.Compute2DCoords(mol)
    atom_ind = val[1]
    highlightAtoms = [atom_ind]
    img = Draw.MolToImage(mol, highlightAtoms=highlightAtoms,highlightColor=(0.0,1.0,1.0))
    mol_file = "img/"+key+"mol"+str(mol_ind)+"_"+str(atom_ind)+".png"
    img.save(mol_file)
    out_d[key].append(mol_file)
