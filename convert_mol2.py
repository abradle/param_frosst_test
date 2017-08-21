
# coding: utf-8
from __future__ import print_function
import atom_types
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from parse_mol2 import run

def pad_atom_name(name, element):

    padded = name
    if (len(element) == 1):
        if (len(name) == 2):
            padded = ' ' + name + ' '
        if (len(name) == 3):
            padded = ' ' + name
    if (len(element) == 2):
        if (len(name) == 2):
            padded = name + '  '
        if (len(name) == 3):
            padded = name + ' '
    return padded

#import pyrogen_swig as pysw
def atom_name_from_atomic_number_and_count(element, count, inc):
    name = element
    name += str(count+inc)
    return name

def add_atom_names(mol):
    nz = {}
    atom_names = []
    for atom in mol.GetAtoms():
       try:
          n = atom.GetProp('name')
          atom_names.append(n)
       except KeyError:

          # we want to generate a name that is not already in the atom_names list
          
          z = atom.GetAtomicNum()
          if z in nz:
             nz[z]  = nz[z] + 1
          else:
             nz[z] = 1;
          ele = atom.GetSymbol().upper()
          # we add inc argument, which gets added to count (nz[z]) in case that we
          # already have a name that matches (previous) return from call to 
          # atom_name_from_atomic_number_and_count()
          #
          inc = 0
          name = atom_name_from_atomic_number_and_count(ele, nz[z], inc)
          p_name = pad_atom_name(name, ele)
          # print('c.f.', name, " atom names", atom_names)
          while p_name in atom_names :
             inc += 1
             name = atom_name_from_atomic_number_and_count(ele, nz[z], inc)
             p_name = pad_atom_name(name, ele)
          # print(atom, 'made-name', p_name, ":")

          atom.SetProp("name", p_name)
          atom_names.append(p_name)
    return atom_names


def write_mols(suppl):
    for idx,mol in enumerate(suppl):
        fn = 'pf-'+str(idx)+'.sdf'
        # print('writing {}'.format(fn))
        print(Chem.MolToMolBlock(mol), file=file(fn, 'w'))

def get_amber_types(m, mol_idx, ref_types=None, check_2_chars_only=True, make_dictionaries=False):
    Chem.AllChem.ComputeGasteigerCharges(m)
    names = add_atom_names(m)
    atom_types.set_atom_types(m)
    atom_types.set_parmfrosst_atom_types(m)
    aa_types=[]
    cif_file_name = "test-pyrogen-"+str(mol_idx)+'.cif'
    if make_dictionaries:
       pysw.mmcif_dict_from_mol("XYZ", m.GetProp('_Name'), m, cif_file_name, False, False, True)
    for idx,atom in enumerate(m.GetAtoms()):
        try:
            # in pathological cases, pf_atom_type is not set
            coot_amber_type = atom.GetProp('pf_atom_type')
            # if we were passed ref_types, compare the atom type
            name = atom.GetProp('name')
            ref_type=ref_types[idx]
            if ref_type == "Nstar":
               ref_type = 'N*'
            check_coot_amber_type = coot_amber_type
            if check_2_chars_only:
                check_coot_amber_type = coot_amber_type[:2]
		# C -> Ca or Cb or Cc etc is a special case
		if coot_amber_type == "Ca":
		    check_coot_amber_type = "C"
		if coot_amber_type == "Cb":
		    check_coot_amber_type = "C"
		if coot_amber_type == "Cc":
		    check_coot_amber_type = "C"
		if coot_amber_type == "Cd":
		    check_coot_amber_type = "C"
		if coot_amber_type == "Ce":
		    check_coot_amber_type = "C"
		if coot_amber_type == "Cf":
		    check_coot_amber_type = "C"
            if check_coot_amber_type != ref_type:
                conf = m.GetConformer()
                coords = conf.GetAtomPosition(idx)
                cc = "(set-rotation-centre " + \
                     str(coords.x) + ' ' + str(coords.y) + ' ' + str(coords.z) + ')'
                s = 'failed to match mol_idx {} atom_idx {:>2} name {} coot_type {:>4} ref_type {:>2} at {}'
                print(s.format(mol_idx, idx, name, coot_amber_type, ref_type, cc))
            else:
                s = 'matched         mol_idx {} atom_idx {:>2} name {} coot_type {:>4} ref_type {}'
                print(s.format(mol_idx, idx, name, coot_amber_type, ref_type))

            p=(atom.GetProp('name'), coot_amber_type)
            aa_types.append(p)
        except TypeError:
            pass
	except KeyError:
	    # pf_atom_type was not set on an atom
	    pass
    return aa_types



if __name__ == '__main__':
   
   mdl_in = 'input.mdl'
   mol2_in = 'input.mol2'
   mol2_out = 'output.mol2'

   # Read in the input mdl file
   mol = Chem.MolFromMolFile(mdl_in,removeHs=False)
   ### Not needed
   # Add names 
   names = add_atom_names(mol)
   atom_types.set_parmfrosst_atom_types(mol)
   repls = []
   for atom in mol.GetAtoms():
      atom_type = atom.GetProp('pf_atom_type')
      repls.append(atom_type)
   run(mol2_in,repls,mol2_out) 
