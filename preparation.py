import os

import openbabel
from htmd.ui import *


def obabel_conversion(input_file, ouput_type, options):
    # get input file type from input file extension
    ftype = input_file.split('.')[-1]

    # set up openbabel conversion and input and output types
    obConv = openbabel.OBConversion()
    obConv.SetInAndOutFormats(ftype, ouput_type)

    # apply all options
    for o in options:
        obConv.AddOption(o)

    mol = openbabel.OBMol()

    # read input, convert, and write output
    obConv.ReadFile(mol, input_file)
    obConv.WriteFile(mol, input_file.replace(str('.' + ftype), str('.' + ouput_type)))

    # return name of the converted file
    return str(input_file.replace(str('.' + ftype), str('.' + ouput_type)))


def prep_protein(docking_dir, protein_pdb):
    # make sure we're in the directory where the docking will be done
    os.chdir(docking_dir)
    # use HTMD to prepare the protein
    mol = Molecule(protein_pdb)
    prot = proteinPrepare(mol)
    prot.remove('resname WAT')

    # write the prepared protein out to a new pdb file
    prot.write(str(protein_pdb).replace('.pdb', '_prepared.pdb'))

    # set the new protein, and convert it with babel to pdbqt
    protein = os.path.join(docking_dir, str(protein_pdb).replace('.pdb', '_prepared.pdb'))
    converted_file = obabel_conversion(input_file=protein, ouput_type='pdbqt', options=['x', 'r'])

    return converted_file


def prep_ligand(docking_dir, ligand_sdf):
    # make sure we're in the directory where the docking will be done
    os.chdir(os.path.join(docking_dir))
    # set the full path of the ligand
    ligand = os.path.join(docking_dir, ligand_sdf)
    # convert to pdbqt, and then mol2
    conv1 = obabel_conversion(input_file=ligand, ouput_type='pdbqt', options=[])
    converted_file = obabel_conversion(input_file=conv1, ouput_type='mol2', options=[])

    return conv1, converted_file