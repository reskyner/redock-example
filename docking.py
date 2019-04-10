import os

from rdkit import Chem
from rdkit.Chem import rdMolTransforms


def write_job(directory, fname, name, exe, options):
    os.chdir(directory)
    job_script = '''#!/bin/bash
cd %s
touch %s.running
%s %s > %s.log
rm %s.running
touch %s.done
''' % (directory, name, exe, options, name, name, name)

    with open(fname, 'wb') as f:
        f.write(job_script)


def prepare_vina_job(docking_dir, prepared_receptor, prepared_ligand, vina_exe, box_size, job_fname, job_name):
    os.chdir(docking_dir)
    # create an rdkit mol from ligand
    mol = Chem.MolFromMol2File(prepared_ligand)

    # get the ligand conformer and find its' centroid
    conf = mol.GetConformer()
    centre = rdMolTransforms.ComputeCentroid(conf)  # out = centre.x, centre.y and centre.z for coords

    vina_out = str(prepared_ligand.split('.')[:-1] + '_vinaout.pdbqt')

    params = [
        '--receptor',
        os.path.join(docking_dir, prepared_receptor),
        '--ligand',
        os.path.join(docking_dir, prepared_ligand),
        '--center_x',
        centre.x,
        '--center_y',
        centre.y,
        '--center_z',
        centre.z,
        '--size_x',
        str(box_size[0]),
        '--size_y',
        str(box_size[1]),
        '--size_z',
        str(box_size[2]),
        '--out',
        vina_out
    ]

    parameters = ' '.join(str(v) for v in params)

    write_job(directory=docking_dir, name=job_name, fname=job_fname, exe=vina_exe, options=parameters)
