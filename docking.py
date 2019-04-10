import os

from rdkit import Chem
from rdkit.Chem import rdMolTransforms

# from preparation import prep_ligand, prep_protein

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








def tmp_stuff():
    root_dir = luigi.Parameter()
    docking_dir = luigi.Parameter(default='comp_chem')
    ligand_pdbqt = luigi.Parameter()
    receptor_pdbqt = luigi.Parameter()
    vina_exe = luigi.Parameter(default='/dls_sw/apps/xchem/autodock_vina_1_1_2_linux_x86/bin/vina')
    box_size = luigi.Parameter(default='[40, 40, 40]')
    job_filename = luigi.Parameter(default='vina.sh')
    job_name = luigi.Parameter(default='vina')

    def requires(self):
        # 1. Prep Ligand (should have already been done for autodock)
        # 2. Prep Protein (should have already been done for autodock)
        # 3. Write vina job script
        # 4. Run vina on cluster
        # 5. Check job output
        return [PrepLigand(docking_dir=self.docking_dir, root_dir=self.root_dir,
                           ligand_sdf=self.ligand_pdbqt.replace('_prepared.pdbqt', '.sdf')),
                PrepProtein(docking_dir=self.docking_dir, root_dir=self.root_dir,
                            protein_pdb=self.receptor_pdbqt.replace('_prepared.pdbqt', '.pdb'))]

        # CheckJobOutput(job_directory=os.path.join(self.root_dir, self.docking_dir), job_output_file=out_name)]

    def run(self):
        # open ligand mol2 file (generated during PrepLigand)
        ligand = os.path.join(self.root_dir, self.docking_dir, self.ligand_pdbqt.replace('_prepared.pdbqt', '.mol2'))

        # create an rdkit mol from ligand
        mol = Chem.MolFromMol2File(ligand)

        if mol is None:
            # convert to mol with obabel
            obConv = openbabel.OBConversion()
            obConv.SetInAndOutFormats('mol2', 'mol')

            mol = openbabel.OBMol()

            # read pdb and write pdbqt
            obConv.ReadFile(mol, ligand)
            obConv.WriteFile(mol, ligand.replace('.mol2', '.mol'))

            ligand = ligand.replace('.mol2', '.mol')
            mol = Chem.MolFromMolFile(ligand)

        # get the ligand conformer and find its' centroid
        conf = mol.GetConformer()
        centre = rdMolTransforms.ComputeCentroid(conf)  # out = centre.x, centre.y and centre.z for coords

        # box size allowed for vina
        box_size = eval(self.box_size)

        # name of output file from vina
        out_name = str(self.ligand_pdbqt).replace('.pdbqt', '_vinaout.pdbqt')

        params = [
            '--receptor',
            os.path.join(self.root_dir, self.docking_dir, self.receptor_pdbqt),
            '--ligand',
            os.path.join(self.root_dir, self.docking_dir, self.ligand_pdbqt),
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
            out_name
        ]

        parameters = ' '.join(str(v) for v in params)

        write_job(job_directory=os.path.join(self.root_dir, self.docking_dir), job_filename=self.job_filename,
                  job_name=self.job_name, job_executable=self.vina_exe, job_options=parameters)

        submit_job(job_directory=os.path.join(self.root_dir, self.docking_dir), job_script=self.job_filename)

    def output(self):
        return luigi.LocalTarget(
            os.path.join(os.path.join(self.root_dir, self.docking_dir),
                         str(str(self.job_filename).replace('.sh', '.job.id'))))