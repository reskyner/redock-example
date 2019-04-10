# Some example files and a first attempt at automated prep and docking

***1. Clone the repo, and cd into it***

***2. Build the docker container***

```
docker build -t docking .
```

***3. Run the container***

```
docker run -it docking /bin/bash
```

***4. change to the code directory and source environment***
```
cd /code
source activate docking-example
```

***5. start ipython***
```
ipython
```

***6. import the preparation and docking modules, and play with the example data. e.g.:***

```
from preparation import *
from docking import *

ligand_files = prep_ligand(docking_dir='/code/example_data/NUDT5A-x0114_1/', ligand_sdf='NUDT5A-x0114_1.sdf') 
protein_file = prep_protein(docking_dir='/code/example_data/NUDT5A-x0114_1/', protein_pdb='NUDT5A-x0114_1.pdb')

prepare_vina_job(docking_dir='/code/example_data/NUDT5A-x0114_1/', prepared_receptor=protein_file,    prepared_ligand=ligand_files[0], vina_exe='/vina/autodock_vina_1_1_2_linux_x86/bin/vina', box_size=[35,35,35] job_fname='vina_new.sh', job_name='vina_test')
```

## 7. run (or look at) the docking job you just generated

```
/code/example_data/NUDT5A-x0114_1/vina_new.sh
``

