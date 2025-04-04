# Usage Guide for proposed_mopac

## Files descriptions
- **main.py** - starts proposed_mopac method
- **orca_calculate_folder.py** - starts ORCA validation calculations
- **saddle.py** - starts SADDLE validation calculations
- **u_functions** - contains auxiliary functions for main
- **estimation.py** - creates plots and calculates KIE
- **new_MOPAC_calculator.py** - contains new methods for ASE, such as: isotopes labels, extraction functions, inputs modifies functions
- **orca_estimation.py** - calculates KIE for ORCA validation
## How to install proposed mopac
1. You should create conda environment. If you don't have conda, go to [install_conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
2. Copy proposed_mopac.yml file to your machine.
3. Type in linux terminal
```
conda env create -f proposed_mopac.yml
```
then, activate environment
```
conda activate proposed_mopac
```

## How to use proposed mopac
1. First, you should copy `proposed_mopac_calc` directory to your computer
2. Second, create **not_opted.xyz** file of your molecule system and rename folder `proposed_mopac_calc/mop/substance_name`
to `proposed_mopac_calc/mop/YOUR_SYSTEM_NAME`
3. Then, paste **not_opted.xyz** to `proposed_mopac_calc/mop/YOUR_SYSTEM_NAME/start`
4. Creat starting command, using example, and replace _XX_ with necessary ones.
```
python3 main.py -x opted_h2o.xyz --eps XX.X --charge X --mult XXXXXXX --step X.XX --isotope X \
--moving_atom_h XX --dest_atom XX --h_neighbour XX --place /mop/YOUR_SYSTEM_NAME/calc --threads 100
```
**necessary parameters**
- **--eps** - dielectric constant, for water 78.4
- **--charge** - total charge of your system, e.g. 0
- **--mult** - total multiplicity of your system, e.g. SINGLET
- **--step** - proton shift length, in Å, e.g. 0.05
- **--isotope** - isotope H,D or T. E.g. D
- **--moving_atom_h** - nummber of exchanging H atom in .xyz file. **STARTS FROM 0**. E.g. 21
- **--dest_atom** - number of the atom in .xyz file with which the proton is exchanged. E.g. 14
- **--h_neighbour** - number of the atom in .xyz file with which the proton is connected. E.g. 13
- **--place** - path to calc folder, e.g. `/mop/YOUR_SYSTEM_NAME/calc`.
- **--threads** - number of threads of your processor, e.g. 100

**additional parameters**
- **--other_isotopes** - numbers of other H atoms in molecules should be isotopes **separated by a space**. E.g. 2 3 4
- **--freez** - numbers of other atoms, that should be freezed in geometry optimization **separated by a space**. E.g. 2 3 4
- **--freez_all_axis** - freezing of all axises. 
5. Paste this starting comand in new `YOUR_SYSTEM_NAME.sh` file. This file should be in `proposed_mopac_calc` folder.
Example YOUR_SYSTEM_NAME.sh file:
```
#!/bin/bash
#export RSH_COMMAND="/usr/bir/ssh -x"
export OMP_NUM_THREADS=70
python3 main.py -x opted_h2o.xyz --eps XX.X --charge X --mult XXXXXXX --step X.XX --isotope X \
--moving_atom_h XX --dest_atom XX --h_neighbour XX --place /mop/YOUR_SYSTEM_NAME/calc --threads 100
```
6. Start your calculation by created command, or if you have `slurm` on your machine,
you can start calculations with num threads and time.
```
sbatch -n 100 -t 1:00:00 YOUR_SYSTEM_NAME.sh
```
7. After calculations done, you can create plots and get KIE.
8. For plots, type:
```
python3 estimation.py --place /mop/YOUR_SYSTEM_NAME/calc --plot True
```
9. Determine **Transition State (TS)**
This will create a plot **Heat of Formation (HOF)** from protons shift/step number.
Path to plot `/mop/YOUR_SYSTEM_NAME/calc/HOFS.pdf` and `/mop/YOUR_SYSTEM_NAME/calc/HOFS.png`. 
**.pdf** HOF from (Stem number), but **.png** HOF from shift,Å.
You shold use **.pdf** for determine **TS** point.
10. After you get **TS** point, e.g. 73 from **.pdf**, and your isotope was T, you can type:
```
python3 estimation.py --place /mop/YOUR_SYSTEM_NAME/calc --ts 73 --isotope T
```
and get KIE value. Also, KIE will be saved in `/mop/YOUR_SYSTEM_NAME/calc/kie.txt`
If your isotope was D, then change on _--isotope D_

## Article Data
`data` - contains all data for deuterium
`data_t` - contains all data for tritium

## We have global file structure:
proposed_mopac
- mop
  - substance_name
    - start
    - calc
- saddle
  - substance_name
    - start
    - end
    - saddle
- for_orca
- main.py
- saddle.py
- u_functions.py
- estimation.py
- new_MOPAC_calculator.py
- orca_estimation.py
- orca_calculate_folder.py
- saddle_example.mop
- substance.sh
- proposed_mopac.yml

