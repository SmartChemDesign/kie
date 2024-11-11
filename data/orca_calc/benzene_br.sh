#!/bin/bash
#export RSH_COMMAND="/usr/bir/ssh -x"
#python3 orca_calculate_folder.py
python3 main.py -x opted_h2o.xyz --eps 0 --step 0.025 --mult DOUBLET --charge 0 --isotope D --moving_atom_h 14 -d 15 -n 11 -p /mop/benzene_br/calc