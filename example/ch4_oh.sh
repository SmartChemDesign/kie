#!/bin/bash
#export RSH_COMMAND="/usr/bir/ssh -x"
#python3 orca_calculate_folder.py
python3 main.py -x opted_h2o.xyz --eps 0 --charge 0 --mult DOUBLET --isotope D --moving_atom_h 1 -d 5 -n 0 -p /mop/ch4_oh/calc
