#!/bin/bash
#export RSH_COMMAND="/usr/bir/ssh -x"
#python3 orca_calculate_folder.py
python3 main.py -x opted_h2o.xyz --eps 78.4 --mult SINGLET --charge -1 --step 0.025 --isotope T --moving_atom_h 15 -d 41 -n 13 -t 100 -p /mop/2-phenylacetophenone/calc
