#!/bin/bash
#export RSH_COMMAND="/usr/bir/ssh -x"
#python3 orca_calculate_folder.py
python3 main.py -x opted_h2o.xyz --eps 32.66 --charge 0 --step 0.025 --isotope D --moving_atom_h 20 -d 32 -n 5 -p /mop/pyrq_meoh/calc