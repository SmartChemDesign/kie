#!/bin/bash
#export RSH_COMMAND="/usr/bir/ssh -x"
#python3 orca_calculate_folder.py
python3 main.py -x opted_h2o.xyz --eps 78.4 --charge -1 --step 0.025 --isotope D --moving_atom_h 5 -d 10 -n 4 -p /mop/ch3ch2no2_oh_new/calc