#!/bin/bash
#export RSH_COMMAND="/usr/bir/ssh -x"
#python3 orca_calculate_folder.py
python3 main.py -x opted_h2o.xyz --eps 78.4 --charge -1 --step 0.025 --isotope D --moving_atom_h 10 -d 20 -n 2 -p /mop/nitrocyclohexane_oh_new/calc