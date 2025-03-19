#!/bin/bash
#export RSH_COMMAND="/usr/bir/ssh -x"
#python3 orca_calculate_folder.py
python3 saddle.py -x opted_h2o.xyz --eps 78.4 --charge -1  --step 0.025  --isotope D --moving_atom_h 2 -d 15 -n 0 -p /mop/ch3no2_ch2clcoo/