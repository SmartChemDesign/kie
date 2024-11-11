#########################################################################
#                                                                       #
#     This module estimates kie from module: "new_MOPAC_calculator"     #
#     You should choose directory, then you will take kie.              #
#                                                                       #
#########################################################################
# based on ase 3.22.1
# necessary inputs
import sys
import csv
import argparse
import numpy as np
import ase
import numpy as np
from u_functions import *
from ase import Atoms
from new_MOPAC_calculator import MYMOPAC
import os
from ase.io import read, write
from ase import Atoms
from ase.units import kcal, mol, Debye
from ase.build import molecule
from ase.calculators.calculator import FileIOCalculator, ReadError, Parameters
from ase.calculators.mopac import MOPAC
#import MOPAC
import pandas as pd
import subprocess
import csv
import matplotlib.pyplot as plt

def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--place', type=str, default="/")
    parser.add_argument('-a', '--all', type=bool, default=False)
    return parser

def find_G(place):
    with open(place, "r") as file:
        data = file.readlines()
    for line in data:
        if "Final Gibbs free energy" in line:
            gibbs = float(line.strip().split("...")[1].replace("Eh", "").strip())*2625.5
            #print(gibbs)
    return gibbs

def check_imag_mods(place):
    sum = 0
    with open(place, "r") as file:
        data = file.readlines()
    for line in data:
        if "***imaginary mode***" in line:
            sum+=1
    return sum

def get_all_sh():
    all_f = os.listdir()
    all_sh = []
    for i in all_f:
        if i[-3:] == ".sh":
            all_sh.append(i)
    return all_sh


if __name__ == '__main__':
    parser = createParser()
    namespace = parser.parse_args(sys.argv[1:])
    cur_dir = find_curent_directory()
    place = f"{cur_dir}{namespace.place}"
    #print(f"{place}/H.csv")

    if namespace.all == True:
        list_all=[]
        open('all_kies.txt', 'w').close()
        shs = get_all_sh()
        for i in shs:
            tmp=[]
            try:
                place = f"{cur_dir}/for_orca/{i[:-3]}"

                start_H = find_G(f"{place}/orca/start_freq_H/start_f_H.out")
                ts_H = find_G(f"{place}/orca/freq_H/f_H.out")
                start_D = find_G(f"{place}/orca/start_freq_D/start_f_D.out")
                ts_D = find_G(f"{place}/orca/freq_D/f_D.out")

                freqs_ts = [f"{place}/orca/freq_H/f_H.out",
                            f"{place}/orca/freq_D/f_D.out"]

                freqs_start = [f"{place}/orca/start_freq_D/start_f_D.out",
                            f"{place}/orca/start_freq_H/start_f_H.out"]

                delta_H_kJ = (ts_H - start_H)
                delta_D_kJ = (ts_D - start_D)
                #print(ts_H, ts_D)
                #print(start_H, start_D)
                #print(f"H_barier = {delta_H_kJ}")
                #print(f"D_barier = {delta_D_kJ}")
                kie = np.exp(-(delta_H_kJ-delta_D_kJ)*1000/8.314/298)


                #with open ("all_kies.txt", "a") as file:
                #    file.write(i[:-3])
                #    file.write("\n")
                tmp.append(i[:-3])

                print(kie)
                #print(place)
                #with open (f"{place}/orca/kie.txt", "w") as file:
                #    file.write(str(kie))
                #print("kie_in_folder_done")
                tmp.append(str(kie))

                
                print("bariers_done")
                for f in freqs_ts:
                    flag = True
                    n_im =check_imag_mods(f)
                    if n_im !=1 :
                        print(f"ts {n_im} im mode")
                        #with open ("all_kies.txt", "a") as file:
                        #    file.write(f"ts {n_im} im mode")
                        #    file.write("\n")
                        flag = False
                    if flag:
                        print("ts_freqs has 1 im mode, OK")
                        #with open ("all_kies.txt", "a") as file:
                        #    file.write("ts_freqs has 1 im mode, OK")
                        #    file.write("\n")
                        #tmp.append("1")
                    tmp.append(n_im)
                print("fr_ts_done")
                for f in freqs_start:
                    flag = True
                    n_im = check_imag_mods(f)
                    if check_imag_mods(f) != 0:
                        print(f"start {n_im} im mode")
                        #with open ("all_kies.txt", "a") as file:
                        #    file.write(f"start {n_im} im mode")
                        #    file.write("\n")
                        flag = False
                        #tmp.append(n_im)
                    if flag:
                        print("start has not im mods, OK")
                    tmp.append(n_im)
                print("fr_st_done")
                #with open ("all_kies.txt", "a") as file:
                #    file.write(str(kie))
                #    file.write("\n")
            except:
                print(f"{i} is broken")
            if len(tmp)>3:
                list_all.append(tmp)
        with open("all_kies.csv", 'w') as myfile:
                    wr = csv.writer(myfile)
                    wr.writerows(list_all)

         
         
         
         
         
         
         
         
         
            




    start_H = find_G(f"{place}/orca/start_freq_H/start_f_H.out")
    ts_H = find_G(f"{place}/orca/freq_H/f_H.out")
    start_D = find_G(f"{place}/orca/start_freq_D/start_f_D.out")
    ts_D = find_G(f"{place}/orca/freq_D/f_D.out")

    freqs_ts = [f"{place}/orca/freq_H/f_H.out",
                f"{place}/orca/freq_D/f_D.out"]

    freqs_start = [f"{place}/orca/start_freq_D/start_f_D.out",
                   f"{place}/orca/start_freq_H/start_f_H.out"]

    delta_H_kJ = (ts_H - start_H)
    delta_D_kJ = (ts_D - start_D)
    #print(ts_H, ts_D)
    #print(start_H, start_D)
    print(f"H_barier = {delta_H_kJ}")
    print(f"D_barier = {delta_D_kJ}")
    kie = np.exp(-(delta_H_kJ-delta_D_kJ)*1000/8.314/298)

    print(kie)
    with open (f"{place}/orca/kie.txt", "w") as file:
        file.write(str(kie))

    for f in freqs_ts:
        flag = True
        n_im =check_imag_mods(f)
        if n_im !=1 :
            print(f"{n_im} im mode")
            flag = False
        if flag:
            print("ts_freqs has 1 im mode, OK")
    for f in freqs_start:
        flag = True
        n_im = check_imag_mods(f)
        if check_imag_mods(f) != 0:
            print(f"{n_im} im mode")
            flag = False
        if flag:
            print("start has not im mods, OK")




