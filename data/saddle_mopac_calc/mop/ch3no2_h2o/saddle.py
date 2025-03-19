#########################################################################
#                                                                       #
#     This module calculates kie proton transfer, it takes              #
#     start xyz geom, num scatering atom H, num host atom,              #
#     eps dielectric constant and moving shag. Atoms from 0!            #
#                                                                       #
#########################################################################
# based on ase 3.22.1
# necessary inputs
import sys
import shutil
from time import time
import argparse
import numpy as np
import ase
from u_functions import *
from ase import Atoms
from new_MOPAC_calculator import MYMOPAC, find_curent_directory
import os
from ase.io import read, write
from ase import Atoms
from ase.units import kcal, mol, Debye
from ase.build import molecule
from ase.calculators.calculator import FileIOCalculator, ReadError, Parameters
from ase.calculators.mopac import MOPAC
#import MOPAC
import subprocess
import csv
import pandas as pd

def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--xyz', type=str)
    parser.add_argument('-m','--moving_atom_h', type=int, default=0)
    parser.add_argument('-d', '--dest_atom', type=int, default=1) # H moving to
    parser.add_argument('-e', '--eps', type=float, default=78.4)
    parser.add_argument('-s', '--step', type=float, default=0.01)
    parser.add_argument('--mult', type=str, default="SINGLET")
    parser.add_argument('--charge', type=int, default=0)
    parser.add_argument('--isotope', type=str, default="D")
    parser.add_argument('--other_isotopes', nargs='+', help='<Required> Set flag', required=False)
    parser.add_argument('-n','--h_neighbour', type=int, default=2)# H moving from
    parser.add_argument('-p', '--place', type=str, default="/")  # H moving from
    parser.add_argument('--freez', nargs='+', help='<Required> Set flag', required=False, default="")
    parser.add_argument('--freez_all_axis', type=bool, default=False)
    return parser


if __name__ == '__main__':
    parser = createParser()
    namespace = parser.parse_args(sys.argv[1:])
    cur_dir = find_curent_directory()
    freez = str_array_to_int_array(namespace.freez)
    if namespace.other_isotopes != None:
        other_isotopes = str_array_to_int_array(namespace.other_isotopes)
        if namespace.moving_atom_h not in other_isotopes:
            other_isotopes.append(namespace.moving_atom_h)
        isotope_positions = other_isotopes
    else:
        isotope_positions = [namespace.moving_atom_h]
    print(f"isotope_positions {isotope_positions}")
    #print(f"other_isotopes {other_isotopes}")

    place = f"{cur_dir}{namespace.place}"
    create_files(place)

    t1 = time()
    #-saddle
    # -----------------------------opt start geom-----------------------------------
    start_geom = "start.xyz"
    start_geom_opted = "start_opt.xyz"

    mol_start_geom = read(f"{place}/start/{start_geom}")
    filename_start_geom = place + "/start/" + start_geom
    filename_start_geom_opted = place + "/start/" + start_geom_opted

    opt_geom_saddle(filename_inp=filename_start_geom, filename_xyz=filename_start_geom_opted,
             charge=namespace.charge, mult=namespace.mult, eps=namespace.eps, freez=freez)  # первая геома и геометрия опта
    # -----------------------------opt end geom-----------------------------------
    end_geom = "end.xyz"
    end_geom_opted = "end_opt.xyz"

    mol_end_geom = read(f"{place}/end/{end_geom}")
    filename_end_geom = place + "/end/" + end_geom
    filename_end_geom_opted = place + "/end/" + end_geom_opted

    opt_geom_saddle(filename_inp=filename_end_geom, filename_xyz=filename_end_geom_opted,
                    charge=namespace.charge, mult=namespace.mult, eps=namespace.eps,
                    freez=freez)  # первая геома и геометрия опта
    # -----------------------------saddle-----------------------------------
    create_folder_my(f"{place}/saddle")
    shutil.copy(f"{place}/saddle_example.mop", f"{place}/saddle/saddle.mop")

    with open(f"{filename_start_geom[:-4]}.mop", "r", encoding="UTF-8") as p:
        par = p.readlines()
        good_par = "SADDLE" + par[3]

    with open(filename_start_geom_opted, "r", encoding="UTF-8") as st:
        xyz_start = st.readlines()
    with open(filename_end_geom_opted, "r", encoding="UTF-8") as en:
        xyz_end = en.readlines()

    with open(f"{place}/saddle/saddle.mop", "r", encoding="UTF-8") as bad_input:
        b_inp = bad_input.readlines()

    with open(f"{place}/saddle/saddle.mop", "w", encoding="UTF-8") as good_input:
        for line in b_inp[:-1]:
            good_input.write(line)

        good_input.write(good_par)

        for line in xyz_start:
            good_input.write(line)

        good_input.write("\n")

        for line in xyz_end:
            good_input.write(line)

    command = f"mopac {place}/saddle/saddle.mop"
    proc = subprocess.Popen(command, shell=True, cwd=f"{dir}/orca/neb/")
    # -----------------------------take xyz from saddle-----------------------------------
    with open(f"{place}/saddle/saddle.mop", "w", encoding="UTF-8") as good_saddle:
        saddle_xyz = good_saddle.readlines()
    flag = 0
    res = []
    for i,line in enumerate(saddle_xyz):
        if "CARTESIAN COORDINATES" in line:
            flag+=1
        if flag == 2:
            start = i
    for j in saddle_xyz[i+2:]:
        while j != " ":
            j_g = j.strip().split()[1:]
            j_good = ""
            for k in j_g:
                j_good+= f" {k}"
        res.append(j_good)
    n_xyz = len(j_good)
    with open(f"{place}/saddle/saddle_opt.xyz", "w", encoding="UTF-8") as saddle_opted:
        saddle_opted.write(f"{n_xyz}\n")
        saddle_opted.write("\n")
        for a in res:
            saddle_opted.write(f"{a}\n")

    # -----------------------------find thermo-----------------------------------
    find_first_thermo(opted_geom=f"{place}/saddle/saddle_opt.xyz", thermo_name=f"H_thermo",
                      charge=namespace.charge, mult=namespace.mult,
                      eps=namespace.eps, place=place)
    find_isotopes_thermo(opted_geom=f"{place}/saddle/saddle_opt.xyz", thermo_name=f"H_thermo",
                         positions=isotope_positions,
                         isotope=namespace.isotope, charge=namespace.charge,
                         mult=namespace.mult, eps=namespace.eps, place=place)

    t2 = time()
    with open(f"{place}/time.txt", "w") as file:
        file.write(f"{place} {t2 - t1}\n")

    df_H = pd.read_csv(f"{place}/H.csv", header=None)
    df_H_column_data = df_H.iloc[:, 0]
    df_H_column_list = df_H_column_data.to_list()

    df_D = pd.read_csv(f"{place}/{namespace.isotope}.csv", header=None)
    df_D_column_data = df_D.iloc[:, 0]
    df_D_column_list = df_D_column_data.to_list()

    start_H = float(df_H_column_list[0])
    ts_H = float(df_H_column_list[1])
    start_D = float(df_D_column_list[0])
    ts_D = float(df_D_column_list[1])

    delta_H_kJ = (ts_H - start_H) * 4.184
    delta_D_kJ = (ts_D - start_D) * 4.184

    kie = np.exp(-(delta_H_kJ - delta_D_kJ) * 1000 / 8.314 / 298)

    print(f"kie = {kie}")
    with open(f"{place}/kie.txt", "w") as file:
        file.write(f"kie = {str(kie)}\n")

