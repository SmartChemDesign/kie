#########################################################################
#                                                                       #
#     This module contains all functions for moving atoms, optimizing   #
#     geometry, rotate molecule, creating files and other.              #
#                                                                       #
#########################################################################
# based on ase 3.22.1
# necessary inputs

import numpy as np
import os
from ase.io import read, write
from ase import Atoms
from ase.calculators.calculator import FileIOCalculator, ReadError,Calculator
from ase.units import kcal, mol, Debye
from ase.build import molecule
from ase.calculators.calculator import FileIOCalculator, ReadError, Parameters
from ase.calculators.mopac import MOPAC
from new_MOPAC_calculator import MYMOPAC
#import MOPAC
import shutil
import subprocess
import csv



def opt_geom_saddle(filename_inp,filename_xyz,
             base_task_opt = "BONDS AUX DUMP=30M ITRY=15000 T=14D",
             charge=0, mult="SINGLET", eps="78.4", place = "/", freez =[]):
    print(filename_inp[:-4])
    xyz = read(filename_inp)
    atoms = Atoms(xyz)
    if eps ==0:
        task_line = base_task_opt + f" CHARGE={charge} {mult}"
    else:
        task_line = base_task_opt + f" CHARGE={charge} {mult} eps={eps}"
    atoms.calc = MYMOPAC(label=filename_inp[:-4],
                         task=task_line)
    atoms.calc.write_input(atoms)
    if len(freez)>0:
        for i in freez:
            add_flag_zero(filename_inp[:-4], i)


    atoms.calc.calculate_not_create_input(atoms)
    
    with open(atoms.calc.label + '.out') as fd:
        lines = fd.readlines()
    res = atoms.calc.read_atoms_from_out_opt(lines)
    res.write(filename_xyz)
    with open(filename_xyz,"r") as file1:
        xyz_data1=file1.readlines()

def str_array_to_int_array(str_arr):
    print(f"str_array = {str_arr}")
    r_arr = []
    for i in str_arr:
        r_arr.append(int(i))
    print(r_arr)
    return r_arr

def find_move_vectore(first_geom, moving_atom_h, dest_atom, step, h_neighbour):
    mol_xyz_1 = read(f"{first_geom}")
    mol_xyz = rotate_xyz_to_x(mol_xyz_1, moving_atom_h, dest_atom)
    distance_h_dest = mol_xyz.get_distance(moving_atom_h, dest_atom, vector=False)
    distance_h_neigh = mol_xyz.get_distance(moving_atom_h, h_neighbour, vector=False)
    distance_vector = np.array(mol_xyz.get_distance(moving_atom_h, dest_atom, vector=True))
    step_vector = distance_vector / (distance_h_dest / step)
    return mol_xyz, step_vector, distance_h_dest, distance_h_neigh

def rotate_xyz_to_x(xyz, h_atom, dest_atom):
    """
    move forward x axis
    :param xyz:
    :return:
    """
    distance_vector = np.array(xyz.get_distance(h_atom, dest_atom, vector=True))
    xyz.rotate(distance_vector, (1,0,0))
    return xyz


def write_zpe(zpe,isotope, place):
    if isotope=="H":
        with open(f"{place}/H.csv" , "a", newline="") as file:
            writer = csv.writer(file, delimiter=' ')
            writer.writerow([zpe])
    if isotope=="D":
        with open(f"{place}/D.csv" , "a", newline="") as file:
            writer = csv.writer(file, delimiter=' ')
            writer.writerow([zpe])
    if isotope=="T":
        with open(f"{place}/T.csv" , "a", newline="") as file:
            writer = csv.writer(file, delimiter=' ')
            writer.writerow([zpe])

def create_files(place):
    fp = open(place + "/" + "H.csv", 'w')
    fp.close()
    fp = open(place + "/" + "D.csv", 'w')
    fp.close()

    fp = open(place + "/" + "H_hof.csv", 'w')
    fp.close()
    fp = open(place + "/" + "D_hof.csv", 'w')
    fp.close()
def create_folder_my(folder):
    try:
        os.mkdir(folder)
    except FileExistsError:
        shutil.rmtree(folder)
        os.mkdir(folder)

def write_HOF(zpe,isotope, place):
    if isotope=="H":
        with open(place + "/" + "H_hof.csv" , "a", newline="") as file:
            writer = csv.writer(file, delimiter=' ')
            writer.writerow([zpe])
    if isotope=="D":
        with open(place + "/" + "D_hof.csv" , "a", newline="") as file:
            writer = csv.writer(file, delimiter=' ')
            writer.writerow([zpe])
    if isotope=="T":
        with open(place + "/" + "T_hof.csv" , "a", newline="") as file:
            writer = csv.writer(file, delimiter=' ')
            writer.writerow([zpe])

def write_GIBBS(zpe,isotope, place):
    if isotope=="H":
        with open(place + "/" + "H_gibbs.csv" , "a", newline="") as file:
            writer = csv.writer(file, delimiter=' ')
            writer.writerow([zpe])
    if isotope=="D":
        with open(place + "/" + "D_gibbs.csv" , "a", newline="") as file:
            writer = csv.writer(file, delimiter=' ')
            writer.writerow([zpe])
    if isotope=="T":
        with open(place + "/" + "T_gibbs.csv" , "a", newline="") as file:
            writer = csv.writer(file, delimiter=' ')
            writer.writerow([zpe])



def find_first_thermo(opted_geom,thermo_name,
                      base_task_opt = "AUX LARGE THERMO(298,298) LET ISOTOPE BONDS DUMP=30M ITRY=15000 T=14D",
                      charge=0, mult="SINGLET", eps="78.4", place = "/", threads = 1):
    xyz = read(opted_geom)
    atoms = Atoms(xyz)
    if eps == 0:
        task_line = base_task_opt + f" CHARGE={charge} {mult} THREADS={threads}"
    else:
        task_line = base_task_opt + f" CHARGE={charge} {mult} eps={eps} THREADS={threads}"
    atoms.calc = MYMOPAC(label=thermo_name,
                         task=task_line)
    atoms.calc.calculate(atoms)

    atoms.calc.read_results()
    write_zpe(atoms.calc.get_zero_point_energy(),"H", place=place)
    write_HOF(atoms.calc.get_my_hof(), "H", place=place)
    write_GIBBS(find_gibbs(f"{thermo_name}.out"), "H", place=place)
    return 0

def find_gibbs(out):
    with open(out, "r") as fd:
        lines = fd.readlines()
    for i,line in enumerate(lines):
        if "TOT." in line:
            tmp = lines[i].strip().split()
    h = float(tmp[2].strip())
    s = float(tmp[4].strip())
    print(h,s)
    return (h-298*s)*4.18


def find_curent_directory():
    if "\\" in os.getcwd():
        return os.getcwd()
    else:
        return os.getcwd()


def change_isotopes_input(mopac,positions,isotope,
                          base_task_opt = "AUX LARGE THERMO(298,298) RESTART BONDS DUMP=30M ITRY=15000 T=14D",
                          charge=0, mult="SINGLET", eps="78.4", threads=1): # с нуля
    with open(mopac + '.mop', "r") as fd:
        lines = fd.readlines()
        #lines[2]="AUX LARGE CHARGE=-1 SINGLET eps=81 THERMO(298,298) RESTART BONDS DUMP=30M ITRY=15000 T=14D  THREADS=40\n"
        if eps == 0:
            lines[2] = base_task_opt + f" CHARGE={charge} {mult} THREADS={threads}\n"
        else:
            lines[2] = base_task_opt + f" CHARGE={charge} {mult} eps={eps} THREADS={threads}\n"
        for atom in positions:
            if "H" in lines[3+atom]:
                #print(1)
                lines[3+atom]=lines[4+atom].replace("H", isotope)
            elif "D" in lines[4+atom]:
                lines[3+atom]=lines[4+atom].replace("D", isotope)
            elif "T" in lines[4+atom]:
                lines[3+atom]=lines[4+atom].replace("T", isotope)
        #print(lines)
    os.remove(f"{mopac}.mop")
    with open(mopac + '.mop', "w") as fd1:
        fd1.writelines(lines)

def add_flag_minus(mopac,number_atom,axis): #c 0
    with open(mopac + '.mop', "r") as fd:
        lines = fd.readlines()
        tmp=lines[3+number_atom].split()
        if axis == "X":
            tmp[2]="-1"
        elif axis == "Y":
            tmp[4]="-1"
        elif axis == "Z":
            tmp[6]="-1"
        line_new=f"{ tmp[0]}  {tmp[1]}  {tmp[2]}  {tmp[3]}  {tmp[4]}  {tmp[5]}  {tmp[6]}  \n"
        lines[4 + number_atom]=line_new
        # print(lines)
    os.remove(f"{mopac}.mop")
    with open(mopac + '.mop', "w") as fd1:
        fd1.writelines(lines)


def add_flag_zero(mopac,number_atom): #c 0
    with open(mopac + '.mop', "r") as fd:
        lines = fd.readlines()
        tmp=lines[3+number_atom].split()
        line_new=f"{ tmp[0]}  {tmp[1]}  {0}  {tmp[3]}  {0}  {tmp[5]}  {0}  \n"
        lines[3 + number_atom]=line_new
        # print(lines)
    os.remove(f"{mopac}.mop")
    with open(mopac + '.mop', "w") as fd1:
        fd1.writelines(lines)

def add_flag_zero_x(mopac,number_atom): #c 0
    with open(mopac + '.mop', "r") as fd:
        lines = fd.readlines()
        tmp=lines[3+number_atom].split()
        line_new=f"{ tmp[0]}  {tmp[1]}  {0}  {tmp[3]}  {tmp[4]}  {tmp[5]}  {tmp[6]}  \n"
        lines[3 + number_atom]=line_new
        # print(lines)
    os.remove(f"{mopac}.mop")
    with open(mopac + '.mop', "w") as fd1:
        fd1.writelines(lines)

def add_tv_to_mop(mopac,tv):
    with open(mopac + '.mop', "r") as fd:
        lines = fd.readlines()
        for l in tv:
            lines.append(l)
    os.remove(f"{mopac}.mop")
    with open(mopac + '.mop', "w") as fd1:
        fd1.writelines(lines)

def find_isotopes_thermo(opted_geom,thermo_name,positions,isotope,
                         base_task_opt = "AUX LARGE THERMO(298,298) RESTART BONDS DUMP=30M ITRY=15000 T=14D",
                         charge=0, mult="SINGLET", eps="78.4", place="/", threads=1):
    xyz = read(opted_geom)

    atoms = Atoms(xyz)
    if eps == 0:
        task_line = base_task_opt + f" CHARGE={charge} {mult} THREADS={threads}"
    else:
        task_line = base_task_opt + f" CHARGE={charge} {mult} eps={eps} THREADS={threads}"
    atoms.calc = MYMOPAC(label=thermo_name,
                         task=task_line)

    atoms.calc.change_isotopes_input(positions,isotope)
    atoms.calc.calculate_not_create_input()



    atoms.calc.read_results()
    write_zpe(atoms.calc.get_zero_point_energy(),isotope, place = place)
    write_HOF(atoms.calc.get_my_hof(), isotope, place = place)
    write_GIBBS(find_gibbs(f"{thermo_name}.out"), isotope, place=place)


def opt_geom(filename_inp,filename_xyz,xyz_trjry, h_num, h_neighbour, dest_num,
             base_task_opt = "BONDS AUX DUMP=30M ITRY=15000 T=14D",
             charge=0, mult="SINGLET", eps="78.4", place = "/", freez =[], freez_all_axis=False, threads=1,
             trj = True, preopt=False):
    print(filename_inp[:-4])
    xyz = read(filename_inp)
    atoms = Atoms(xyz)
    if eps ==0:
        task_line = base_task_opt + f" CHARGE={charge} {mult} THREADS={threads}"
    else:
        task_line = base_task_opt + f" CHARGE={charge} {mult} eps={eps} THREADS={threads}"
    atoms.calc = MYMOPAC(label=filename_inp[:-4],
                         task=task_line)
    atoms.calc.write_input(atoms)
    if freez_all_axis:
        add_flag_zero(filename_inp[:-4], h_num)
    else:
        add_flag_zero_x(filename_inp[:-4],h_num)
    if not preopt:
        add_flag_zero(filename_inp[:-4], h_neighbour)
        add_flag_zero(filename_inp[:-4], dest_num)
        if len(freez)>0:
            for i in freez:
                add_flag_zero(filename_inp[:-4], i)

    atoms.calc.calculate_not_create_input(atoms)
    
    with open(atoms.calc.label + '.out') as fd:
        lines = fd.readlines()
    res = atoms.calc.read_atoms_from_out_opt(lines)
    res.write(filename_xyz)
    if trj:
        with open(filename_xyz,"r") as file1:
            xyz_data1=file1.readlines()
        with open(xyz_trjry,"a") as file2:
            file2.writelines(xyz_data1)

def move_atom(opted_geom, atom, vector, first_geom):
    xyz = read(opted_geom)
    xyz.positions[atom] += vector
    #atoms = Atoms(xyz,cell=cell_xyz.cell)
    atoms = Atoms(xyz)
    atoms.write(f"{first_geom[:-4]}.xyz")


#1)ОПТ в circle
#2) Прочесть из аута геометрию и сохранить
#3) создать инпут для термо и опсчитать его
#4) заменить атомы ерех change_isotopes
#5) Создать инпут с restart и посчитать 1 изотоп. Сохранить ZPE в список
#6) создать инпут под второй изотоп и посчитать термо м


#first_geom_cell="1008775_cell.cif"
#first_geom="opted_h2o.xyz"
#opted_geom="opted_h2o_x2.xyz"
#thermo_name="opted_h2o_thermo"
#xyz_trjry="p.xyz"
#cell_xyz = read(first_geom_cell)# новую закинуть
