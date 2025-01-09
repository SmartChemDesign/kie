#********************************************************
#
# this code takes XYZ files from directories,
# creates orca inputs and calculates them
#
#********************************************************
# define imports
# ase 3.22.1
from ase.calculators.orca import ORCA
from ase.io import read, write
from ase.atoms import Atoms
from ase.optimize import BFGS
import os
import subprocess
import re
import shutil
import numpy as np
import datetime

from time import time
from ase.units import Hartree, Bohr
from ase.io.orca import write_orca
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError

#thre = 10

def create_folder_my(folder):
    try:
        os.mkdir(folder)
    except FileExistsError:
        shutil.rmtree(folder)
        os.mkdir(folder)





# define my Orca class, because ase has biiiiiiiiig troubles

class My_ORCA(ORCA):
    def __init__(self, restart=None,
                 ignore_bad_restart_file=FileIOCalculator._deprecated,
                 label='ab', atoms=None,  **kwargs):

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)
        self.pcpot = None

    def calculate_not_create_input(self, path_to_orca):
        command = f"{path_to_orca} {self.label}.inp > {self.label}.out"
        #module_load = r"module load ompi/4.1.0"
        #mod = subprocess.Popen(module_load, shell=True, cwd=self.directory)
        #mod.wait()
        print(command)
        try:
            proc = subprocess.Popen(command, shell=True, cwd=self.directory)
        except OSError as err:
            # Actually this may never happen with shell=True, since
            # probably the shell launches successfully.  But we soon want
            # to allow calling the subprocess directly, and then this
            # distinction (failed to launch vs failed to run) is useful.
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err
        errorcode = proc.wait()

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters
        p.write(self.label + '.ase')
        p['label'] = self.label
        if self.pcpot:  # also write point charge file and add things to input
            p['pcpot'] = self.pcpot

        write_orca(atoms, **p)

        with open(self.label + '.inp', "r") as bad_input:
            inp = bad_input.readlines()
            inp[0] = inp[0].replace("engrad", "")

        with open(self.label + '.inp', "w") as good_input:
            for line in inp:
                good_input.write(line)


# define directories
def orca_calculate_xyz(file, label, lom_path_orca,
        change_isotope = False,
        mopac_isotope = "none.mop",
        sh_file = "none.sh",
        isotope = "2",
        charge=-1,
        mult=1,
        eps = 78.4,
        maxiter=2000,
        orcasimpleinput='opt B3LYP def2-TZVP CPCM D3',
        orcablocks='%pal nprocs 102 end'):


    orc_eps = f"%cpcm epsilon {eps} end"
    orcablocks = orcablocks + "\n" + orc_eps

    calc_orca = My_ORCA(label=label,
        charge=charge,
        mult=mult,
        maxiter=maxiter,
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks)
    mo = read(file)
    at = Atoms(mo)
    calc_orca.write_input(at)
    if change_isotope:
        isotope_num = find_isotope_number_from_sh(sh_file)
        #isotope_num = find_isotope_number_from_mopac(mopac_isotope)
        orca_change_isotope(label+".inp", isotope_num, isotope)

    print("start_calc")
    calc_orca.calculate_not_create_input(lom_path_orca)


def check_norm_end(file):
    with open(file, "r", encoding="UTF-8") as out_res:
        res = out_res.readlines()
    for i in res:
        if "NORMALLY" in i:
            return True
    return False


def find_isotope_number_from_sh(sh_file):
    with open(sh_file, "r", encoding="UTF-8") as bad_input:
        inp = bad_input.readlines()

    for i, line in enumerate(inp):
        if "main.py" in line:
            st = i
    data_inp = inp[i].strip().split(" ")

    for k,j in enumerate(data_inp):
        if "--moving_atom_h" in j.lower():
            d = data_inp[k+1]
    print(int(d))
    return(int(d))


def find_isotope_number_from_mopac(mopac_isotope):
    with open(mopac_isotope, "r", encoding="UTF-8") as bad_input:
        inp = bad_input.readlines()
    for i in inp:
        if ("move_atom" in i) and("[" in i):
            d = i.split(",")
    print(int(d[1]))
    return(int(d[1]))


def orca_change_isotope(file, position, isotope):
    with open(file, "r", encoding="UTF-8") as bad_input:
        start = 0
        inp = bad_input.readlines()
    for i,j in enumerate(inp):
        if "*xyz" in j:
            start = i+1

    print(int(start+position))
    print( inp[int(start+position)])
    print(inp[int(start+position)][0])

    inp[int(start + position)] = inp[int(start + position)].strip() + f" M = {isotope}" + "\n"
    with open(file, "w") as good_input:
        for line in inp:
            good_input.write(line)


    return 0
def orca_read_mult_from_mopac_input(mopac_file):
    with open(mopac_file, "r", encoding="UTF-8") as bad_input:
        inp = bad_input.readlines()
    #print(inp)
    for line in inp[3:]:
        if "CHARGE" in line:
            #print(line)
            charge_l = line.split("=")
            for i,j in enumerate(charge_l):
                if "CHARGE" in j:
                    charge = charge_l[i+1].split(" ")[0]
        if "SINGLET" in line.upper():
            mult = 1
        if "DOUBLET" in line.upper():
            mult = 2
        if "EPS" in line.upper():
            #print(line)
            eps_l = line.split("=")
            #print(eps_l)
            for k, l in enumerate(eps_l):
                if "EPS" in l.upper():
                    #print(l)
                    eps = eps_l[k+1].split(" ")[0]
                    #print(eps)


    #print(charge,mult)
    return int(charge), int(mult), float(eps)


def orca_read_mult_from_sh(sh_file):
    with open(sh_file, "r", encoding="UTF-8") as bad_input:
        inp = bad_input.readlines()
    #print(inp)

    eps = 78.4
    charge = 0
    mult = "SINGLET"
    for i, line in enumerate(inp):
        if "main.py" in line:
            st = i
    data_inp  = inp[i].strip().split(" ")
    for k,j in enumerate(data_inp):
        if j.lower() == "--eps":
            eps = data_inp[k+1]
        if j.lower() == "--charge":
            charge = data_inp[k+1]
        if j.lower() == "--mult":
            mult = data_inp[k+1]
    mult_dict = {"SINGLET":1, "DOUBLET":2}


    #print(charge,mult)
    return int(charge), mult_dict[str(mult)], float(eps)




def take_xyz_from_dir(dir, flag):
    files = os.listdir(dir)
    for i in files:
        if flag in i:
            return f"{dir}/{i}"

#***************Calculating***************
def calculate_dir(dir, lom_path_orca, sh_file):
# Calculate start_geometry
    first_xyz = take_xyz_from_dir(f"{dir}/orca/orca_start", ".xyz")
    start_charge, start_mult, start_eps = orca_read_mult_from_sh(f"{sh_file}")
    orca_calculate_xyz(first_xyz,lom_path_orca=lom_path_orca, label = f"{dir}/orca/orca_start/start", charge=start_charge,mult=start_mult)
# Read out xyz from start
    create_folder_my(f"{dir}/orca/neb")
    shutil.copy(f"{dir}/orca/orca_start/start.xyz", f"{dir}/orca/neb/start.xyz")
    shutil.copy(f"{dir}/orca/orca_start/start.out", f"{dir}/orca/neb/start.out")

# read cpcm
# Calculate end_geometry
    end_xyz = take_xyz_from_dir(f"{dir}/orca/orca_end", ".xyz")
    orca_calculate_xyz(end_xyz,lom_path_orca=lom_path_orca, label=f"{dir}/orca/orca_end/end", charge=start_charge, mult=start_mult)
    shutil.copy(f"{dir}/orca/orca_end/end.xyz", f"{dir}/orca/neb/end.xyz")
    shutil.copy(f"{dir}/orca/orca_end/end.out", f"{dir}/orca/neb/end.out")


# Create NEB file
    shutil.copy(f"neb.inp", f"{dir}/orca/neb/neb.inp")
# Add mult, charge and eps
    with open(f"{dir}/orca/neb/neb.inp", "r") as bad_neb:
        neb_data = bad_neb.readlines()
        for i,j in enumerate(neb_data):
            if "CHARGE" in j:
                neb_data[i] = neb_data[i].replace("CHARGE", str(start_charge))
            if "MULT" in j:
                neb_data[i] = neb_data[i].replace("MULT", str(start_mult))
            if "EPS" in j:
                neb_data[i] = neb_data[i].replace("EPS", str(start_eps))


    with open(f"{dir}/orca/neb/neb.inp", "w") as good_neb:
        for line in neb_data:
            good_neb.write(line)

# Run NEB
    command = f"{lom_path_orca} neb.inp > neb.out"
    print(check_norm_end(f"{dir}/orca/neb/end.out") and check_norm_end(f"{dir}/orca/neb/start.out"))
    if check_norm_end(f"{dir}/orca/neb/end.out") and check_norm_end(f"{dir}/orca/neb/start.out"):
        print("running neb")
        try:
            proc = subprocess.Popen(command, shell=True, cwd=f"{dir}/orca/neb/")
        except:
            print("neb error")
        errorcode = proc.wait()

# Create freq_folders thermo_ts
    create_folder_my(f"{dir}/orca/freq_H")
    create_folder_my(f"{dir}/orca/freq_D")
    shutil.copy(f"{dir}/orca/neb/neb.xyz", f"{dir}/orca/freq_H/neb.xyz")
    shutil.copy(f"{dir}/orca/neb/neb.xyz", f"{dir}/orca/freq_D/neb.xyz")
# Run thermo freq_H
    freq_H_xyz = take_xyz_from_dir(f"{dir}/orca/freq_H", ".xyz")
    orca_calculate_xyz(freq_H_xyz, lom_path_orca=lom_path_orca, label = f"{dir}/orca/freq_H/f_H",
                       charge=start_charge,mult=start_mult, orcasimpleinput='NUMFREQ B3LYP def2-TZVP CPCM D3')

# Run thermo freq_D
    freq_D_xyz = take_xyz_from_dir(f"{dir}/orca/freq_D", ".xyz")
    orca_calculate_xyz(freq_D_xyz, lom_path_orca=lom_path_orca,
                       label = f"{dir}/orca/freq_D/f_D", charge=start_charge,mult=start_mult,
                       orcasimpleinput='NUMFREQ B3LYP def2-TZVP CPCM D3', isotope="2",
                       change_isotope=True, sh_file= f"{sh_file}")


# Create freq_folders thermo_ts
    create_folder_my(f"{dir}/orca/start_freq_H")
    create_folder_my(f"{dir}/orca/start_freq_D")
    shutil.copy(f"{dir}/orca/orca_start/start.xyz", f"{dir}/orca/start_freq_H/neb.xyz")
    shutil.copy(f"{dir}/orca/orca_start/start.xyz", f"{dir}/orca/start_freq_D/neb.xyz")

# Run thermo start_freq_H
    start_freq_H_xyz = take_xyz_from_dir(f"{dir}/orca/start_freq_H", ".xyz")
    orca_calculate_xyz(start_freq_H_xyz, lom_path_orca=lom_path_orca, label = f"{dir}/orca/start_freq_H/start_f_H",
                       charge=start_charge,mult=start_mult, orcasimpleinput='NUMFREQ B3LYP def2-TZVP CPCM D3')

# Run thermo start_freq_D
    start_freq_D_xyz = take_xyz_from_dir(f"{dir}/orca/start_freq_D", ".xyz")
    orca_calculate_xyz(start_freq_D_xyz, lom_path_orca=lom_path_orca,
                       label = f"{dir}/orca/start_freq_D/start_f_D", charge=start_charge,mult=start_mult,
                       orcasimpleinput='NUMFREQ B3LYP def2-TZVP CPCM D3', isotope="2",
                       change_isotope=True, sh_file= f"{sh_file}")






# all files
# запуск из папки со скриптами

lom_path_orca = '/home/rudenko/orca/orca'
#lom_path_orca = '/opt/orca/orca'

home_path = os.getcwd()

path_to_orca_dirs = home_path + "/for_orca"
orca_dirs = os.listdir(path_to_orca_dirs)

# Calculate circle

with open ("times.txt", "w") as file:
    file.write(f"{datetime.datetime.now()}\n")

for i in orca_dirs:
    sh_file = f"{i}.sh"
    #calculate_dir(f"{path_to_orca_dirs}/{i}", lom_path_orca, sh_file=sh_file)
    t1 = time()
    try:
        calculate_dir(f"{path_to_orca_dirs}/{i}", lom_path_orca,sh_file=sh_file)
        t2 = time()
        with open("times.txt", "a") as file:
            file.write(f"{i} {t2 - t1}\n")
    except:
        print(f"directory {i} error")
        t2 = time()
        with open("times.txt", "a") as file:
            file.write(f"b_{i} {t2 - t1}\n")

#dir = "/home/rudenko/mopac_kie/KIE_MOPAC_PYHON/CH3NO2_H2O_near"
#calculate_dir(dir)
#orca_read_mult_from_mopac_input("opt.mop")
#orca_calculate_xyz("opted_h2o.xyz")

