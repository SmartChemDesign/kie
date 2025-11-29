#!/usr/bin/env python3
#########################################################################
#                                                                       #
#  This script reoptimizes TS from existing NEB calculations           #
#  and recalculates frequencies for KIE estimation                     #
#                                                                       #
#########################################################################

import sys
import argparse
import numpy as np
import os
import shutil
import subprocess
import datetime
from time import time
from ase.io import read, write
from ase.atoms import Atoms
from ase.units import Hartree, Bohr
from ase.io.orca import write_orca
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError
from ase.calculators.orca import ORCA


def create_folder_my(folder):
    """Create folder, remove if exists"""
    try:
        os.mkdir(folder)
    except FileExistsError:
        shutil.rmtree(folder)
        os.mkdir(folder)


class My_ORCA(ORCA):
    """Custom ORCA calculator class"""
    def __init__(self, restart=None,
                 ignore_bad_restart_file=FileIOCalculator._deprecated,
                 label='ab', atoms=None, **kwargs):
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)
        self.pcpot = None

    def calculate_not_create_input(self, path_to_orca):
        """Run ORCA calculation without creating input"""
        command = f"{path_to_orca} {self.label}.inp > {self.label}.out"
        print(f"Running: {command}")
        try:
            proc = subprocess.Popen(command, shell=True, cwd=self.directory)
        except OSError as err:
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err
        errorcode = proc.wait()

    def write_input(self, atoms, properties=None, system_changes=None):
        """Write ORCA input file"""
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters
        p.write(self.label + '.ase')
        p['label'] = self.label
        if self.pcpot:
            p['pcpot'] = self.pcpot

        write_orca(atoms, **p)

        with open(self.label + '.inp', "r") as bad_input:
            inp = bad_input.readlines()
            inp[0] = inp[0].replace("engrad", "")

        with open(self.label + '.inp', "w") as good_input:
            for line in inp:
                good_input.write(line)


def check_norm_end(file):
    """Check if ORCA calculation finished normally"""
    try:
        with open(file, "r", encoding="UTF-8") as out_res:
            res = out_res.readlines()
        for i in res:
            if "NORMALLY" in i:
                return True
    except FileNotFoundError:
        print(f"File {file} not found")
    return False


def find_isotope_number_from_sh(sh_file):
    """Find isotope atom number from .sh file"""
    try:
        with open(sh_file, "r", encoding="UTF-8") as bad_input:
            inp = bad_input.readlines()

        for i, line in enumerate(inp):
            if "main.py" in line:
                st = i
        data_inp = inp[st].strip().split(" ")

        for k, j in enumerate(data_inp):
            if "--moving_atom_h" in j.lower():
                d = data_inp[k+1]
        print(f"Found isotope atom number: {int(d)}")
        return int(d)
    except:
        print(f"Warning: Could not find isotope number from {sh_file}")
        return None


def orca_read_mult_from_sh(sh_file):
    """Read charge, multiplicity and epsilon from .sh file"""
    with open(sh_file, "r", encoding="UTF-8") as bad_input:
        inp = bad_input.readlines()

    eps = 78.4
    charge = 0
    mult = "SINGLET"
    
    for i, line in enumerate(inp):
        if "main.py" in line:
            st = i
    data_inp = inp[st].strip().split(" ")
    
    for k, j in enumerate(data_inp):
        if j.lower() == "--eps":
            eps = data_inp[k+1]
        if j.lower() == "--charge":
            charge = data_inp[k+1]
        if j.lower() == "--mult":
            mult = data_inp[k+1]
    
    mult_dict = {"SINGLET": 1, "DOUBLET": 2}
    
    return int(charge), mult_dict[str(mult)], float(eps)


def orca_change_isotope(file, position, isotope):
    """Change isotope in ORCA input file"""
    with open(file, "r", encoding="UTF-8") as bad_input:
        start = 0
        inp = bad_input.readlines()
    
    for i, j in enumerate(inp):
        if "*xyz" in j:
            start = i+1

    print(f"Changing isotope at position {int(start+position)}")
    inp[int(start + position)] = inp[int(start + position)].strip() + f" M = {isotope}" + "\n"
    
    with open(file, "w") as good_input:
        for line in inp:
            good_input.write(line)
    
    return 0


def orca_calculate_xyz(file, label, lom_path_orca,
                       change_isotope=False,
                       sh_file="none.sh",
                       isotope="2",
                       charge=-1,
                       mult=1,
                       eps=78.4,
                       maxiter=2000,
                       orcasimpleinput='opt B3LYP def2-TZVP CPCM D3',
                       orcablocks='%pal nprocs 102 end'):
    """Calculate using ORCA"""
    
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
        if isotope_num is not None:
            orca_change_isotope(label+".inp", isotope_num, isotope)

    print("Starting calculation...")
    calc_orca.calculate_not_create_input(lom_path_orca)


def take_xyz_from_dir(dir, flag):
    """Find xyz file in directory"""
    files = os.listdir(dir)
    for i in files:
        if flag in i:
            return f"{dir}/{i}"
    return None


def reoptimize_and_recalc_freqs(dir, lom_path_orca, sh_file, nprocs=102):
    """
    Reoptimize TS from existing NEB and recalculate frequencies
    
    Args:
        dir: directory with calculation results
        lom_path_orca: path to ORCA executable
        sh_file: corresponding .sh file
        nprocs: number of processors
    """
    
    print(f"\n{'='*70}")
    print(f"Processing directory: {dir}")
    print(f"{'='*70}\n")
    
    # Check if NEB results exist
    neb_xyz = f"{dir}/orca/neb/neb.xyz"
    if not os.path.exists(neb_xyz):
        print(f"ERROR: NEB results not found at {neb_xyz}")
        return False
    
    # Read parameters from .sh file
    try:
        start_charge, start_mult, start_eps = orca_read_mult_from_sh(sh_file)
        print(f"Parameters: charge={start_charge}, mult={start_mult}, eps={start_eps}")
    except:
        print(f"ERROR: Could not read parameters from {sh_file}")
        return False
    
    # Prepare orcablocks with specified number of processors
    orcablocks = f'%pal nprocs {nprocs} end'
    
    # ============= TS Optimization =============
    print("\n--- Step 1: TS Optimization ---")
    create_folder_my(f"{dir}/orca/ts_opt")
    shutil.copy(neb_xyz, f"{dir}/orca/ts_opt/ts_init.xyz")
    
    ts_init_xyz = take_xyz_from_dir(f"{dir}/orca/ts_opt", ".xyz")
    orca_calculate_xyz(ts_init_xyz, 
                       lom_path_orca=lom_path_orca, 
                       label=f"{dir}/orca/ts_opt/ts_opt",
                       charge=start_charge,
                       mult=start_mult,
                       eps=start_eps,
                       orcasimpleinput='OptTS B3LYP def2-TZVP CPCM D3',
                       orcablocks=orcablocks)
    
    # Check if TS optimization finished normally
    if not check_norm_end(f"{dir}/orca/ts_opt/ts_opt.out"):
        print("WARNING: TS optimization did not finish normally!")
        return False
    
    print("TS optimization completed successfully!")
    
    # ============= Frequency calculations for TS =============
    print("\n--- Step 2: Frequency calculations for TS ---")
    
    # Create freq_folders
    create_folder_my(f"{dir}/orca/freq_H")
    create_folder_my(f"{dir}/orca/freq_D")
    
    # Copy optimized TS structure
    shutil.copy(f"{dir}/orca/ts_opt/ts_opt.xyz", f"{dir}/orca/freq_H/ts_opt.xyz")
    shutil.copy(f"{dir}/orca/ts_opt/ts_opt.xyz", f"{dir}/orca/freq_D/ts_opt.xyz")
    
    # Run freq_H
    print("\nCalculating frequencies for H isotope...")
    freq_H_xyz = take_xyz_from_dir(f"{dir}/orca/freq_H", ".xyz")
    orca_calculate_xyz(freq_H_xyz, 
                       lom_path_orca=lom_path_orca, 
                       label=f"{dir}/orca/freq_H/f_H",
                       charge=start_charge,
                       mult=start_mult, 
                       eps=start_eps,
                       orcasimpleinput='NUMFREQ B3LYP def2-TZVP CPCM D3',
                       orcablocks=orcablocks)
    
    # Run freq_D
    print("\nCalculating frequencies for D isotope...")
    freq_D_xyz = take_xyz_from_dir(f"{dir}/orca/freq_D", ".xyz")
    orca_calculate_xyz(freq_D_xyz, 
                       lom_path_orca=lom_path_orca,
                       label=f"{dir}/orca/freq_D/f_D", 
                       charge=start_charge,
                       mult=start_mult,
                       eps=start_eps,
                       orcasimpleinput='NUMFREQ B3LYP def2-TZVP CPCM D3', 
                       isotope="2",
                       change_isotope=True, 
                       sh_file=sh_file,
                       orcablocks=orcablocks)
    
    # ============= Frequency calculations for Start =============
    print("\n--- Step 3: Frequency calculations for Start ---")
    
    # Check if start structure exists
    start_xyz = f"{dir}/orca/orca_start/start.xyz"
    if not os.path.exists(start_xyz):
        print(f"ERROR: Start structure not found at {start_xyz}")
        return False
    
    # Create freq_folders
    create_folder_my(f"{dir}/orca/start_freq_H")
    create_folder_my(f"{dir}/orca/start_freq_D")
    shutil.copy(start_xyz, f"{dir}/orca/start_freq_H/neb.xyz")
    shutil.copy(start_xyz, f"{dir}/orca/start_freq_D/neb.xyz")
    
    # Run start_freq_H
    print("\nCalculating start frequencies for H isotope...")
    start_freq_H_xyz = take_xyz_from_dir(f"{dir}/orca/start_freq_H", ".xyz")
    orca_calculate_xyz(start_freq_H_xyz, 
                       lom_path_orca=lom_path_orca, 
                       label=f"{dir}/orca/start_freq_H/start_f_H",
                       charge=start_charge,
                       mult=start_mult, 
                       eps=start_eps,
                       orcasimpleinput='NUMFREQ B3LYP def2-TZVP CPCM D3',
                       orcablocks=orcablocks)
    
    # Run start_freq_D
    print("\nCalculating start frequencies for D isotope...")
    start_freq_D_xyz = take_xyz_from_dir(f"{dir}/orca/start_freq_D", ".xyz")
    orca_calculate_xyz(start_freq_D_xyz, 
                       lom_path_orca=lom_path_orca,
                       label=f"{dir}/orca/start_freq_D/start_f_D", 
                       charge=start_charge,
                       mult=start_mult,
                       eps=start_eps,
                       orcasimpleinput='NUMFREQ B3LYP def2-TZVP CPCM D3', 
                       isotope="2",
                       change_isotope=True, 
                       sh_file=sh_file,
                       orcablocks=orcablocks)
    
    print(f"\n{'='*70}")
    print(f"Successfully completed reoptimization for: {dir}")
    print(f"{'='*70}\n")
    
    return True


def createParser():
    """Create argument parser"""
    parser = argparse.ArgumentParser(
        description='Reoptimize TS from existing NEB calculations and recalculate frequencies')
    parser.add_argument('-p', '--path', type=str, required=True,
                       help='Path to directory with for_orca folder')
    parser.add_argument('-o', '--orca', type=str, 
                       default='/home/rudenko/orca/orca',
                       help='Path to ORCA executable (default: /home/rudenko/orca/orca)')
    parser.add_argument('-n', '--nprocs', type=int, default=102,
                       help='Number of processors for ORCA (default: 102)')
    parser.add_argument('-d', '--dir', type=str, default=None,
                       help='Specific directory to process (if not specified, process all)')
    return parser


if __name__ == '__main__':
    parser = createParser()
    args = parser.parse_args()
    
    # Setup paths
    home_path = args.path
    lom_path_orca = args.orca
    nprocs = args.nprocs
    
    path_to_orca_dirs = os.path.join(home_path, "for_orca")
    
    if not os.path.exists(path_to_orca_dirs):
        print(f"ERROR: Directory {path_to_orca_dirs} does not exist!")
        sys.exit(1)
    
    # Get list of directories to process
    if args.dir:
        orca_dirs = [args.dir]
        print(f"Processing single directory: {args.dir}")
    else:
        orca_dirs = os.listdir(path_to_orca_dirs)
        print(f"Found {len(orca_dirs)} directories to process")
    
    # Start processing
    log_file = os.path.join(home_path, "reoptimization_times.txt")
    with open(log_file, "w") as file:
        file.write(f"Reoptimization started at: {datetime.datetime.now()}\n")
        file.write(f"ORCA path: {lom_path_orca}\n")
        file.write(f"Number of processors: {nprocs}\n")
        file.write(f"{'='*70}\n")
    
    success_count = 0
    failed_dirs = []
    
    for i in orca_dirs:
        sh_file = os.path.join(home_path, f"{i}.sh")
        dir_path = os.path.join(path_to_orca_dirs, i)
        
        # Check if .sh file exists
        if not os.path.exists(sh_file):
            print(f"WARNING: .sh file not found for {i}, skipping...")
            continue
        
        # Check if directory exists
        if not os.path.exists(dir_path):
            print(f"WARNING: Directory {dir_path} not found, skipping...")
            continue
        
        t1 = time()
        try:
            success = reoptimize_and_recalc_freqs(dir_path, lom_path_orca, sh_file, nprocs)
            t2 = time()
            
            if success:
                success_count += 1
                with open(log_file, "a") as file:
                    file.write(f"SUCCESS: {i} - {t2 - t1:.2f} seconds\n")
            else:
                failed_dirs.append(i)
                with open(log_file, "a") as file:
                    file.write(f"FAILED: {i} - {t2 - t1:.2f} seconds\n")
        except Exception as e:
            print(f"ERROR processing directory {i}: {str(e)}")
            t2 = time()
            failed_dirs.append(i)
            with open(log_file, "a") as file:
                file.write(f"ERROR: {i} - {str(e)} - {t2 - t1:.2f} seconds\n")
    
    # Write summary
    with open(log_file, "a") as file:
        file.write(f"\n{'='*70}\n")
        file.write(f"Reoptimization completed at: {datetime.datetime.now()}\n")
        file.write(f"Total directories processed: {len(orca_dirs)}\n")
        file.write(f"Successful: {success_count}\n")
        file.write(f"Failed: {len(failed_dirs)}\n")
        if failed_dirs:
            file.write(f"\nFailed directories:\n")
            for d in failed_dirs:
                file.write(f"  - {d}\n")
    
    print(f"\n{'='*70}")
    print(f"SUMMARY:")
    print(f"Total directories: {len(orca_dirs)}")
    print(f"Successful: {success_count}")
    print(f"Failed: {len(failed_dirs)}")
    print(f"Log file: {log_file}")
    print(f"{'='*70}\n")
