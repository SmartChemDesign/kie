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

def createParser():
    """
    Function for creating parser for command line arguments

    Parameters
    ----------
    -x, --xyz : str
        Path to xyz file with optimized structure
    -m, --moving_atom_h : int
        Number of moving H atom
    -d, --dest_atom : int
        Number of atom, to which H is moving
    -e, --eps : float
        Dielectric constant
    -s, --step : float
        Step for H moving atom
    --mult : str
        Multiplicity of molecule
    --charge : int
        Charge of molecule
    --isotope : str
        Isotope which is used for estimation
    --other_isotopes : str
        List of other isotopes atoms
    -n, --h_neighbour : int
        Number of atom, from which H is moving
    -p, --place : str
        Path to directory, where molecule is located
    --freez : str
        Atoms, which should be freezed
    --freez_all_axis : bool
        If True, then all atoms will be freezed
    """
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
    parser.add_argument('-p', '--place', type=str, default="proposed_mopac/")  # H moving from
    parser.add_argument('--freez', nargs='+', help='<Required> Set flag', required=False, default="")
    parser.add_argument('--freez_all_axis', type=bool, default=False)
    parser.add_argument('-t', '--threads', type=bool, default=100)
    return parser


if __name__ == '__main__':
    # Parsing command line arguments
    parser = createParser()
    namespace = parser.parse_args(sys.argv[1:])
    cur_dir = find_curent_directory()
    freez = str_array_to_int_array(namespace.freez)
    place = f"{cur_dir}{namespace.place}"
    
    # Define isotopes atoms
    if namespace.other_isotopes != None:
        other_isotopes = str_array_to_int_array(namespace.other_isotopes)
        if namespace.moving_atom_h not in other_isotopes:
            other_isotopes.append(namespace.moving_atom_h)
        isotope_positions = other_isotopes
    else:
        isotope_positions = [namespace.moving_atom_h]
    print(f"isotope_positions {isotope_positions}")
    #print(f"other_isotopes {other_isotopes}")

    # Preopt geom
    subs = namespace.place.split('/')[-2]
    not_opted_geom = place[:-5] + "/" + "start" + "/" + "not_opted.xyz"
    opted_h2o = place[:-5] + "/" + "start" + "/" + "opted_h2o.xyz"
    opt_geom(filename_inp=not_opted_geom, filename_xyz = opted_h2o, trj=False, xyz_trjry="error.txt",
             charge = namespace.charge, mult=namespace.mult, eps=namespace.eps,
             h_neighbour=namespace.h_neighbour, h_num=namespace.moving_atom_h,
             dest_num=namespace.dest_atom, freez=freez, freez_all_axis=namespace.freez_all_axis,
             threads=namespace.threads)
    print("preopt_done")
    # Copy opted geom to /calc
    shutil.copyfile(f"{opted_h2o}", f"{place}/opted_h2o.xyz")


    # Read optimized structure 
    # Rotate molecule and define H atom direction
    mol_xyz_1 = read(f"{place}/{namespace.xyz}")
    mol_xyz = rotate_xyz_to_x(mol_xyz_1,namespace.moving_atom_h, namespace.dest_atom)
    distance = mol_xyz.get_distance(namespace.moving_atom_h, namespace.dest_atom, vector=False)
    distance_vector = np.array(mol_xyz.get_distance(namespace.moving_atom_h, namespace.dest_atom, vector=True))
    step_vector = distance_vector / (distance/namespace.step)
    
    # Create files
    first_geom = place + "/" + namespace.xyz
    mol_xyz.write(f"{first_geom[:-4]}.xyz")
    opted_geom = place + "/" + "opted_in_circle.xyz"
    thermo_name = place + "/" + "thermo_in_circle.xyz"
    xyz_trjry = place + "/" + "trj.xyz"
    create_files(place)
    
    # Calculate time
    t1 = time()
    distance_n = mol_xyz.get_distance(namespace.moving_atom_h, namespace.dest_atom, vector=False)
    print(distance_n)
    n = 0
    
    
    # Start optimizing circle
    
    while distance > 1:
        n += 1 
        print(n, " dist", distance)
        # Optimize geometry
        opt_geom(filename_inp=first_geom, filename_xyz=opted_geom, xyz_trjry=xyz_trjry,
                 charge = namespace.charge, mult=namespace.mult, eps=namespace.eps,
                 h_neighbour=namespace.h_neighbour, h_num=namespace.moving_atom_h,
                 dest_num=namespace.dest_atom, freez=freez, freez_all_axis=namespace.freez_all_axis,
                 threads=namespace.threads)
        # Find thermo for optimized geometry
        find_first_thermo(opted_geom = opted_geom, thermo_name = thermo_name,
                          charge = namespace.charge, mult=namespace.mult,
                          eps=namespace.eps, place=place, threads=namespace.threads)
        # Find thermo for isotope
        find_isotopes_thermo(opted_geom=opted_geom, thermo_name=thermo_name,
                             positions= isotope_positions,
                             isotope=namespace.isotope, charge = namespace.charge,
                             mult=namespace.mult, eps=namespace.eps, place=place,
                             threads=namespace.threads)
        #find_isotopes_thermo(opted_geom=opted_geom, thermo_name=thermo_name, positions= [namespace.moving_atom_h],
        # isotope="T", charge = namespace.charge, mult=namespace.mult, eps=namespace.eps)

        # Rotate geometry
        mol_xyz, step_vector, distance, distance_n = find_move_vectore(first_geom=first_geom, moving_atom_h=namespace.moving_atom_h,
                                                           dest_atom=namespace.dest_atom, step=namespace.step, h_neighbour=namespace.h_neighbour)
        mol_xyz.write(f"{first_geom[:-4]}.xyz")
        #mol_xyz.write(f"{opted_geom[:-4]}.xyz")
        # Move atom H 
        move_atom(opted_geom=opted_geom, atom =namespace.moving_atom_h,
                  vector=step_vector, first_geom=first_geom)
        
        # Next step

    # Calculate time
    t2 = time()
    with open(f"{place}/time.txt", "w") as file:
        file.write(f"{place} {t2 - t1}\n")


