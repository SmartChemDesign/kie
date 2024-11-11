#########################################################################
#                                                                       #
#     This module estimates kie from module: "new_MOPAC_calculator"     #
#     You should choose directory, then you will take kie.              #
#                                                                       #
#########################################################################
# based on ase 3.22.1
# necessary inputs
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import sys
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
import pandas as pd
import subprocess
import csv
import matplotlib.pyplot as plt

def createParser():
    """
    Function for creating parser for command line arguments

    Parameters
    ----------
    -p, --place : str
        Path to directory "calc"
        example: "mop/ch32chno2_ch3coo_new/calc"
    --plot : bool
        If True, then plot 2D graph of HOF
    --start_p : int
        From which point start drawing
        default: 0
    --ts : int
        Point number in transition state
    --isotope : str
        Isotope which is used for estimation
        default: "D"
    --end_p : int
        The point up to which the assessment will be made
        default: -1
    --step : int
        Step between points
        default: 0.25

    Returns
    -------
    "HOFS.pdf" : file with HOFS plot for TS determinating
    "HOFS.png" : file with HOFS plot
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--place', type=str, default="/")
    parser.add_argument( '--plot', type=bool, default=False)
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--ts', type=int)
    parser.add_argument('--isotope', type=str, default="D")
    parser.add_argument('--end_p', type=int, default=-1)
    parser.add_argument('--start_p', type=int, default=0)
    parser.add_argument('--step', type=int, default=0.25)
    return parser



if __name__ == '__main__':
    parser = createParser()
    namespace = parser.parse_args(sys.argv[1:])
    cur_dir = find_curent_directory()
    place = f"{cur_dir}{namespace.place}"
    if namespace.plot == True:
#--------------------------------Plot HOFS by point, for TS determination--------------------
        df = pd.read_csv(f"{place}/H_hof.csv", header=None)
        column_data = df.iloc[:, 0]
        column_list = column_data.to_list()[namespace.start_p:namespace.end_p]
        print(column_list)
        plt.scatter(range(len(column_list)),column_list)
        plt.title('HOFs')
        plt.ylabel("heat of formation, kcal")
        plt.xlabel("shift, A")
        plt.grid(which='major', linewidth = 1, alpha=0.2)
        plt.savefig(f"{place}/HOFS.pdf")
        plt.savefig(f"HOFS.pdf")
        
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        
        # Major ticks every 20, minor ticks every 5
        major_ticks_x = np.arange(0, len(column_list), 10)
        minor_ticks_x = np.arange(0, len(column_list), 2)
        major_ticks_y = np.arange(min(column_list) - 10, max(column_list) + 10, int((abs(max(column_list))+abs(min(column_list)))+1/10) )
        minor_ticks_y = np.arange(min(column_list) - 10,max(column_list) + 10, int((abs(max(column_list))+abs(min(column_list)))+1/50))
        
        
        ax.set_xticks(major_ticks_x)
        ax.set_xticks(minor_ticks_x, minor=True)
        ax.set_yticks(major_ticks_y)
        ax.set_yticks(minor_ticks_y, minor=True)
        
        # And a corresponding grid
        ax.grid(which='both')
        plt.scatter(range(len(column_list)), column_list)
        # Or if you want different settings for the grids:
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)
        plt.savefig(f"{place}/HOFS.pdf")
        plt.savefig(f"HOFS.pdf")
        plt.show()
        
        fig = plt.subplots()
        ax = plt.gca()
        plt.tight_layout(pad=3.5)
        y = column_list
        x = [namespace.step*i for i in range(1,1+len(y))]
        plt.scatter(x, y, linewidth=0)
        
        ax.set_title('HOFs')
        ax.set_ylabel("heat of formation, kcal")
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xlabel("shift, A")
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_major_formatter('{x:.1f}')

        plt.savefig(f"{place}HOFS.png", dpi=400)
        plt.savefig(f"HOFS.png", dpi=400)

        plt.close()
#---------------------------------------Plot Gibbs Energy-----------------------------------------------
        df_G = pd.read_csv(f"{place}/H_gibbs.csv", header=None)
        column_data_G = df_G.iloc[:, 0]
        column_list_G = column_data_G.to_list()[namespace.start_p:namespace.end_p]
        print(column_list_G)


        fig1 = plt.figure()
        ax1 = fig1.add_subplot(1, 1, 1)

        # Major ticks every 20, minor ticks every 5
        major_ticks_x_G = np.arange(0, len(column_list_G), 10)
        minor_ticks_x_G = np.arange(0, len(column_list_G), 2)
        major_ticks_y_G = np.arange(min(column_list_G) - 10, max(column_list_G) + 10,
                                  int((abs(max(column_list_G)) + abs(min(column_list_G))) + 1 / 10))
        minor_ticks_y_G = np.arange(min(column_list_G) - 10, max(column_list_G) + 10,
                                  int((abs(max(column_list_G)) + abs(min(column_list_G))) + 1 / 50))

        ax1.set_xticks(major_ticks_x_G)
        ax1.set_xticks(minor_ticks_x_G, minor=True)
        ax1.set_yticks(major_ticks_y_G)
        ax1.set_yticks(minor_ticks_y_G, minor=True)

        # And a corresponding grid
        ax1.grid(which='both')
        plt.scatter(range(len(column_list_G)), column_list_G)
        # Or if you want different settings for the grids:
        ax1.grid(which='minor', alpha=0.2)
        ax1.grid(which='major', alpha=0.5)
        ax1.set_title('HOFS_G')
        ax1.set_ylabel("heat of formation, kcal")
        ax1.xaxis.set_ticks_position('bottom')
        ax1.set_xlabel("shift, points")
        
        plt.savefig(f"{place}/HOFS_G.pdf")
        plt.savefig(f"HOFS_G.pdf")
        plt.show()



#-------------------Calculating KIE----------------------------------

    else:
        df_H = pd.read_csv(f"{place}/H.csv", header=None)
        df_H_column_data = df_H.iloc[:, 0]
        df_H_column_list = df_H_column_data.to_list()

        df_H_hof = pd.read_csv(f"{place}/H_hof.csv", header=None)
        df_H_hof_column_data = df_H_hof.iloc[:, 0]
        df_H_hof_column_list = df_H_hof_column_data.to_list()

        df_D = pd.read_csv(f"{place}/{namespace.isotope}.csv", header=None)
        df_D_column_data = df_D.iloc[:, 0]
        df_D_column_list = df_D_column_data.to_list()

        start_H = float(df_H_column_list[namespace.start_p])
        ts_H = float(df_H_column_list[namespace.ts])
        start_D = float(df_D_column_list[namespace.start_p])
        ts_D = float(df_D_column_list[namespace.ts])

        start_H_hof = float(df_H_hof_column_list[namespace.start_p])
        ts_H_hof = float(df_H_hof_column_list[namespace.ts])
        delta_H_hof = (ts_H_hof - start_H_hof) * 4.184

        print(ts_H)
        delta_H_kJ = (ts_H - start_H) * 4.184
        delta_D_kJ = (ts_D - start_D) * 4.184

        kie = np.exp(-(delta_H_kJ-delta_D_kJ)*1000/8.314/298)

        df_H_G = pd.read_csv(f"{place}/H_gibbs.csv", header=None)
        df_H_G_column_data = df_H_G.iloc[:, 0]
        df_H_G_column_list = df_H_G_column_data.to_list()

        df_D_G = pd.read_csv(f"{place}/{namespace.isotope}_gibbs.csv", header=None)
        df_D_G_column_data = df_D_G.iloc[:, 0]
        df_D_G_column_list = df_D_G_column_data.to_list()

        start_H_G = float(df_H_G_column_list[namespace.start_p])
        start_D_G = float(df_D_G_column_list[namespace.start_p])
        ts_H_G = float(df_H_G_column_list[namespace.ts])
        ts_D_G = float(df_D_G_column_list[namespace.ts])
        delta_H_kJ_G = (ts_H_G - start_H_G) / 1000
        delta_D_kJ_G = (ts_D_G - start_D_G) / 1000

        kie_G = np.exp(-(delta_H_kJ_G - delta_D_kJ_G) * 1000 / 8.314 / 298)
        print(f"kie_G = {kie_G}")

        print(f"kie = {kie}")
        with open (f"{place}/kie.txt", "w") as file:
            file.write(f"kie = {str(kie)}\n")
            file.write(f"kie_G = {str(kie_G)}\n")
            file.write(f"H_delta_zpe = {delta_H_kJ}\n")
            file.write(f"{namespace.isotope}_delta_zpe = {delta_D_kJ}\n")
            file.write(f"H_barier = {delta_H_hof}\n")

            file.write(f"H_barier_G = {delta_H_kJ_G}\n")
            file.write(f"{namespace.isotope}_barier_G = {delta_D_kJ_G}")






