from ase.calculators.mopac import MOPAC
from ase.calculators.calculator import FileIOCalculator, ReadError,Calculator
import subprocess
import os
import re
from typing import Sequence
from warnings import warn

import numpy as np
from packaging import version

from ase import Atoms
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError
from ase.units import Debye, kcal, mol

def find_curent_directory():
    if "\\" in os.getcwd():
        return os.getcwd()
    else:
        return os.getcwd()

class MYMOPAC(MOPAC):
    def get_zero_point_energy(self):
        return self.zero_point_energy
    def get_my_hof(self):
        return self.my_hof

    def change_isotopes_input(self, positions, isotope):  # с нуля
        with open(self.label + '.mop', "r") as fd:
            lines = fd.readlines()
            for atom in positions:
                if "H" in lines[3 + atom]:
                    # print(1)
                    lines[3 + atom] = lines[4 + atom].replace("H", isotope)
                elif "D" in lines[4 + atom]:
                    lines[3 + atom] = lines[4 + atom].replace("D", isotope)
                elif "T" in lines[4 + atom]:
                    lines[3 + atom] = lines[4 + atom].replace("T", isotope)
            # print(lines)
        os.remove(f"{self.label}.mop")
        with open(self.label + '.mop', "w") as fd1:
            fd1.writelines(lines)

    def calculate_not_create_input(self, atoms=None, properties=['energy']):
        Calculator.calculate(self, atoms, properties)
        if self.command is None:
            raise CalculatorSetupError(
                'Please set ${} environment variable '
                .format('ASE_' + self.name.upper() + '_COMMAND') +
                'or supply the command keyword')
        command = self.command
        if 'PREFIX' in command:
            command = command.replace('PREFIX', self.prefix)
        print(f"---{command}---")
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

        if errorcode:
            path = os.path.abspath(self.directory)
            msg = ('Calculator "{}" failed with command "{}" failed in '
                   '{} with error code {}'.format(self.name, command,
                                                  path, errorcode))
            raise CalculationFailed(msg)

        self.read_results()



    def change_isotopes_input(self, positions, isotope):  # с нуля
        with open(self.label + '.mop', 'r') as fd:
            lines = fd.readlines()
            for atom in positions:
                if "H" in lines[3 + atom]:
                    # print(1)
                    lines[3 + atom] = lines[3 + atom].replace("H", isotope)
                elif "D" in lines[3 + atom]:
                    lines[3 + atom] = lines[3 + atom].replace("D", isotope)
                elif "T" in lines[3 + atom]:
                    lines[3 + atom] = lines[3 + atom].replace("T", isotope)
            # print(lines)
        os.remove(f"{self.label}.mop")
        with open(self.label + '.mop', "w") as fd1:
            fd1.writelines(lines)

        return 0

    default_parameters = dict(
        method='PM7',
        task='1SCF GRADIENTS',
        charge=0,
        relscf=1)
    def read_results(self):
        """Read the results, such as energy, forces, eigenvalues, etc.
        """
        FileIOCalculator.read(self, self.label)
        if not os.path.isfile(self.label + '.out'):
            raise ReadError

        with open(self.label + '.out') as fd:
            lines = fd.readlines()

        for i, line in enumerate(lines):
            if line.find('TOTAL ENERGY') != -1:
                self.results['energy'] = float(line.split()[3])
            if line.find('ZERO POINT ENERGY    ') != -1:
                self.zero_point_energy = float(line.split()[3])
            if line.find('FINAL HEAT OF FORMATION') != -1:
                self.final_hof = float(line.split()[5]) * kcal / mol
            elif line.find('HEAT OF FORMATION') != -1:
               	print(line.split())
                try:
                    self.my_hof = float(line.strip().split()[4])
                except:
                    self.my_hof = 0.000
            elif line.find('NO. OF FILLED LEVELS') != -1:
                self.nspins = 1
                self.no_occ_levels = int(line.split()[-1])
            elif line.find('NO. OF ALPHA ELECTRON') != -1:
                self.nspins = 2
                self.no_alpha_electrons = int(line.split()[-1])
                self.no_beta_electrons = int(lines[i+1].split()[-1])
                self.results['magmom'] = abs(self.no_alpha_electrons -
                                             self.no_beta_electrons)
            elif line.find('EIGENVALUES') != -1:
                if line.find('ALPHA') != -1:
                    j = i + 1
                    eigs_alpha = []
                    while not lines[j].isspace():
                        eigs_alpha += [float(eps) for eps in lines[j].split()]
                        j += 1
                elif line.find('BETA') != -1:
                    j = i + 1
                    eigs_beta = []
                    while not lines[j].isspace():
                        eigs_beta += [float(eps) for eps in lines[j].split()]
                        j += 1
                    eigs = np.array([eigs_alpha, eigs_beta]).reshape(2, 1, -1)
                    self.eigenvalues = eigs
                else:
                    eigs = []
                    j = i + 1
                    while not lines[j].isspace():
                        eigs += [float(e) for e in lines[j].split()]
                        j += 1
                    self.eigenvalues = np.array(eigs).reshape(1, 1, -1)
            #elif line.find('DIPOLE   ') != -1:
                #self.results['dipole'] = np.array(
                    #lines[i + 3].split()[1:1 + 3], float) * Debye

    def read_atoms_from_out(self, lines):
        """Read the Atoms from the output file stored as list of str in lines.
        Parameters:

            lines: list of str
        """
        # first try to read from final point (last image)
        i = self.get_index(lines, '   NUMBER   SYMBOL  FORCE CONSTANT  FORCE CONSTANT  FORCE CONSTANT')
        if i is None:  # XXX should we read it from the input file?
            assert 0, 'Not implemented'

        lines1 = lines[i:]
        j = 2
        symbols = []
        positions = []
        while not lines1[j].isspace():  # continue until we hit a blank line
            l = lines1[j].split()
            symbols.append(l[1])
            positions.append([float(c) for c in l[2: 2 + 3]])
            j += 1

        return Atoms(symbols=symbols, positions=positions)

    def read_atoms_from_out_opt(self, lines):
        """Read the Atoms from the output file stored as list of str in lines.
        Parameters:

            lines: list of str
        """
        # first try to read from final point (last image)
        i = self.get_index(lines, '                             CARTESIAN COORDINATES')
        if i is None:  # XXX should we read it from the input file?
            assert 0, 'Not implemented'

        lines1 = lines[i:]
        j = 2
        symbols = []
        positions = []
        while not lines1[j].isspace():  # continue until we hit a blank line
            l = lines1[j].split()
            symbols.append(l[1])
            print(l)
            positions.append([float(c) for c in l[2: 2 + 3]])
            j += 1


        return Atoms(symbols=symbols, positions=positions)

