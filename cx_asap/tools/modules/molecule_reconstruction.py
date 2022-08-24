#!/usr/bin/env python3

###################################################################################################
# ---------------------------------CX-ASAP: Molecule Reconstruction--------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config
import logging
import numpy as np
import pathlib
import os
from typing import Tuple

# ----------Class Definition----------#


class Molecule_Reconstruction:
    def __init__(
        self,
        reference_res: str,
        atoms: list,
        internal_structure: list,
        starting_atom: str,
        starting_coordinates: list,
        a_gradient: float,
        b_gradient: float,
        c_gradient: float,
        alpha_gradient: float,
        beta_gradient: float,
        gamma_gradient: float,
        max_position: int,
        min_position: int,
        step_size: int,
    ) -> None:

        """Initialises the class

        Sets up the yaml parameters input by the user

        Also defines the location for the system yaml file

        which stores a yaml of code-only parameters accessible throughout

        the software package

        Args:
            reference_res (str): full path to the reference .res file
            atoms (list): list of atoms present in the asymmetric unit
            internal_structure (list): list of bonds present (ie Cu1-O1)
            starting_atom (str): atom to start reconstruction with
            starting_coordinates (list): coordinates of starting atom (for special pos only)
            a_gradient (float): rate of change of a-axis
            b_gradient (float): rate of change of b-axis
            c_gradient (float): rate of change of c-axis
            alpha_gradient (float): rate of change of alpha angle
            beta_gradient (float): rate of change of beta angle
            gamma_gradient (float): rate of change of gamma angle
            max_position (int): maximum position to extrapolate out to
            min_position (int): minimum position to extrapoalte out to
            step_size (int): step size between extrapolated positions
        """

        # Setup yaml files and logger

        config = Config()

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

        self.path = reference_res
        self.atoms = atoms
        self.internal_structure = {}

        for item in internal_structure:
            self.internal_structure[item] = []

        self.starting_atom = starting_atom
        self.starting_coordinates = starting_coordinates
        self.a_gradient = float(a_gradient)
        self.b_gradient = float(b_gradient)
        self.c_gradient = float(c_gradient)
        self.alpha_gradient = float(alpha_gradient)
        self.beta_gradient = float(beta_gradient)
        self.gamma_gradient = float(gamma_gradient)

        try:
            self.positions = list(
                range(int(min_position), int(max_position), int(step_size))
            )
        except:
            print(
                "Invalid positions entered into conf.yaml. Please make sure the min_position, max_position and position_step size are all integers, and that the stepsize is smaller than the total range of values."
            )
            logging.critical(
                "Invalid positions entered into conf.yaml. Please make sure the min_position, max_position and position_step size are all integers, and that the stepsize is smaller than the total range of values."
            )
            exit()

        os.chdir(pathlib.Path(self.path).parent)

    def read_in_reference(self) -> None:

        """This function reads in the defined reference

        file and extracts the neutral unit cell and the

        neutral atomic coordinates

        """

        cell = []
        coordinates = {}
        with open(self.path, "rt") as f:
            lines = f.readlines()
        for i in lines:

            # Unit Cell extraction - saved to variable 'cell'

            if i.startswith("CELL"):
                for j, item in enumerate(i.split(" ")):
                    if j >= 2 and "\n" in item:
                        cell.append(float(item.strip("\n")))
                    elif j >= 2 and "\n" not in item and item != "":
                        cell.append(float(item))

            # Neutral coordinates of atoms extraction - saved to variable 'coordinates'

            for j in self.atoms:
                if i.startswith(j):
                    temp = i.split(" ")
                    temp2 = []
                    for item in temp:

                        # This extracts the values from the line that are just 'numbers' for easier identification of coordinates

                        try:
                            float(item)
                        except:
                            pass
                        else:
                            temp2.append(float(item))
                    coordinates[temp[0]] = temp2[1:4]

        self.neutral_cell = cell[1:]
        self.neutral_fractional_coordinates = coordinates

    def check(self, to_check: str, item: str, atom_list: str) -> list:

        """This function identifies atoms that are bonded to a

        specified atom and appends them to a list

        Args:
            to_check (str): an atom
            item (str): a bond
            atom_list (list): list of atoms bonded to the atom in 'to_check'

        Returns:
            atom_list (list): list of atoms bonded to the atom in 'to_check'
        """

        if to_check in item and to_check == item.split("-")[0]:
            atom_list.append(item.split("-")[1])
        elif to_check in item and to_check == item.split("-")[1]:
            atom_list.append(item.split("-")[0])

        return atom_list

    def find_construction_order(self) -> None:

        """This function determines the order that the atom

        will be reconstructed in the unit cell

        This is important, as each bond is represented as a vector,

        where the direction is important

        Each vector needs to be defined relative to the construction order,

        otherwise the atom will not be reconstructed correctly

        Ie, the direction of the vector will be incorrect,

        therefore the bonds will be in the incorrect position

        This is an artefact of the defined 'bonds' in the atom not

        necessarily being the same as the order of construction

        Ie the bond may be 'C1-C2', but the reconstruction will go from 'C2-C1'

        Rather than entering the bonds in a specific way,

        this function checks the construction order

        It will switch the order of vector subtraction to account for this

        """

        atoms_searched = [self.starting_atom]
        next_atoms = []
        new_next_atoms = []
        self.atomic_order = {}
        j = 0
        self.atomic_order[str(j)] = [self.starting_atom]
        j += 1

        # This loop iterates through all items in the list of bonds and finds any that contain the starting atom using the function 'check'

        # The 'atomic_order' variable is a dictionary that lists all of the atoms present, grouped together by when they are constructed

        # Ie, the first entry is the starting atom, the second entry is all atoms connected to the starting atom, etc

        for item in self.internal_structure:
            next_atoms = self.check(self.starting_atom, item, next_atoms)
            self.atomic_order[str(j)] = next_atoms
        for item in next_atoms:

            # This variable stores a list of all atoms searched - this will be used later to avoid duplicate atoms

            atoms_searched.append(item)

        # This loop functions essentially the same as the previous checking loop

        while len(atoms_searched) < len(self.neutral_fractional_coordinates):
            j += 1
            for i in next_atoms:
                for item in self.internal_structure:
                    new_next_atoms = self.check(i, item, new_next_atoms)

            # This makes sure that duplicate atoms are not present in the list (ie if one atom is bonded to two different ones)

            discrete_new_next_atoms = list(dict.fromkeys(new_next_atoms))

            # This removes any atoms that have already been searched

            # This ensures the correct 'direction' of molecule reconstruction

            # Ie, from the starting atom outwards, rather than back towards the starting atom

            # Lists are then set up for the next iteration of the loop

            for item in discrete_new_next_atoms:
                if item in atoms_searched:
                    discrete_new_next_atoms.remove(item)

            # Reapeating the above loop because if too many in a row it misses them?

            for item in discrete_new_next_atoms:
                if item in atoms_searched:
                    discrete_new_next_atoms.remove(item)

            for item in discrete_new_next_atoms:
                atoms_searched.append(item)
            next_atoms = discrete_new_next_atoms
            self.atomic_order[str(j)] = next_atoms
            new_next_atoms = []
            discrete_new_next_atoms = []

            print(self.atomic_order)

    def ordering(self, atom1, atom2) -> Tuple[str, str]:

        """This function iterates through the 'atomic_order'

        list and compares it to the atoms defined in the input

        It checks which atom comes first in the construction order,

        so that the coordinates can be subtracted in the correct order

        For vector subtraction, order is important!

        Args:
            atom1 (str): first atom for comparison
            atom2 (str): second atom for comparison

        Returns:
            first_atom (str): the first atom in the construction order
            second_atom (str): the second atom in the construction order

        """

        for i, j in enumerate(self.atomic_order):
            if atom1 in self.atomic_order[str(i)]:
                position1 = i
            elif atom2 in self.atomic_order[str(i)]:
                position2 = i

        if position1 < position2:
            first_atom = atom1
            second_atom = atom2
        else:
            first_atom = atom2
            second_atom = atom1

        return first_atom, second_atom

    def find_internal_vectors(self) -> None:

        """This function calculates an internal definition of the neutral atom using vectors"""

        neutral_coordinates = {}

        # The neutral coordinates need to be converted into real space first

        # This is done by multiplying the fractional coordinates by the cell parameters

        for item in self.neutral_fractional_coordinates:
            neutral_coordinates[item] = np.array(
                [
                    [
                        self.neutral_fractional_coordinates[item][0]
                        * self.neutral_cell[0],
                        self.neutral_fractional_coordinates[item][1]
                        * self.neutral_cell[1],
                        self.neutral_fractional_coordinates[item][2]
                        * self.neutral_cell[2],
                    ]
                ]
            )

        # Vector subtraction to define the bonds based on the construction order

        for item in self.internal_structure:
            atoms = item.split("-")
            A1, A2 = self.ordering(atoms[0], atoms[1])
            atom_1 = neutral_coordinates[A1]
            atom_2 = neutral_coordinates[A2]
            self.internal_structure[item] = np.subtract(atom_2, atom_1)

    def calculate_cell_parameter(
        self, gradient: float, position: int, neutral: float
    ) -> float:

        """This function will calculate a new cell parameter

        based on a gradient and an arbitrary 'position'

        Args:
            gradient (float): the rate of change of the input parameter
            position (int): the position number
            neutral (float): the neutral value of the input parameter

        Returns:
            parameter (float): the new value taking position and rate of change into account
        """

        deformation = gradient * position
        parameter = ((deformation * neutral) / 100) + neutral

        return abs(parameter)

    def calculate_new_cell(self, position: int) -> None:

        """This function will define the new cell,

        calling the calculate_cell_parameter function for each variable

        Args:
            position (int): The position number
        """

        #

        new_cell_a = self.calculate_cell_parameter(
            self.a_gradient, position, self.neutral_cell[0]
        )
        new_cell_b = self.calculate_cell_parameter(
            self.b_gradient, position, self.neutral_cell[1]
        )
        new_cell_c = self.calculate_cell_parameter(
            self.c_gradient, position, self.neutral_cell[2]
        )
        new_cell_alpha = self.calculate_cell_parameter(
            self.alpha_gradient, position, self.neutral_cell[3]
        )
        new_cell_beta = self.calculate_cell_parameter(
            self.beta_gradient, position, self.neutral_cell[4]
        )
        new_cell_gamma = self.calculate_cell_parameter(
            self.gamma_gradient, position, self.neutral_cell[5]
        )

        self.new_cell = [
            new_cell_a,
            new_cell_b,
            new_cell_c,
            new_cell_alpha,
            new_cell_beta,
            new_cell_gamma,
        ]

    def calculate_fractional_coordinates(self) -> None:

        """This function will reconstruct the molecule in the specified new unit cell"""

        new_molecule = {}
        self.new_fractional_coordinates = {}
        self.new_fractional_coordinates[self.starting_atom] = self.starting_coordinates

        # First, the internal vectors are converted into fractional coordinates in the new cell

        for item in self.internal_structure:
            new_molecule[item] = []
            new_molecule[item].append(
                self.internal_structure[item][0][0] / self.new_cell[0]
            )
            new_molecule[item].append(
                self.internal_structure[item][0][1] / self.new_cell[1]
            )
            new_molecule[item].append(
                self.internal_structure[item][0][2] / self.new_cell[2]
            )

        # This series of loops reconstructs the molecule based on the order previously determined using vector addition

        for i in self.atomic_order:
            for j in self.atomic_order[i]:
                for item in self.internal_structure:
                    if j in item.split("-"):
                        try:
                            self.atomic_order[str(int(i) + 1)]
                        except KeyError:
                            pass
                        else:
                            for k in self.atomic_order[str(int(i) + 1)]:
                                if k in item.split("-"):
                                    self.new_fractional_coordinates[k] = np.add(
                                        self.new_fractional_coordinates[j],
                                        new_molecule[item],
                                    )

    def write_res(self, item: int) -> None:

        """This function writes a .res file as an output

        based on the reference - only the cell and coordinates are changed

        This should not be used to refine structures!

        This is purely to provide files for hypothetical calculations

        Args:
            item (int): The position number
        """

        with open(self.path, "rt") as f:
            lines = f.readlines()

        try:
            os.mkdir("Extrapolated_Files")
        except:
            pass

        os.chdir("Extrapolated_Files")

        with open("Position_" + str(item) + ".res", "w") as f:
            for i in lines:
                split = i.split(" ")
                if i.startswith("CELL"):

                    if split[1] != "":

                        new_line = (
                            split[0]
                            + " "
                            + split[1]
                            + " "
                            + str(round(self.new_cell[0], 4))
                            + " "
                            + str(round(self.new_cell[1], 4))
                            + " "
                            + str(round(self.new_cell[2], 4))
                            + " "
                            + str(round(self.new_cell[3], 4))
                            + " "
                            + str(round(self.new_cell[4], 4))
                            + " "
                            + str(round(self.new_cell[5], 4))
                            + " \n"
                        )

                    else:
                        new_line = (
                            split[0]
                            + " "
                            + split[2]
                            + " "
                            + str(round(self.new_cell[0], 4))
                            + " "
                            + str(round(self.new_cell[1], 4))
                            + " "
                            + str(round(self.new_cell[2], 4))
                            + " "
                            + str(round(self.new_cell[3], 4))
                            + " "
                            + str(round(self.new_cell[4], 4))
                            + " "
                            + str(round(self.new_cell[5], 4))
                            + " \n"
                        )

                    f.write(new_line)

                elif split[0] in self.new_fractional_coordinates:
                    j = 0
                    for index, char in enumerate(split):
                        if char != "" and j in range(2, 5):
                            if j == 2:
                                if split[0] == self.starting_atom:
                                    split[index] = str(
                                        self.new_fractional_coordinates[split[0]][0]
                                    )
                                else:
                                    split[index] = str(
                                        round(
                                            self.new_fractional_coordinates[split[0]][
                                                0
                                            ],
                                            7,
                                        )
                                    )
                            elif j == 3:
                                if split[0] == self.starting_atom:
                                    split[index] = str(
                                        self.new_fractional_coordinates[split[0]][1]
                                    )
                                else:
                                    split[index] = str(
                                        round(
                                            self.new_fractional_coordinates[split[0]][
                                                1
                                            ],
                                            7,
                                        )
                                    )
                            else:
                                if split[0] == self.starting_atom:
                                    split[index] = str(
                                        self.new_fractional_coordinates[split[0]][2]
                                    )
                                else:
                                    split[index] = split[index] = str(
                                        round(
                                            self.new_fractional_coordinates[split[0]][
                                                2
                                            ],
                                            7,
                                        )
                                    )
                            j += 1
                        elif char != "" and j != 2:
                            j += 1
                    new_line = " ".join(split)
                    f.write(new_line)
                else:
                    f.write(i)

        os.chdir(pathlib.Path(self.path).parent)
