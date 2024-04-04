#!/usr/bin/env python3

###################################################################################################
# -------------------------------------CX-ASAP: rotation_planes------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config
import pathlib
import os
import math
import pandas as pd
import logging
import numpy as np

# ----------Class Definition----------#


class Rotation:
    def __init__(self, test_mode: bool = False) -> None:
        """Initialises the class

        Sets up the yaml parameters input by the user

        Also defines the location for the system yaml file

        which stores a yaml of code-only parameters accessible throughout

        the software package

        Args:
            test_mode (bool): Automatically false, if true it will

                            make the functions compatible with the testing script
        """

        # Setup yaml files and logger

        self.test_mode = test_mode

        config = Config(self.test_mode)

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    def configure(self, ref_plane: list) -> None:
        """Checks that the user has input a valid reference plane

        Args:
            ref_plane (list): plane to compare MPLA against as a list
                            ie [1,0,0] corresponds to the (100) plane
        """

        self.ref_plane = []

        for item in ref_plane:
            if str(item).isdigit():
                self.ref_plane.append(int(item))

        if len(self.ref_plane) != 3:
            logging.critical(
                __name__ + " : Reference plane in config file not understood by system"
            )
            print("Error - check logs!")
            exit()

    def grab_cell(self, file_name: str) -> None:
        """Collects the unit cell for later calculations

        Args:
            file_name (str): full path to the file containing the unit cell info
        """

        # collects the cell from that run to calculate the rotation against

        with open(file_name, "rt") as ins_file:
            split_line = []
            for line in ins_file:
                if "CELL" in line:
                    split_line = line.split()

        if len(split_line) != 0:
            ref_parameters = [
                "ref_INS_a",
                "ref_INS_b",
                "ref_INS_c",
                "ref_INS_alpha",
                "ref_INS_beta",
                "ref_INS_gamma",
            ]

            self.ref_values = [0, 0, 0, 0, 0, 0]

            for index, item in enumerate(ref_parameters):
                self.ref_values[index] = float(split_line[index + 2])
        else:
            self.bad_flag = True

    def calculate_planes(self, data: str, ref_plane_i: list, ref_values: list) -> float:
        """Finds the results from the MPLA command in the .lst file

        Converts into fractional coordinates using the unit cell from the

        above function

        Calculates the angle between it and the reference plane

        Note that the requires a number of calculations to be peformed on the unit cell. Some of these are planned to be moved to the unit cell class in the future.

        The unit cell class may then also be moved to the system_files folder.

        Args:
            data (str): .lst file with MPLA info as a string
            ref_plane_i (lst): the reference plane for comparison

        Returns:
            angle (float): the resulting angle from the calculations
        """

        index = 3
        flag = 0

        angle = 0

        for line in data:
            if "Least-squares planes" in line:
                plane = data[index]
                flag = 1
            index += 1

        data = []

        # signs do not matter = just make sure angle closest to the 0 instead of 180
        index = 0

        if flag == 1:
            for item in plane.split():
                try:
                    float(item)
                except ValueError:
                    pass
                else:
                    data.append(float(item))

            # calculate G the metric matrix
            alpha = np.radians(ref_values[3])
            beta = np.radians(ref_values[4])
            gamma = np.radians(ref_values[5])
            a = ref_values[0]
            b = ref_values[1]
            c = ref_values[2]
            G = np.zeros((3, 3))
            G[0, 0] = a**2
            G[0, 1] = a * b * np.cos(gamma)
            G[0, 2] = a * c * np.cos(beta)
            G[1, 0] = a * b * np.cos(gamma)
            G[1, 1] = b**2
            G[1, 2] = b * c * np.cos(alpha)
            G[2, 0] = a * c * np.cos(beta)
            G[2, 1] = b * c * np.cos(alpha)
            G[2, 2] = c**2
            detG = np.linalg.det(G)
            V = np.sqrt(detG)
            V_star = 1 / V
            a_star = b * c * np.sin(alpha) * V_star
            b_star = a * c * np.sin(beta) * V_star
            c_star = a * b * np.sin(gamma) * V_star
            alpha_star = np.arccos(
                (np.cos(beta) * np.cos(gamma) - np.cos(alpha))
                / (np.sin(beta) * np.sin(gamma))
            )
            beta_star = np.arccos(
                (np.cos(alpha) * np.cos(gamma) - np.cos(beta))
                / (np.sin(alpha) * np.sin(gamma))
            )
            gamma_star = np.arccos(
                (np.cos(alpha) * np.cos(beta) - np.cos(gamma))
                / (np.sin(alpha) * np.sin(beta))
            )
            # cacluates the orthonormalisation matrix M
            M = np.zeros((3, 3))
            M[0, 0] = a
            M[0, 1] = b * np.cos(gamma)
            M[0, 2] = c * np.cos(beta)
            M[1, 0] = 0
            M[1, 1] = b * np.sin(gamma)
            M[1, 2] = -c * np.sin(beta) * np.cos(alpha_star)
            M[2, 0] = 0
            M[2, 1] = 0
            M[2, 2] = c * np.sin(beta) * np.sin(alpha_star)
            M_star = np.linalg.inv(M)

            ## converst the molecule plane into fractional coordinates
            cart_coords = np.zeros((1, 3))
            cart_coords[0, 0] = data[0]
            cart_coords[0, 1] = data[1]
            cart_coords[0, 2] = data[2]
            frac_coords = np.dot(cart_coords, M_star)
            molecule_plane = frac_coords

            # convert the reference plane into fractional coordinates

            ref_plane = np.zeros((1, 3))
            ref_plane[0, 0] = ref_plane_i[0]
            ref_plane[0, 1] = ref_plane_i[1]
            ref_plane[0, 2] = ref_plane_i[2]
            ref_frac_coords = np.dot(ref_plane, M_star)

            ## Calculates the difference in angle between atom plane and reference plane

            angle = np.degrees(
                np.arccos(
                    np.dot(molecule_plane, ref_frac_coords.T)
                    / (
                        np.dot(
                            (np.linalg.norm(molecule_plane)),
                            (np.linalg.norm(ref_frac_coords)),
                        )
                    )
                )
            )
            # there is probably a better way to do this it seems to be required for the plotting function in pipeline to get a single value without double square brackets around it!
            angle = angle[0]
            angle = angle[0]
            if 180 - angle < 90:
                angle = 180 - angle
        return angle

    def find_planes(self, file_name: str) -> float:
        """Imports .lst file and then runs the function to calculate the rotation angle

        Args:
            file_name (str): full path to the .lst file with MPLA info

        Returns:
            angle (float): the resulting angle from the calculations
        """

        # Finds the calculation in the .lst file

        with open(file_name, "rt") as lst_file:
            data = lst_file.readlines()

        angle = self.calculate_planes(data, self.ref_plane, self.ref_values)

        return angle

    def analysis(self, lst_name: str, structure_number: int, results_path: str) -> None:
        """Runs the previous functions and also outputs the results to a .csv file
        Args:
            lst_name (str): full path to the .lst file for analysis
            structure_number (int): gives the structure number as independent variable
            results_path (str): full path to the location of the output results
        """

        self.bad_flag = False

        if lst_name == "":
            logging.info(__name__ + " : Refinement failed, no mean plane to analyse")
        else:
            self.grab_cell(pathlib.Path(lst_name))

            if self.bad_flag == False:
                rot_angle = self.find_planes(pathlib.Path(lst_name))
                self.df = pd.DataFrame(
                    {"Structure": [structure_number], "Rotation Angle": [rot_angle]}
                )
                os.chdir(results_path)
                try:
                    old_data = pd.read_csv("rotation_angles.csv")
                except FileNotFoundError:
                    self.df.to_csv("rotation_angles.csv", index=None)
                else:
                    new_df = old_data.append(self.df)
                    new_df.to_csv("rotation_angles.csv", index=None)
