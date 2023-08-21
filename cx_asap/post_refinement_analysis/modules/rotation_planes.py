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

            ref_parameters = ["ref_INS_a", "ref_INS_b", "ref_INS_c"]

            self.ref_values = [0, 0, 0]

            for index, item in enumerate(ref_parameters):
                self.ref_values[index] = float(split_line[index + 2])
        else:
            self.bad_flag = True

    def calculate_planes(self, data: str, ref_plane: list, ref_values: list) -> float:

        """Finds the results from the MPLA command in the .lst file

        Converts into fractional coordinates using the unit cell from the

        above function

        Calculates the angle between it and the reference plane

        Args:
            data (str): .lst file with MPLA info as a string
            ref_plane (lst): the reference plane for comparison

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

            # Calculates the difference in angle between atom plane and reference plane

            x = data[0] 
            y = data[1] 
            z = data[2]

            molecule_plane = [x, y, z]

            dot_product = (
                molecule_plane[0] * ref_plane[0]
                + molecule_plane[1] * ref_plane[1]
                + molecule_plane[2] * ref_plane[2]
            )

            normal_molecule = math.sqrt(
                (molecule_plane[0] ** 2)
                + (molecule_plane[1] ** 2)
                + (molecule_plane[2] ** 2)
            )

            normal_reference = math.sqrt(
                (ref_plane[0] ** 2) + (ref_plane[1] ** 2) + (ref_plane[2] ** 2)
            )

            angle = math.acos(dot_product / (normal_molecule * normal_reference)) * (
                180 / math.pi
            )

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
