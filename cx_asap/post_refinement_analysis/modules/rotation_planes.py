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
from system_files.crystal_math import orthonorm_matrix
from post_refinement_analysis.modules.lst_read import LST_Read
import pathlib
import os
import pandas as pd
import logging
import numpy as np
import re

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

        self.lst_reader = LST_Read(self.test_mode)

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

    def calculate_planes(self, data: list, ref_plane: list, ref_values: list) -> list:
        """Calculates the angle between the reference plane and every MPLA plane
        found in the .lst file.

        Args:
            data (list): lines of the .lst file
            ref_plane (list): reference crystallographic plane as [h, k, l]
            ref_values (list): unit cell parameters [a, b, c, alpha, beta, gamma]

        Returns:
            angles (list): angle in degrees between the reference plane and each
                           MPLA plane in order; empty list if no planes found
        """

        if ref_values:
            self.ref_values = ref_values

        if ref_plane:
            self.ref_plane = ref_plane

        plane_normals = self.lst_reader.extract_plane_normals(data)

        if not plane_normals:
            logging.warning(__name__ + " : No MPLA planes found in .lst file")
            return []

        M = orthonorm_matrix(self.ref_values)
        M_star = np.linalg.inv(M)

        # Convert reference plane vector to fractional space
        ref = np.array(
            [[self.ref_plane[0], self.ref_plane[1], self.ref_plane[2]]], dtype=float
        )
        ref_frac = np.dot(ref, M_star)

        angles = []
        for normal in plane_normals:
            cart = np.array([[normal[0], normal[1], normal[2]]], dtype=float)
            frac = np.dot(cart, M_star)

            angle = float(
                np.degrees(
                    np.arccos(
                        np.dot(frac, ref_frac.T)
                        / (np.linalg.norm(frac) * np.linalg.norm(ref_frac))
                    )
                )[0][0]
            )
            if 180 - angle < 90:
                angle = 180 - angle
            angles.append(angle)

        return angles

    def find_planes(self, file_name: str) -> list:
        """Reads a .lst file and calculates rotation angles for all MPLA planes.

        Args:
            file_name (str): full path to the .lst file with MPLA info

        Returns:
            angles (list): angle in degrees between the reference plane and each
                           MPLA plane in order
        """

        with open(file_name, "rt") as lst_file:
            data = lst_file.readlines()

        return self.calculate_planes(data, self.ref_plane, self.ref_values)

    def analysis(self, lst_name: str, structure_number: int, results_path: str) -> None:
        """Calculates rotation angles for all MPLA planes and appends to a .csv.

        One column per MPLA plane is written, named 'MPLA_1_Rotation_Angle',
        'MPLA_2_Rotation_Angle', etc.

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
                rot_angles = self.find_planes(pathlib.Path(lst_name))
                if not rot_angles:
                    logging.warning(
                        __name__
                        + " : No rotation angles calculated for "
                        + str(lst_name)
                    )
                    return

                row = {"Structure": [structure_number]}
                for i, angle in enumerate(rot_angles, start=1):
                    row[f"MPLA_{i}_Rotation_Angle"] = [angle]

                self.df = pd.DataFrame(row)
                os.chdir(results_path)
                try:
                    old_data = pd.read_csv("rotation_angles.csv")
                except FileNotFoundError:
                    self.df.to_csv("rotation_angles.csv", index=None)
                else:
                    new_df = pd.concat([old_data, self.df])
                    new_df.to_csv("rotation_angles.csv", index=None)

    def find_interplane_angle(self, file_name: str) -> "float | None":
        """Reads a .lst file and extracts the SHELXL inter-plane angle.

        Requires two MPLA commands in the .ins file so that SHELXL reports
        'Angle to previous plane' in the Least-squares planes section.

        Args:
            file_name (str): full path to the .lst file

        Returns:
            angle (float): angle in degrees between the two planes,
                           or None if not found
        """

        with open(file_name, "rt") as f:
            for line in f:
                if "Angle to previous plane" in line:
                    match = re.search(r"=\s*([\d.]+)", line)
                    if match:
                        return float(match.group(1))
        return None

    def analyse_interplane_angle(
        self, lst_name: str, structure_number: int, results_path: str
    ) -> None:
        """Extracts the SHELXL inter-plane angle and appends it to a .csv file.

        Requires two MPLA commands in the .ins file.

        Args:
            lst_name (str): full path to the .lst file with two MPLA commands
            structure_number (int): structure number as the independent variable
            results_path (str): full path to the output results directory
        """

        if lst_name == "":
            logging.info(__name__ + " : Refinement failed, no mean planes to analyse")
            return

        angle = self.find_interplane_angle(pathlib.Path(lst_name))
        if angle is None:
            logging.warning(
                __name__
                + " : Could not find 'Angle to previous plane' in "
                + str(lst_name)
                + ". This usually means fewer than two MPLA commands were present; "
                + "using fallback value 0.0."
            )
            angle = 0.0

        df = pd.DataFrame(
            {"Structure": [structure_number], "Interplane Angle": [angle]}
        )

        os.chdir(results_path)
        try:
            old_data = pd.read_csv("interplane_angles.csv")
        except FileNotFoundError:
            df.to_csv("interplane_angles.csv", index=None)
        else:
            new_df = pd.concat([old_data, df])
            new_df.to_csv("interplane_angles.csv", index=None)

