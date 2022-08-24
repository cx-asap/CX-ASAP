#!/usr/bin/env python3

###################################################################################################
# ---------------------------------CX-ASAP: XDS_Cell_Transformation--------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

import numpy as np
import logging
from system_files.utils import Nice_YAML_Dumper, Config

# ----------Class Definition----------#


class XDS_Cell_Transformation:
    def __init__(self) -> None:

        """Initialises the class

        Sets up the yaml parameters input by the user

        Also defines the location for the system yaml file

        which stores a yaml of code-only parameters accessible throughout

        the software package
        """

        # Setup yaml files and logger

        config = Config()

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    def find_axis_coordinates(self, file_path: str) -> "np.array":

        """Extracts the coordiantes for the unit cell axes from XDS_ASCII file

        Args:
            file_path (str): Full path to the XDS_ASCII file to extract the axis coordinates

        Returns:
            coordinates(np.array): returns the unit cell axes coordinates as an array
        """

        with open(file_path, "rt") as f:
            good_lines = {
                "UNIT_CELL_A-AXIS": [],
                "UNIT_CELL_B-AXIS": [],
                "UNIT_CELL_C-AXIS": [],
            }
            for lines in f:
                for param in good_lines:
                    if param in lines:
                        split = lines.strip("\n").split(" ")
                        for item in split:
                            try:
                                value = float(item)
                            except ValueError:
                                pass
                            else:
                                good_lines[param].append(value)

            coordinates = np.array(
                [
                    good_lines["UNIT_CELL_A-AXIS"],
                    good_lines["UNIT_CELL_B-AXIS"],
                    good_lines["UNIT_CELL_C-AXIS"],
                ]
            )

        return coordinates

    def find_transformation(self, file_1: str, file_2: str) -> str:

        """Calculates the transformation matrix between two XDS_ASCII files

        Args:
            file_1 (str): Full path to the first XDS_ASCII file
            file_2 (str): Full path to the second XDS_ASCII file
        Returns:
            M_xprep (str): transformation matrix expressed as a string
                            suitable for input into xprep
        """

        first_coordinates = self.find_axis_coordinates(file_1)
        second_coordinates = self.find_axis_coordinates(file_2)

        M = np.linalg.solve(
            first_coordinates.transpose(), second_coordinates.transpose()
        ).transpose()
        M = np.round(M, 1)
        M_xprep = ""
        for i in M:
            for j in i:
                M_xprep += str(j) + " "
        M_xprep.strip(" ")

        return M_xprep
