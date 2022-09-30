#!/usr/bin/env python3

###################################################################################################
# ---------------------------------------CX-ASAP: XPREP_module-------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config
import os
import subprocess
import logging

# ----------Class Definition----------#


class XPREP:
    def __init__(self) -> None:

        """Initialises the class

        Sets up the yaml parameters input by the user

        Also defines the location for the system yaml file

        which stores a yaml of code-only parameters accessible throughout

        the software package
        """

        # Set up yaml files and logger

        config = Config()

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    def run(
        self,
        matrix: str,
        location: str,
        file_name: str,
        space_group: str,
        chemical_formula: str,
        output_name: str,
    ) -> None:

        """This function runs xprep for a single structure

        Args:
            matrix (str): the matrix transformation to be input
            location (str): full path to the folder of the file to be run
            file_name (str): name of the file to be run through xprep
            space_group (str): space group to be input
            chemical_formula (str): chemical formula to be input
            output_name (str): name of the output .hkl/.ins file from xprep
        """

        os.chdir(location)

        xprep = subprocess.Popen(["xprep"], stdin=subprocess.PIPE, encoding="utf8")
        xprep.stdin.write(file_name + "\n")
        xprep.stdin.write("X\n")
        xprep.stdin.write("Y\n")
        xprep.stdin.write("P\n")
        xprep.stdin.write("U\n")
        xprep.stdin.write("O\n")
        xprep.stdin.write(matrix + "\n")
        xprep.stdin.write("S\n")
        xprep.stdin.write("I\n")
        xprep.stdin.write(space_group + "\n")
        xprep.stdin.write("Y\n")
        xprep.stdin.write("D\n")
        xprep.stdin.write("S\n")
        xprep.stdin.write("A\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("E\n")
        xprep.stdin.write("C\n")
        xprep.stdin.write(chemical_formula + "\n")
        xprep.stdin.write("E\n")
        xprep.stdin.write("F\n")
        xprep.stdin.write(output_name + "\n")
        xprep.stdin.write("S\n")
        xprep.stdin.write("Y\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("Q\n")
        xprep.stdin.close()
        xprep.wait(20)

        # Below code for windows... when add into main release need a better solution for closing xprep

        # try:
        # xprep.wait(120)
        # except subprocess.TimeoutExpired:
        # xprep.terminate()

    def asdefaults(self, location: str, formula: str = "C40H30N6FeCl2O8") -> None:

        """This function runs xprep for a single aussynchrotron

        structure accepting all defaults (AS has default name XDS_ASCII.HKL_p1)

        Args:
            location (str): full path to the folder of the file to be run
            formula (str): chemical formula
        """
        os.chdir(location)
        xprep = subprocess.Popen(["xprep"], stdin=subprocess.PIPE, encoding="utf8")
        xprep.stdin.write("XDS_ASCII.HKL_p1\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write(formula + "\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("cxasap\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("\n")
        xprep.stdin.write("Q\n")
        xprep.stdin.close()
        xprep.wait(20)

        # Below code for windows... when add into main release need a better solution for closing xprep

        ###try:
        ###xprep.wait(15)
        ###except subprocess.TimeoutExpired:
        ###xprep.terminate()
