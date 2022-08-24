#!/usr/bin/env python3

###################################################################################################
# -------------------------------------CX-ASAP: XPREP_pipeline-------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config, Directory_Browse
from data_reduction.modules.xprep_module import XPREP
import logging

# ----------Class Definition----------#


class XPREP_Pipeline:
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

    def multiple_xprep(
        self,
        matrix: str,
        location: str,
        file_name: str,
        space_group: str,
        chemical_formula: str,
    ) -> None:

        """This function runs xprep for a range of structures

        There should be a folder which contains folders of different data sets

        Each data set should have a reflection file in it for analysis

        Args:
            matrix (str): the matrix transformation to be input
            location (str): full path to the folder containing separate folders
                            each with a reflection file to run through xprep
            file_name (str): name of the file to be run through xprep
            space_group (str): space group to be input
            chemical_formula (str): chemical formula to be input
        """

        # Browses through the directories and runs xprep in each

        self.tree = Directory_Browse(location)
        self.xprep = XPREP()
        for index, item in enumerate(self.tree.directories):

            # Files needed for eventual checking haven't been created yet, so '.ins' is just a dummy variable here - it doesn't really do anything

            self.tree.enter_directory(item, ".ins")
            self.xprep.run(
                matrix,
                item,
                file_name,
                space_group,
                chemical_formula,
                "structure_" + str(index + 1),
            )

            # Sets the item_file to the new ins file to make sure that xprep has worked

            self.tree.item_file = "shelx.ins"
            self.tree.check_file_contents()
            self.tree.exit_directory()

    def multiple_asdefaults(self, location: str) -> None:

        """This function runs brute xprep for a range of structures

        Brute xprep means that it accepts all the defaults

        There should be a folder which contains folders of different data sets

        Each data set should have a reflection file in it for analysis

        Args:
            location (str): full path to the folder containing separate folders
        """

        self.tree = Directory_Browse(location)
        self.xprep = XPREP()
        for index, item in enumerate(self.tree.directories):
            self.tree.enter_directory(item, ".HKL_p1", ignore_check=True)
            if self.tree.item_file != "":
                self.xprep.asdefaults(item)
            # Sets the item_file to the new ins file to make sure that xprep has worked
            self.tree.exit_directory()
