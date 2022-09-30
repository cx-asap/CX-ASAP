#!/usr/bin/env python3

###################################################################################################
# ----------------------------------------CX-ASAP: SHELXT----------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config, Directory_Browse

# from system_files.yaml_configuration import Nice_YAML_Dumper, Config
# from system_files.directory_browse import Directory_Browse
from tools.modules.shelx_t import SHELXT
import os
import logging
import pathlib
import shutil
import time

# ----------Class Definition----------#


class SHELXT_Pipeline:
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

    def multi_shelxt(self, location: str, file_name: str, flags: bool = False):

        """Runs SHELXT on a seires of datasets in individual folders contained

        within a single parent folders

        Each dataset folder requires a .ins and a .hkl

        Args:
            location (str): full path to the file location
            file_name (str): name of file to run through shelxt
            flags (bool): additional flags to add into shelxt
                            NOT CURRENTLY IMPLEMENTED - FUTURE FRAMEWORK
        """

        # Browses through the directories and runs shelxt in each
        self.tree = Directory_Browse(location, self.test_mode)
        self.shelxt = SHELXT()

        # This is because the first few datasets weren't running through shelxt - xprep was still open or something
        time.sleep(15)

        for item in self.tree.directories:
            self.tree.enter_directory(item, ".ins", ignore_check=True)
            if self.tree.item_file != "":
                self.shelxt.run_shelxt(item, self.tree.item_file.stem, flags)
            self.tree.exit_directory()


# --------------------#


# def main():
# pass


# --------------------#


# if __name__ == "__main__":
# main()
