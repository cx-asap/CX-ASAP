#!/usr/bin/env python3

###################################################################################################
# --------------------------------------CX-ASAP: cif_combine---------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config, File_Sorter
from cif_validation.modules.cif_merge import Cif_Merge
from CifFile import ReadCif
import os
import logging

# ----------Class Definition----------#


class Cif_Combine:
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

    def combine_cifs_single_folder(self, location: str) -> None:
        """Combines all the CIF files that are found within a single folder

        Args:
            location (str): Full path to the folder with CIFs for combining
        """

        os.chdir(location)
        validate = Cif_Merge()
        sorter = File_Sorter()

        for item in sorter.sorted_properly(os.listdir()):
            if item.endswith(".cif"):
                self.cif = ReadCif(item)
                with open("combined.cif", "a") as f:
                    f.write(self.cif.WriteOut())
                validate.validate_CIFs(item)
                with open("combined_check_CIF.chk", "a") as f:
                    for line in validate.validation:
                        f.write(line)
