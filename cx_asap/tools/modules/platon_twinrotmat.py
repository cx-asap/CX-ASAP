#!/usr/bin/env python3

###################################################################################################
# ----------------------------------------CX-ASAP: TEMPLATE----------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config
import logging
import subprocess
import pathlib
import os

# ----------Class Definition----------#


class Platon_Twin:
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

    def run_twinrotmat(self, file_name: str) -> None:
        """Runs the PLATON twinrotmat function on a specified cif (and accomanying List 4 or list 8 .fcf) they must have the same file name stem
        you will be required to click on hklf5 generate if you want one...

        Args:
            file_name (str): full path to the cif file to run through PLATON
        """

        file_name_path = pathlib.Path(file_name)

        if file_name_path.name != file_name:
            os.chdir(file_name_path.parent)
            twinrotmat = subprocess.Popen(
                ["platon", "-T1", file_name_path.name],
                stdin=subprocess.PIPE,
                encoding="utf8",
            )
            twinrotmat.stdin.close()
            twinrotmat.wait()
        else:
            twinrotmat = subprocess.Popen(
                ["platon", "-T1", file_name], stdin=subprocess.PIPE, encoding="utf8"
            )
            twinrotmat.stdin.close()
            twinrotmat.wait()
