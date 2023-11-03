#!/usr/bin/env python3

###################################################################################################
# ----------------------------------------CX-ASAP: SHELXT----------------------------------------#
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
import signal
import sys

# ----------Class Definition----------#


class SHELXT:
    def __init__(self) -> None:
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

        config = Config()

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    def run_shelxt(self, location: str, file_name: str, flags: bool = False) -> None:
        """Runs shelxt on a specific .ins file with a .hkl in the same folder

        Args:
            location (str): full path to the file location
            file_name (str): name of file to run through shelxt
            flags (bool): additional flags to add into shelxt
            NOT CURRENTLY IMPLEMENTED - FUTURE FRAMEWORK
        """

        os.chdir(location)
        if flags == False:
            shelxt = subprocess.Popen(
                ["shelxt", file_name], stdin=subprocess.PIPE, encoding="utf8"
            )
            shelxt.stdin.close()
            shelxt.wait()
        else:
            shelxt = subprocess.Popen(
                ["shelxt", file_name, flags], stdin=subprocess.PIPE, encoding="utf8"
            )
            shelxt.stdin.close()
            shelxt.wait()

    def run_shelxt_as(self, location: str, file_name: str, flags: bool = False) -> None:
        """Runs shelxt on a specific .ins file with a .hkl in the same folder
        this has a built in time out for Brute use at the Australian Synchrotron beamliine
        it has a hard coded time out of two minutes, but this can be changed in the future
        Args:
            location (str): full path to the file location
            file_name (str): name of file to run through shelxt
            flags (bool): additional flags to add into shelxt
                            NOT CURRENTLY IMPLEMENTED - FUTURE FRAMEWORK
        """
        os.chdir(location)
        timeout_s = 0.5
        if flags == False:
            try:
                shelxt = subprocess.Popen(
                    ["shelxt", file_name],
                    stdin=subprocess.PIPE,
                    encoding="utf8",
                    start_new_session=True,
                )
                shelxt.stdin.close()
                shelxt.wait(timeout=timeout_s)
            except subprocess.TimeoutExpired:
                print("timeout worked :)")
                os.killpg(os.getpgid(shelxt.pid), signal.SIGINT)
        else:
            shelxt = subprocess.Popen(
                ["shelxt", file_name, flags], stdin=subprocess.PIPE, encoding="utf8"
            )
            shelxt.stdin.close()
            shelxt.wait()
