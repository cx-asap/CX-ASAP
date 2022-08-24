#!/usr/bin/env python3

###################################################################################################
# --------------------------------------CX-ASAP: xds_pipeline--------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config, Directory_Browse
from data_reduction.modules.xds_reprocess import XDS_reprocess
import logging

# ----------Class Definition----------#


class XDS_Pipeline:
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

    def multiple_xds(self, location: str, XDS_INP_path: str) -> None:

        """Browses through the directories in a given location

        Copyings a reference XDS.INP file into it and runs xds in each

        Args:
            location (str): path to the folder containing folders of XDS data sets
            XDS_INP_path (str): full path to the XDS.INP file
                                for copying prior to running xds
        """

        self.tree = Directory_Browse(location)
        self.xds = XDS_reprocess()
        for index, item in enumerate(self.tree.directories):

            # Files needed for eventual checking haven't been created yet, so '.ins' is just a dummy variable here - it doesn't really do anything

            self.tree.enter_directory(item, ".ins")
            self.xds.run_xds(item.stem, XDS_INP_path, item, str(index + 1))
            self.tree.exit_directory()
