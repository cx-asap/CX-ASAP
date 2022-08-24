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

# ----------Class Definition----------#


class SHELXT_Pipeline_auto:
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

    def new_shelxt_directory(self, location: str) -> None:
        """Makes a new directory for shelxt to be performed in
        Args:
        location (str): full path to the location containing folders of
                        datasets for shelxt to be applied
        """
        parent_location = pathlib.Path(location).parent
        os.chdir(parent_location)
        try:
            os.mkdir(pathlib.Path(location).stem + "_Shelxt")
        except FileExistsError:
            pass
        new_location = os.getcwd() + "/" + str(pathlib.Path(location).stem + "_Shelxt")
        os.chdir(location)
        self.tree = Directory_Browse(location, self.test_mode)
        files_to_copy = [".ins", ".hkl", ".cif_od", ".cif"]
        for item in self.tree.directories:
            new_folder = new_location + "/" + str(item.stem)
            try:
                os.mkdir(new_folder)
            except FileExistsError:
                pass
            for files in files_to_copy:
                self.tree.enter_directory(item, files)
                if self.tree.item_file != "":
                    shutil.copyfile(
                        self.tree.item_file, new_folder + "/" + self.tree.item_file.name
                    )
                self.tree.exit_directory()
        self.new_location = new_location

    def multi_shelxt(self, location: str, flags: bool = False) -> None:
        """Runs SHELXT on a seires of datasets in individual folders contained

        within a single parent folders

        Each dataset folder requires a .ins and a .hkl

        Args:
            location (str): full path to the file location
            flags (bool): additional flags to add into shelxt
                            NOT CURRENTLY IMPLEMENTED - FUTURE FRAMEWORK
        """

        # Browses through the directories and runs shelxt in each
        self.tree = Directory_Browse(location, self.test_mode)
        self.shelxt = SHELXT()
        for item in self.tree.directories:
            self.tree.enter_directory(item, ".ins", ignore_check=True)
            if self.tree.item_file != "":
                self.shelxt.run_shelxt(item, self.tree.item_file.stem, flags)

        # note this function assumes only the _a files are correct and removes the others including original
            for j in os.listdir():
                if j.endswith("_a.hkl") or j.endswith("_a.res") or j.endswith(".cif") or j.endswith(".cif_od"):
                    pass
                else:
                    os.remove(j)
                if j.endswith("_a.res"):
                    os.rename(j, j.strip("_a.res") + ".ins")
                if j.endswith("_a.hkl"):
                    os.rename(j, j.strip("_a.hkl") + ".hkl")

            self.tree.exit_directory()


# --------------------#


# def main():
# pass


# --------------------#


# if __name__ == "__main__":
# main()
