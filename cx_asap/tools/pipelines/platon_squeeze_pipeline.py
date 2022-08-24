#!/usr/bin/env python3

###################################################################################################
# ----------------------------------------CX-ASAP: TEMPLATE----------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config, Directory_Browse

# from system_files.yaml_configuration import Nice_YAML_Dumper, Config
# from system_files.directory_browse import Directory_Browse
from tools.modules.platon_squeeze import Platon_Squeeze
import os
import logging
import pathlib
import shutil

# ----------Class Definition----------#


class Squeeze_Pipeline:
    def __init__(self, test_mode: bool = False) -> None:

        """Initialises the class

        Sets up the yaml parameters input by the user

        Also defines the location for the system yaml file

        which stores a yaml of code-only parameters accessible throughout

        the software package

        Args:
            test_mode (bool): Automatically false, if true it will

        """

        # Setup yaml files and logger

        self.test_mode = test_mode

        config = Config(self.test_mode)

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    def new_squeeze_directory(self, location: str) -> None:

        """Makes a new directory for squeeze to be performed in

        Args:
            location (str): full path to the location containing folders of
                            datasets for squeeze to be applied

        """

        parent_location = pathlib.Path(location).parent

        os.chdir(parent_location)
        try:
            os.mkdir(pathlib.Path(location).stem + "_Squeeze")
        except FileExistsError:
            pass
        new_location = os.getcwd() + "/" + str(pathlib.Path(location).stem + "_Squeeze")

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

    # AKA SQUEEZY

    def multi_squeeze(self, location: str) -> None:

        """Runs the PLATON Squeeze function on a series of datasets

        contained in individual folders in a common parent location

        Args:
            location (str): full path to the location containing folders of
                            datasets for squeeze to be applied

        """

        self.tree = Directory_Browse(location, self.test_mode)
        self.squeeze = Platon_Squeeze()
        for item in self.tree.directories:
            self.tree.enter_directory(item, ".ins")
            if self.tree.item_file != "":
                self.squeeze.run_squeeze(self.tree.item_file.name)
                os.remove(self.tree.item_file)
                os.remove(self.tree.item_file.stem + ".hkl")
                try:
                    os.remove(self.tree.item_file.stem + ".cif")
                except:
                    pass
                try:
                    os.rename(
                        self.tree.item_file.stem + ".cif_od",
                        self.tree.item_file.stem + "_sq.cif_od",
                    )
                    # os.remove(self.tree.item_file.stem + '.cif_od')
                except:
                    pass

            crystals_squeeze_yuck_files = ["sqd.hkl", "sqd.ins"]

            # Might need to rename 'sqd.sqf' and 'sqd.sqz' to make file names consistent

            for j in os.listdir():
                for k in crystals_squeeze_yuck_files:
                    if j.endswith(k):
                        os.remove(j)

            self.tree.exit_directory()


# --------------------#


# def main():
# pass


# --------------------#

# if __name__ == "__main__":
# main()
