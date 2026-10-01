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
from tools.modules.platon_twinrotmat import Platon_Twin
import os
import logging
import pathlib
import shutil

# ----------Class Definition----------#


class Twin_Pipeline:
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

    def new_twinrotmat_directory(self, location: str) -> None:
        """Makes a new directory for twinrotmat to be performed in

        Args:
            location (str): full path to the location containing folders of
                            datasets for twinrotmat to be applied

        """

        parent_location = pathlib.Path(location).parent

        os.chdir(parent_location)
        try:
            os.mkdir(pathlib.Path(location).stem + "_Twinrotmat")
        except FileExistsError:
            pass
        new_location = (
            os.getcwd() + "/" + str(pathlib.Path(location).stem + "_Twinrotmat")
        )

        os.chdir(location)

        self.tree = Directory_Browse(location, self.test_mode)

        # files_to_copy = [".ins", ".hkl", ".cif_od", ".cif"]
        files_to_copy = [".cif", ".fcf", ".cif_od"]
        # check this - it will need the .fcf

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

    # AKA TWINNY

    def multi_twinrotmat(self, location: str) -> None:
        """Runs the PLATON Twinrotmat function on a series of datasets

        contained in individual folders in a common parent location

        Args:
            location (str): full path to the location containing folders of
                            datasets for twinrotmat to be applied

        """

        self.tree = Directory_Browse(location, self.test_mode)
        self.twinrotmat = Platon_Twin()
        for item in self.tree.directories:
            self.tree.enter_directory(item, ".cif")
            if self.tree.item_file != "":
                self.twinrotmat.run_twinrotmat(self.tree.item_file.name)
                os.remove(self.tree.item_file)
                os.remove(self.tree.item_file.stem + "_pl.spf")
                os.remove(self.tree.item_file.stem + ".lis")
                os.remove(self.tree.item_file.stem + ".fcf")
                os.remove(self.tree.item_file.stem + ".ckf")

                try:
                    os.rename(
                        self.tree.item_file.stem + ".cif_od",
                        self.tree.item_file.stem + "_tw.cif_od",
                    )
                except:
                    pass

            self.tree.exit_directory()


# --------------------#


# def main():
# pass


# --------------------#

# if __name__ == "__main__":
# main()
