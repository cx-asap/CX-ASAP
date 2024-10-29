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
import os
import logging
import pathlib
import shutil

# ----------Class Definition----------#


class XPREP_Pipeline:
    def __init__(self, test_mode: bool = False) -> None:
        """Initialises the class

        Sets up the yaml parameters input by the user

        Also defines the location for the system yaml file

        which stores a yaml of code-only parameters accessible throughout

        the software package
        """

        # Setup yaml files and logger

        self.test_mode = test_mode
        config = Config(self.test_mode)

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    def new_xprep_directory(self, location: str) -> None:
        """Makes a new directory for xprep transformation to be performed in
        Args:
            location (str): full path to the location containing folders of
                            datasets for xprep to be applied
        """
        parent_location = pathlib.Path(location).parent
        os.chdir(parent_location)
        try:
            os.mkdir(pathlib.Path(location).stem + "_Xprep")
        except FileExistsError:
            pass
        new_location = os.getcwd() + "/" + str(pathlib.Path(location).stem + "_Xprep")

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

    def multiple_asdefaults(
        self, location: str, formula: str = "C40H30N6FeCl2O8"
    ) -> None:
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
                self.xprep.asdefaults(item, formula)
            # Sets the item_file to the new ins file to make sure that xprep has worked
            self.tree.exit_directory()

    def multiple_xprep_custom(
        self,
        location: str,
        data_type: str,
        lattice_type: str,
        matrix: str,
        file_name: str,
        space_group: str,
        chemical_formula: str,
        output_name: str,
    ) -> None:
        """This function runs the transformation function of xprep for a range of structures
        There should be a folder which contains folders of different data sets
        Each data set should have a reflection file in it for analysis
        Args:
           data_type (str): the type of data to be input (e.g. XDS, Shelx etc)
            lattice_type (str): the type of lattice to be input (e.g. P, C, I etc)
            matrix (str): the matrix transformation to be input
            location (str): full path to the folder of the file to be run
            file_name (str): name of the file to be run through xprep
            space_group (str): space group to be input
            chemical_formula (str): chemical formula to be input
            output_name (str): name of the output .hkl/.ins file from xprep
        """
        # Browses through the directories and runs xprep in each
        self.tree = Directory_Browse(location)
        self.xprep = XPREP()
        for index, item in enumerate(self.tree.directories):
            # Files needed for eventual checking haven't been created yet, so '.ins' is just a dummy variable here - it doesn't really do anything
            self.tree.enter_directory(item, ".ins")
            # converts ins to p4p by removing wavelength from CELL line
            if self.tree.item_file != "":
                with open(self.tree.item_file, "r") as f:
                    lines = f.readlines()
                with open(self.tree.item_file.stem + ".p4p", "w") as f:
                    for line in lines:
                        if line.startswith("CELL"):
                            parts = line.split()
                            parts.pop(1)
                            line = " ".join(parts) + "\n"
                        f.write(line)
                    else:
                        pass

                self.xprep.xprep_custom(
                    location,
                    data_type,
                    lattice_type,
                    matrix,
                    self.tree.item_file.stem,
                    space_group,
                    chemical_formula,
                    output_name + str(index + 1),
                )
                # tidy up
            # try:
            #     os.remove(self.tree.item_file.stem + ".p4p")
            # except:
            #     pass
            try:
                os.remove(self.tree.item_file.stem + ".ins")
            except:
                pass
            try:
                os.remove(self.tree.item_file.stem + ".hkl")
            except:
                pass
            try:
                os.remove(self.tree.item_file.stem + ".cif")
            except:
                pass
            try:
                os.rename(
                    self.tree.item_file.stem + ".cif_od",
                    output_name + str(index + 1) + ".cif_od",
                )
            except:
                pass

            # Sets the item_file to the new ins file to make sure that xprep has worked
            # self.tree.item_file = "shelx.ins"
            # self.tree.check_file_contents()
            self.tree.exit_directory()
