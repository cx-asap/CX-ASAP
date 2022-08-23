#!/usr/bin/env python3

###################################################################################################
# ---------------------------------------CX-ASAP: CIF pipeline-------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

import os
import pathlib
from cif_validation.modules.cif_merge import Cif_Merge
import logging
from system_files.utils import Nice_YAML_Dumper, Config, Directory_Browse

# ----------Class Definition----------#


class CIF_Compile_Pipeline:
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

    def configure(
        self,
        location: str,
        solution_program: str,
        chemical_formula: str,
        crystal_habit: str,
        crystal_colour: str,
        max_dimension: str,
        middle_dimension: str,
        min_dimension: str,
        instrument_ending: str = False,
        instrument_file: str = False,
        additional_user_parameters: list = False,
    ) -> None:

        """Allows input of variables to be saved to the class

        Args:
            location(str): full path to the folder containing folders with CIF files in them
            solution_program(str): the structure solution program used
            chemical_formula(str): the chemical formula
            crystal_habit(str): the crystal habit
            crystal_colour(str): the colour of the crystal
            max_dimension(str): the largest dimension of the crystal
            middle_dimension(str): the middle dimension of the crystal
            min_dimension(str): the smallest dimension of the crystal
            instrument_ending(str): ending of the instrument file if they are
                                    all consistent with different names
            instrument_file(str): name of the instrument file if they are all consistent
            additional_user_parameters(list): list of extra cif parameters the user has edited
                                            CURRENTLY UNDER DEVELOPMENT AS OPTION NOT IMPLEMENTED
        """

        self.tree = Directory_Browse(location, self.test_mode)
        self.finalise = Cif_Merge(self.test_mode)
        self.location = location
        self.solution_program = solution_program
        self.chemical_formula = chemical_formula
        self.crystal_habit = crystal_habit
        self.crystal_colour = crystal_colour
        self.max_dimension = max_dimension
        self.middle_dimension = middle_dimension
        self.min_dimension = min_dimension
        self.instrument_ending = instrument_ending
        self.instrument_file = instrument_file
        self.additional_user_parameters = additional_user_parameters

    def compile_cifs(self, output_location: str, ignored_folders: list = []) -> None:

        """Goes to a defined experiment location (self.location)

        Iterates through all the folders present in this location

        Imports the instrument CIF and the structure CIF present in the folder

        Makes edits and merges them together

        Runs platon checkCIF and also outputs a checkCIF file

        Outputs checkCIF and merged CIF into the output_location

        Args:
            output_location(str): full path to the output location
            ignored_folders(list): list of any folders which should be ignored during the iteration
        """

        # Compiles the cif based on the directory browse for a set of data

        for item in self.tree.directories:

            if item not in ignored_folders:

                self.tree.enter_directory(
                    item, ".cif", self.instrument_file, ignored_folders
                )
                if self.instrument_file == False:
                    self.finalise.import_CIFs(
                        pathlib.Path(self.tree.item_name + self.instrument_ending),
                        self.tree.item_file,
                    )
                if (
                    self.instrument_ending == False
                    and self.instrument_file == "autoprocess.cif"
                ):
                    self.finalise.synchrotron_cif_edit(
                        pathlib.Path(self.instrument_file)
                    )
                    self.finalise.import_CIFs(
                        pathlib.Path(self.instrument_file), self.tree.item_file
                    )
                if (
                    self.instrument_ending == False
                    and self.instrument_file != "autoprocess.cif"
                ):
                    self.finalise.import_CIFs(
                        pathlib.Path(self.instrument_file), self.tree.item_file
                    )
                self.finalise.merge_CIFs()
                if self.additional_user_parameters == False:
                    self.finalise.user_edits(
                        self.solution_program,
                        self.chemical_formula,
                        self.crystal_habit,
                        self.crystal_colour,
                        self.max_dimension,
                        self.middle_dimension,
                        self.min_dimension,
                    )

                if self.tree.item_file != "":
                    logging.info(__name__ + " : Adding Cif..." + str(item))
                    self.finalise.single_write_out(self.tree.item_file.name)
                    self.finalise.validate_CIFs(self.tree.item_file.name)
                    self.finalise.write_out(
                        output_location,
                        "combined.cif",
                        "check_CIF.chk",
                        self.tree.item_file,
                    )
                    self.tree.check_file_contents()
                self.tree.exit_directory()

        os.chdir(output_location)
