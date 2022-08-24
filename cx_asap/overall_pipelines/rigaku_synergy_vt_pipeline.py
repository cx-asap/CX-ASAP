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
from data_refinement.pipelines.refine_pipeline import Refinement_Pipeline
from cif_validation.pipelines.cif_pipeline import CIF_Compile_Pipeline
from post_refinement_analysis.pipelines.variable_temperature_analysis import (
    VT_Analysis_Pipeline,
)
import os
import pathlib
import shutil
import logging

# ----------Class Definition----------#


class Synergy_VT:
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

    def initialise(self, experiment_location) -> None:

        """Makes two directories:

        A directory for refinement graphs and result report

        A directory for compiled CIFs and CIF analysis

        Args:
            experiment_location (str): full path to the folder containing folders
                                        of data for analysis
        """

        os.chdir(experiment_location)

        self.stats_location = (
            pathlib.Path(experiment_location) / "Refinement_Statistics"
        )
        self.results_location = pathlib.Path(experiment_location) / "CIF_Analysis"

        if os.path.exists("Refinement_Statistics") != True:
            os.mkdir("Refinement_Statistics")

        os.chdir(self.stats_location)
        tmp = os.listdir()
        tmp2 = []
        for item in tmp:
            if os.path.isdir(item):
                tmp2.append(int(item))
        tmp2.sort()

        try:
            self.new_folder = int(tmp2[-1]) + 1
        except IndexError:
            self.new_folder = 1

        individual_stats_location = pathlib.Path(self.stats_location) / str(
            self.new_folder
        )

        os.mkdir(individual_stats_location)

        self.stats_location = individual_stats_location

        os.chdir(experiment_location)

        if os.path.exists("CIF_Analysis") != True:
            os.mkdir("CIF_Analysis")

        os.chdir(self.results_location)
        tmp = os.listdir()
        tmp2 = []
        for item in tmp:
            if os.path.isdir(item):
                tmp2.append(int(item))
        tmp2.sort()

        try:
            self.new_folder = int(tmp2[-1]) + 1
        except IndexError:
            self.new_folder = 1

        individual_results_location = pathlib.Path(self.results_location) / str(
            self.new_folder
        )

        os.mkdir(individual_results_location)

        self.results_location = individual_results_location

        os.chdir(experiment_location)

    def process(
        self,
        location: str,
        reference: str,
        graph_output_location: str,
        refinements_to_check: int,
        tolerance: float,
        max_cycles: int,
    ):

        """Runs SHELXL on a series of structures based on one reference

        Args:
            location (str): full path to the folder containing folders of .ins files
            reference (str): full path to the reference .ins/.res file
            graph_output_location (str): full path to the location of output files
            refinements_to_check (int): number of refinements to check for shift convergence
            tolerance (float): target shift value
            max_cycles (int): maximum cycles SHELXL can run before stopping
        """

        shelxl = Refinement_Pipeline()
        shelxl.multiple_refinement(
            location,
            reference,
            graph_output_location,
            refinements_to_check,
            tolerance,
            max_cycles,
        )

    def analyse(
        self,
        reference: str,
        data_location: str,
        results_location: str,
        solution_program: str,
        chemical_formula: str,
        crystal_habit: str,
        crystal_colour: str,
        max_dimension: str,
        middle_dimension: str,
        min_dimension: str,
        bonds: bool,
        angles: bool,
        torsions: bool,
        cif_parameters: list,
        atoms_for_analysis: list,
        adps: bool,
        instrument_ending: str = False,
        instrument_file: str = False,
    ) -> None:

        """Compiles all of the output CIFs and runs an analysis pipeline on these files

        Args:
            reference (str): full path to the reference .ins/.res file
            data_location (str): full path to the folder containing folders of .ins files
            results_location (str): full path to the folder that will contain results
            solution_program(str): the structure solution program used
            chemical_formula(str): the chemical formula
            crystal_habit(str): the crystal habit
            crystal_colour(str): the colour of the crystal
            max_dimension(str): the largest dimension of the crystal
            middle_dimension(str): the middle dimension of the crystal
            min_dimension(str): the smallest dimension of the crystal
            bonds (bool): whether or not bond analysis should be run
            angles (bool): whether or not angle analysis should be run
            torsions (bool): whether or not torsion analysis should be run
            cif_parameters (list): list of parameters in the CIF for tabulation
            atoms_for_analysis (list): list of atoms for geometrical analysis
            adps (bool): whether or not ADP analysis should be run
            instrument_ending(str): ending of the instrument file if they are
                                    all consistent with different names
            instrument_file(str): name of the instrument file if they are all consistent

        """

        cif = CIF_Compile_Pipeline()
        cif.configure(
            data_location,
            solution_program,
            chemical_formula,
            crystal_habit,
            crystal_colour,
            max_dimension,
            middle_dimension,
            min_dimension,
            instrument_ending,
            instrument_file,
        )
        cif.compile_cifs(results_location, [self.stats_location, self.results_location])
        analysis = VT_Analysis_Pipeline()
        analysis.analyse_data(
            reference,
            results_location,
            cif_parameters,
            atoms_for_analysis,
            bonds,
            angles,
            torsions,
            adps,
        )
