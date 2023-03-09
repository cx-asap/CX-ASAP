#!/usr/bin/env python3

###################################################################################################
# -------------------------------------CX-ASAP: General Pipeline-----------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config, File_Sorter, Cell_Import
from data_refinement.pipelines.refine_pipeline import Refinement_Pipeline
from cif_validation.pipelines.cif_pipeline import CIF_Compile_Pipeline
from cif_validation.modules.instrument_cif_generation import Instrument_CIF
from post_refinement_analysis.pipelines.variable_cif_parameter import (
    Variable_Analysis_Pipeline,
)
import os
import pathlib
import shutil
import subprocess
from CifFile import ReadCif
from CifFile import CifFile
from CifFile import CifBlock
import logging

# ----------Class Definition----------#


class General_Pipeline:
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

        config = Config(test_mode)

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

        self.sorter = File_Sorter()

    def make_dirs(self, experiment_location: str) -> None:

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

    def reference_extract(
        self,
        cif_location: str,
        working_directory: str,
        varying_data: list,
        varying_param: str,
    ):

        """Runs shredCIF on a reference .cif

        Extracts the .ins/.res file to use as a reference structure

        Also makes an instrument CIF based on the reference cif

        Args:
            cif_location (str): full path to the reference .cif
            working_directory (str): full path to the folder containing folders of data
            varying_data (list): list of data for the changing param
                                EG if variable temperature, list all temperatures (with errors)
            varying_param (str): the parameter that is varying through the experiment
                                in correct CIF format
                                IE "_diffrn_ambient_temperature" for a variable temp experiment
        """

        cif_path = pathlib.Path(cif_location)

        os.chdir(cif_path.parent)

        shredCIF = subprocess.call(["shredcif", cif_path.name])

        for item in os.listdir():
            if item.endswith(".ins") or item.endswith(".res"):
                self.reference_res = pathlib.Path.cwd() / item

        try:
            test = self.reference_res
        except AttributeError:
            print("Error! Check Error Logs")
            logging.info(
                __name__
                + " : Could not extract ins/res from reference CIF using shredCIF"
            )
            exit()

        ref_cell = Cell_Import(self.test_mode)
        ref_cell.cell_import(self.reference_res)

        make_instrument_cif = Instrument_CIF(self.test_mode)
        make_instrument_cif.read_reference_cif(cif_location)
        make_instrument_cif.make_instrument_cif()

        instrument_path = pathlib.Path.cwd() / "instrument.cif"

        os.chdir(working_directory)

        # DON'T use enumerate here, because stats_location and results_location would also contribute to numbers

        index = 0

        for item in self.sorter.sorted_properly(os.listdir()):
            if (
                item != pathlib.Path(self.stats_location.parent).stem
                and item != pathlib.Path(self.results_location.parent).stem
                and os.path.isdir(item) == True
                and item != 'error_logs'
            ):

                print(item)

                shutil.copy(instrument_path, item)

                os.chdir(item)

                cif = ReadCif(instrument_path.name)
                data_block = cif.first_block()
                data_block[varying_param] = varying_data[index]

                with open("instrument.cif", "w") as f:
                    f.write(cif.WriteOut())

                index += 1
                os.chdir("..")

        self.chemical_formula = ""
        self.crystal_habit = ""
        self.crystal_colour = ""
        self.max_dimension = ""
        self.mid_dimension = ""
        self.min_dimension = ""

    def file_tree_setup(
        self,
        res_location: str,
        working_directory: str,
        varying_data: list,
        varying_param: str,
    ) -> None:

        """Similar to 'reference_extract', but assuming no reference .cif is available

        Instead, starts with a reference .ins/.res

        YAML will ask for instrument information so this function will make

        the instrument cif from that

        This function will also create the required file structure (ie one folder for

        each dataset in a common directory)

        It will then pause and prompt the user to put all their data in these newly

        created folders

        Args:
            res_location (str): full path to the reference .res
            working_directory (str): full path to the folder containing folders of data
            varying_data (list): list of data for the changing param
                                EG if variable temperature, list all temperatures (with errors)
            varying_param (str): the parameter that is varying through the experiment
                                in correct CIF format
                                IE "_diffrn_ambient_temperature" for a variable temp experiment
        """

        self.reference_res = pathlib.Path(res_location)
        ref_cell = Cell_Import(self.test_mode)
        ref_cell.cell_import(self.reference_res)

        os.chdir(self.reference_res.parent)

        instrument_cif = CifFile()
        instrument_data_block = CifBlock()
        instrument_cif["instrument_information"] = instrument_data_block

        for item in self.cfg:
            if item.startswith("_"):
                instrument_cif["instrument_information"][item] = self.cfg[item]

        with open("instrument.cif", "w") as f:
            f.write(instrument_cif.WriteOut())

        instrument_path = pathlib.Path.cwd() / "instrument.cif"

        os.chdir(working_directory)

        for i, item in enumerate(varying_data):
            os.mkdir(str(i) + "_" + item)

        input(
            "Please put a .ins and .hkl into each folder. After all files have been moved, press any key to continue..."
        )

        # DON'T use enumerate here, because stats_location and results_location would also contribute to numbers

        index = 0

        for item in self.sorter.sorted_properly(os.listdir()):
            if (
                item != self.stats_location.name
                and item != self.results_location.name
                and os.path.isdir(item) == True
            ):

                try:
                    shutil.copy(instrument_path, item)
                except shutil.SameFileError:
                    print(
                        "Error! Please put your reference outside of the working directory"
                    )
                    exit()

                os.chdir(item)

                cif = ReadCif(instrument_path.name)
                data_block = cif.first_block()
                data_block[varying_param] = varying_data[index]

                with open("instrument.cif", "w") as f:
                    f.write(cif.WriteOut())

                index += 1
                os.chdir("..")

        self.chemical_formula = ""
        self.crystal_habit = ""
        self.crystal_colour = ""
        self.max_dimension = ""
        self.mid_dimension = ""
        self.min_dimension = ""

    def process(
        self,
        location: str,
        reference: str,
        graph_output_location: str,
        refinements_to_check: int,
        tolerance: float,
        max_cycles: int,
    ) -> None:

        """Runs SHELXL on a series of structures based on one reference

        Args:
            location (str): full path to the folder containing folders of .ins files
            reference (str): full path to the reference .ins/.res file
            graph_output_location (str): full path to the location of output files
            refinements_to_check (int): number of refinements to check for shift convergence
            tolerance (float): target shift value
            max_cycles (int): maximum cycles SHELXL can run before stopping
        """

        shelxl = Refinement_Pipeline(self.test_mode)
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
        varying_param: str,
        instrument_ending: str,
        instrument_file: str,
        additional_params: list,
        adps: bool,
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
            varying_param (str): which parameter in the CIF is varying over the experiments
                                IE _diffrn_ambient_temperature for a VT experiment
            instrument_ending(str): ending of the instrument file if they are
                                    all consistent with different names
            instrument_file(str): name of the instrument file if they are all consistent
            additional_user_parameters(list): list of extra cif parameters the user has edited
                                            CURRENTLY UNDER DEVELOPMENT AS OPTION NOT IMPLEMENTED
            adps (bool): whether or not ADP analysis should be run
        """

        cif = CIF_Compile_Pipeline(self.test_mode)
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
            additional_params,
        )
        cif.compile_cifs(results_location, [self.stats_location, self.results_location])
        analysis = Variable_Analysis_Pipeline(self.test_mode)
        analysis.analyse_data(
            reference,
            results_location,
            cif_parameters,
            atoms_for_analysis,
            varying_param,
            bonds,
            angles,
            torsions,
            adps,
        )
