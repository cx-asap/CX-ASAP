#!/usr/bin/env python3

###################################################################################################
# ------------------------------CX-ASAP: variable_temperature_pipeline-----------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config, Autoprocess_Setup, Cell_Import
from data_reduction.pipelines.xprep_pipeline import XPREP_Pipeline
from data_refinement.pipelines.refine_pipeline import Refinement_Pipeline
from cif_validation.pipelines.cif_pipeline import CIF_Compile_Pipeline
from post_refinement_analysis.pipelines.variable_temperature_analysis import (
    VT_Analysis_Pipeline,
)
from post_refinement_analysis.pipelines.rotation_pipeline import Rotation_Pipeline
import yaml
import logging

# ----------Class Definition----------#


class VT_Pipeline:
    def __init__(self):

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

        self.config = Config()

        self.cfg = self.config.cfg
        self.sys = self.config.sys
        self.conf_path = self.config.conf_path
        self.sys_path = self.config.sys_path

    def setup(
        self,
        location: str,
        experiment_name: str,
        experiment_type: str,
        reference_location: str,
        MPLA_atoms: str,
    ) -> None:

        """Organises the directory tree for this experiment

        Saves new locations into the system yaml file

        Imports the reference .ins/.res file

        Also inserts the MPLA command into the reference according to listed atoms

        Args:
            location (str): full path to the location of all the autoprocessed datasets
            experiment_name (str): common name between autoprocessed datasets
            experiment_type (str): type of experiment (only used for folder name)
            reference_location (str): full path to the reference .ins/.res
            MPLA_atoms (str): Atoms for mean plane analysis

        """

        setup = Autoprocess_Setup()
        setup.config(location, experiment_name, experiment_type)
        setup.Organise_Directory_Tree(reference_location, True)
        reference_setup = Cell_Import()
        reference_setup.cell_import(reference_location)
        reference_setup.ref_edit(reference_location, MPLA_atoms)
        self.cfg, self.sys = self.config.yaml_reload()

    def process(
        self,
        matrix: str,
        location: str,
        file_name: str,
        space_group: str,
        chemical_formula: str,
        reference: str,
        graph_output_location: str,
        refinements_to_check: int,
        tolerance: float,
        max_cycles: int,
        reference_plane: list,
    ) -> None:

        """This is the function that actually calls XPREP and SHELXL

        Args:
            matrix (str): matrix transformation for input into xprep
            location (str): full path to the folder containing folders of .ins files
            file_name (str): filename for xprep (ie 'XDS_ASCII.HKL_p1')
            space_group (str): space group as a string for input into xprep
            chemical_formula (str): chemical formula of the crystal
            reference (str): full path to the reference .ins/.res file
            graph_output_location (str): full path to the location of output files
            refinements_to_check (int): number of refinements to check for shift convergence
            tolerance (float): target shift value
            max_cycles (int): maximum cycles SHELXL can run before stopping
            reference_plane (list): crystallographic plane to complare MPLA with
                                    Eg [1,0,0] to compare MPLA to the (100) plane
        """

        xprep = XPREP_Pipeline()
        xprep.multiple_xprep(matrix, location, file_name, space_group, chemical_formula)
        shelxl = Refinement_Pipeline()
        shelxl.multiple_refinement(
            location,
            reference,
            graph_output_location,
            refinements_to_check,
            tolerance,
            max_cycles,
        )
        rotation = Rotation_Pipeline()
        rotation.analysis(location, reference_plane, graph_output_location)

    def analyse(
        self,
        ref_cell: list,
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
    ):

        """Compiles all of the CIFs and performs analysis on them

        Args:
            ref_cell (list): reference unit cell
            reference (str): full path to the reference .ins/.res file
            data_location (str): full path to the folder containing folders of CIFs
            results_location (str): full path to the location of output files
            solution_program(str): the structure solution program used
            cheical_formula (str): the chemical formula of the crystal
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
        cif.compile_cifs(results_location)
        analysis = VT_Analysis_Pipeline()
        analysis.analyse_data(
            ref_cell,
            reference,
            results_location,
            cif_parameters,
            atoms_for_analysis,
            bonds,
            angles,
            torsions,
            adps,
        )
