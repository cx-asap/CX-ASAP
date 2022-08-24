#!/usr/bin/env python3

###################################################################################################
# ------------------------------CX-ASAP: variable_position_pipeline-----------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import (
    Nice_YAML_Dumper,
    Config,
    Reprocess_Setup,
    Configure_Flexible,
    Cell_Import,
    XDS_File_Edit,
)
from data_reduction.pipelines.xds_pipeline import XDS_Pipeline
from data_reduction.pipelines.xprep_pipeline import XPREP_Pipeline
from data_refinement.pipelines.refine_pipeline import Refinement_Pipeline
from cif_validation.pipelines.cif_pipeline import CIF_Compile_Pipeline
from post_refinement_analysis.pipelines.rotation_pipeline import Rotation_Pipeline
from post_refinement_analysis.pipelines.variable_position_analysis import (
    VP_Analysis_Pipeline,
)
import yaml
import pathlib
import os
import logging

# ----------Class Definition----------#


class VP_Pipeline:
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

        self.xds_edit = XDS_File_Edit()
        self.sys["process_counter"] = 0
        with open(self.sys_path, "w") as f:

            yaml.dump(
                self.sys,
                f,
                default_flow_style=False,
                Dumper=Nice_YAML_Dumper,
                sort_keys=False,
            )

        self.cfg, self.sys = self.config.yaml_reload()

    def setup(
        self,
        location: str,
        experiment_name: str,
        experiment_type: str,
        reference_location: str,
        xds_inp_location: str,
        background_reference: str,
        instrument_parameters_path: str,
        max_processors: int,
        neggia_library: str,
        space_group_number: int,
        MPLA_atoms: str,
        instrument_cif_path: str,
        total_angle: int,
    ) -> None:

        """Organises the directory tree for this experiment

        Saves new locations into the system yaml file

        Also inserts the MPLA command into the reference according to listed atoms

        Args:
            location (str): full path to the location of all the raw frames
            experiment_name (str): common name between frames of different datasets
            experiment_type (str): type of experiment (only used for folder name)
            reference_location (str): full path to the reference .ins/.res
            xds_inp_location (str): full path to the reference XDS.INP
            background_reference (str): full path to the reference XDS background files
                                        BKGINIT.cbf
                                        BLANK.cbf
                                        GAIN.cbf
                                        X-CORRECTIONS.cbf
                                        Y-CORRECTIONS.cbf
            instrument_parameters_reference (str): full path to the reference GXPARM.XDS
            max_processors (int): maximum number of processors for XDS to use
            neggia_library (str): full path to the dectris-neggia.so file
            space_group_number (int): expected space group number
            MPLA_atoms (str): Atoms for mean plane analysis
            instrument_cif_path (str): full path to the instrument CIF
            total_angle (int): total wedge angle per experiment

        """

        setup = Reprocess_Setup()
        setup.config(location, experiment_name, experiment_type)
        setup.Organise_Directory_Tree(reference_location, xds_inp_location)
        configure = Configure_Flexible()
        configure.file_movement(
            background_reference, instrument_parameters_path, instrument_cif_path
        )
        configure.flexible_setup(
            max_processors, neggia_library, space_group_number, MPLA_atoms, total_angle
        )
        self.cfg, self.sys = self.config.yaml_reload()

    def change_angle(self, total_angle: int, start_angle: int, angle: int) -> None:

        """Calculates frames per degree of data

        Also calculates what the middle angle is from the total angle

        Does this so that can only analyse (for example) 20 degrees

        from a 40 degree wedge, and it centres this around the middle

        This is required for flexible experiments

        Args:
            total_angle (int): total wedge angle of data per dataset
            start_angle (int): starting angle of data collection
                                (would be read in from XDS.INP)
            angle (int): desired wedge angle for analysis
        """

        angles_to_edit = [
            "STARTING_ANGLE",
            "STARTING_FRAME",
            "DATA_RANGE",
            "SPOT_RANGE",
        ]

        middle_angle = total_angle / 2 + start_angle

        middle_frame = round((int(self.sys["total_frames"]) + 1) / 2)

        if angle == total_angle:
            starting_frame = "1"
            ending_frame = self.sys["total_frames"]
            starting_angle = str(start_angle)
        else:
            half_angle = angle / 2
            frames_per_degree = float(self.sys["total_frames"]) / float(total_angle)
            frames_either_side_of_middle = round(half_angle * frames_per_degree)

            starting_frame = str(middle_frame - frames_either_side_of_middle)
            ending_frame = str(middle_frame + frames_either_side_of_middle)

            starting_angle = str((start_angle + (total_angle - angle)) - half_angle)

        self.xds_edit.change(
            self.sys["XDS_inp_organised"], "STARTING_ANGLE", starting_angle
        )
        self.xds_edit.change(
            self.sys["XDS_inp_organised"], "STARTING_FRAME", starting_frame
        )
        self.xds_edit.change(
            self.sys["XDS_inp_organised"],
            "DATA_RANGE",
            starting_frame + " " + ending_frame,
        )
        self.xds_edit.change(
            self.sys["XDS_inp_organised"],
            "SPOT_RANGE",
            starting_frame + " " + ending_frame,
        )

    def change_parameters(self, a: int, b: int, c: int, d: int, e: int) -> None:

        """Changes the parameters in the XDS.INP file to those

        chosen in the conf.YAML file by the user

        Args:
            a (int): desired SPOT_MAXIMUM-CENTROID value
            b (int): desired MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT value
            c (int): desired STRONG_PIXEL value
            d (int): desired SEPMIN value
            e (int): desired wedge angle for data analysis
        """

        # Changing parameter function

        self.xds_edit.change(self.sys["XDS_inp_organised"], "SPOT_MAXIMUM-CENTROID", a)
        self.xds_edit.change(
            self.sys["XDS_inp_organised"], "MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT", b
        )
        self.xds_edit.change(self.sys["XDS_inp_organised"], "STRONG_PIXEL", c)
        self.xds_edit.change(self.sys["XDS_inp_organised"], "SEPMIN", d)
        self.xds_edit.change(
            self.sys["XDS_inp_organised"], "CLUSTER_RADIUS", int(d) / 2
        )
        self.change_angle(self.sys["total_angle"], self.sys["start_angle"], e)

    def flexible_parameter_loops(
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
        xds_INP_path: str,
        angles: list,
        min_pixels: list,
        sepmin: list,
        spot_MC: list,
        strong_pixels: list,
        solution_program: str,
        crystal_habit: str,
        crystal_colour: str,
        max_dimension: str,
        middle_dimension: str,
        min_dimension: str,
        instrument_ending: str,
        instrument_file: str,
    ) -> None:

        """Loops through a range of XDS parameters and runs a full process

        and refinement for each possible combination

        The values for testing are defined by the user

        Currently, 5 parameters are automated for testing in the XDS.INP file:

         - Wedge Angle
         - MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT
         - SPOT_MAXIMUM-CENTROID
         - STRONG_PIXEL
         - SEPMIN (this also changes CLUSTER_RADIUS which is set to half of SEPMIN)

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
            xds_INP_path (str): full path to the reference XDS.INP file
            angles (list): list of wedge angles to test
            min_pixels (list): list of MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT values to test (XDS param)
            sepmin (list): list of SEPMIN values to test (XDS param)
            spot_MC (list): list of SPOT_MAXIMUM-CENTROID values to test (XDS param)
            strong_pixels (list): list of STRONG_PIXEL values to test (XDS param)
            solution_program(str): the structure solution program used
            crystal_habit(str): the crystal habit
            crystal_colour(str): the colour of the crystal
            max_dimension(str): the largest dimension of the crystal
            middle_dimension(str): the middle dimension of the crystal
            min_dimension(str): the smallest dimension of the crystal
            instrument_ending(str): ending of the instrument file if they are
                                    all consistent with different names
            instrument_file(str): name of the instrument file if they are all consistent
        """

        # The current possible parameters to vary

        parameters = [angles, min_pixels, sepmin, spot_MC, strong_pixels]

        # Checks that the user is not trying to test an angle of 0 or an angle larger than the collected data

        for item in angles:

            if item == 0:
                logging.critical(__name__ + " : Cannot test angle of 0!")
                print("Error! Check logs")
                exit()
            elif item > self.sys["total_angle"]:
                logging.critical(
                    __name__
                    + " : Cannot test angle larger than wedge of data collected!"
                )
                print("Error! Check logs")
                exit()

        # Checks that the values entered by the user are actually numbers

        for parameter in parameters:
            for value in parameter:
                if type(value) != int and type(value) != float:
                    logging.critical(
                        __name__
                        + " : Please enter numbers only for angles/parameters to test"
                    )
                    print("Error! Check logs")
                    exit()

        # Cycle through all combinations of parameters to explore options

        for item1 in spot_MC:
            for item2 in min_pixels:
                for item3 in strong_pixels:
                    for item4 in sepmin:
                        for item5 in angles:
                            self.sys["process_counter"] += 1

                            with open(self.sys_path, "w") as f:

                                yaml.dump(
                                    self.sys,
                                    f,
                                    default_flow_style=False,
                                    Dumper=Nice_YAML_Dumper,
                                    sort_keys=False,
                                )

                            self.cfg, self.sys = self.config.yaml_reload()

                            self.change_parameters(item1, item2, item3, item4, item5)

                            # Execute processing pipeline here

                            self.process(
                                matrix,
                                location,
                                file_name,
                                space_group,
                                chemical_formula,
                                reference,
                                graph_output_location,
                                refinements_to_check,
                                tolerance,
                                max_cycles,
                                reference_plane,
                                xds_INP_path,
                                solution_program,
                                crystal_habit,
                                crystal_colour,
                                max_dimension,
                                middle_dimension,
                                min_dimension,
                                instrument_ending,
                                instrument_file,
                            )

                            # Rename Files!

                            os.chdir(graph_output_location)

                            for item in os.listdir():
                                if item.split("_")[0].isdigit() == False:
                                    os.rename(
                                        item,
                                        str(self.sys["process_counter"]) + "_" + item,
                                    )

        # Resets back to the starting values to keep the file good

        self.xds_edit.change(
            self.sys["XDS_inp_organised"], "STARTING_ANGLE", self.sys["start_angle"]
        )
        self.xds_edit.change(self.sys["XDS_inp_organised"], "STARTING_FRAME", 1)
        self.xds_edit.change(
            self.sys["XDS_inp_organised"],
            "DATA_RANGE",
            "1 " + str(self.sys["total_frames"]),
        )
        self.xds_edit.change(
            self.sys["XDS_inp_organised"],
            "SPOT_RANGE",
            "1 " + str(self.sys["total_frames"]),
        )

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
        xds_INP_path: str,
        solution_program: str,
        crystal_habit: str,
        crystal_colour: str,
        max_dimension: str,
        middle_dimension: str,
        min_dimension: str,
        instrument_ending: str,
        instrument_file: str,
    ) -> None:

        """This is the function that actually calls XDS, XPREP and SHELXL

        and also compiles the CIFS from a single set of XDS parameters

        This function is called from the above "flexible_parameter_loops" function

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
            xds_INP_path (str): full path to the reference XDS.INP file
            solution_program(str): the structure solution program used
            crystal_habit(str): the crystal habit
            crystal_colour(str): the colour of the crystal
            max_dimension(str): the largest dimension of the crystal
            middle_dimension(str): the middle dimension of the crystal
            min_dimension(str): the smallest dimension of the crystal
            instrument_ending(str): ending of the instrument file if they are
                                    all consistent with different names
            instrument_file(str): name of the instrument file if they are all consistent
        """

        xds = XDS_Pipeline()
        xds.multiple_xds(location, xds_INP_path)
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
        cif = CIF_Compile_Pipeline()
        cif.configure(
            location,
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
        cif.compile_cifs(graph_output_location)

    def analyse(
        self,
        results_location: str,
        bonds: bool,
        angles: bool,
        torsions: bool,
        adps: bool,
        cif_parameters: list,
        atoms_for_analysis: list,
        ref_cell: str,
        atoms_for_plane: str,
        step_size: int,
        spot_MC: list,
        min_pixels: list,
        strong_pixels: list,
        sepmin: list,
        wedge_angles: list,
        reference_plane: list,
    ) -> None:

        """Analyses the combined CIF files of all the different XDS tests

        Ie, each mapping experiment will have a combined CIF of all data sets

        There will be a combined CIF for each XDS parameter tested

        Outputs a series of graphs looking at changes

        Args:
            results_location (str): full path to the folder that will contain results
            bonds (bool): whether or not bond analysis should be run
            angles (bool): whether or not angle analysis should be run
            torsions (bool): whether or not torsion analysis should be run
            adps (bool): whether or not ADP analysis should be run
            cif_parameters (list): list of parameters in the CIF for tabulation
            atoms_for_analysis (list): list of atoms for geometrical analysis
            ref_cell (str): full path to the reference file for the unit cell
            atoms_for_plane (str): atoms used for mean plane analysis
            step_size (int): step size used in mapping experiment
            spot_MC (list): list of SPOT_MAXIMUM-CENTROID values to test (XDS param)
            min_pixels (list): list of MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT values to test (XDS param)
            strong_pixels (list): list of STRONG_PIXEL values to test (XDS param)
            sepmin (list): list of SEPMIN values to test (XDS param)
            wedge_angles (list): list of wedge angles to test
            reference_plane (list): crystallographic plane to complare MPLA with
                                    Eg [1,0,0] to compare MPLA to the (100) plane
        """

        analysis = VP_Analysis_Pipeline()
        analysis.analyse_data(
            ref_cell,
            results_location,
            cif_parameters,
            atoms_for_analysis,
            atoms_for_plane,
            step_size,
            spot_MC,
            min_pixels,
            strong_pixels,
            sepmin,
            wedge_angles,
            reference_plane,
            bonds,
            angles,
            torsions,
            adps,
        )
