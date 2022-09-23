#!/usr/bin/env python3

###################################################################################################
# -----------------------------------CX-ASAP: rotation_pipeline------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

import logging
import os

import pandas as pd
from post_refinement_analysis.modules.rotation_planes import Rotation
from system_files.utils import Config, Directory_Browse, Grapher, Nice_YAML_Dumper

# ----------Class Definition----------#


class Rotation_Pipeline:
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

    def analysis(self, working_directory: str, reference_plane: list, results_directory: str) -> None:

        """Initialises the class

        Analyses the angle between the MPLA output and a user defined

        reference plane for a series of .lst files in separate folders

        in a common parent folder

        Makes a graph at the end showing how this angle changes across

        the different datasets

        Args:
            working_directory (str): full path to the parent folder containing
                                    folders with .lst files
            reference_plane (list): plane to compare MPLA against as a list
                            ie [1,0,0] corresponds to the (100) plane
            results_directory (str): full path to the location for output graphs
        """
        plane = Rotation()
        plane.configure(reference_plane)
        tree = Directory_Browse(working_directory)
        for index, item in enumerate(tree.directories):
            tree.enter_directory(item, ".lst")
            plane.analysis(tree.item_file, index + 1, results_directory)
            tree.exit_directory()
        os.chdir(results_directory)

        full_data = pd.read_csv("rotation_angles.csv")
        x = full_data["Structure"]
        angle = full_data["Rotation Angle"]

        graph = Grapher()
        graph.single_scatter_graph(
            x,
            angle,
            "Structure Number",
            "Angle($^\circ$C)",
            "Rotation Angles",
            "rotation_angles.png",
        )
