#!/usr/bin/env python3

###################################################################################################
# -----------------------------------CX-ASAP: rotation_pipeline------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config, Directory_Browse, Grapher
from post_refinement_analysis.modules.rotation_planes import Rotation
import os
import pandas as pd
import logging

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

    def analysis(
        self, working_directory: str, reference_plane: list, results_directory: str
    ) -> None:
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

        try:
            full_data = pd.read_csv("rotation_angles.csv")
            x = full_data["Structure"]
            angle_cols = [
                c
                for c in full_data.columns
                if c.startswith("MPLA_") and c.endswith("_Rotation_Angle")
            ]
            y_data = [list(full_data[c]) for c in angle_cols]
            labels = [c.replace("_", " ") for c in angle_cols]
            graph = Grapher()
            graph.single_scatter_graph(
                x,
                y_data,
                "Structure Number",
                "Angle($^\circ$)",
                "Rotation Angles",
                "rotation_angles.png",
                y_series_title=labels if len(labels) > 1 else None,
            )
        except FileNotFoundError:
            logging.error("No rotation angles file found...")
            pass

    def interplane_angle_analysis(
        self, working_directory: str, results_directory: str
    ) -> None:
        """Extracts the SHELXL inter-plane angle for a series of .lst files
        in separate folders within a common parent folder.

        Requires two MPLA commands in the .ins file for each refinement.

        Produces a scatter graph of inter-plane angle vs structure number.

        Args:
            working_directory (str): full path to the parent folder containing
                                     folders with .lst files
            results_directory (str): full path to the output directory
        """

        plane = Rotation()
        tree = Directory_Browse(working_directory)

        for index, item in enumerate(tree.directories):
            tree.enter_directory(item, ".lst")
            plane.analysis_interplane_angle(
                tree.item_file,
                index + 1,
                results_directory,
            )
            tree.exit_directory()

        os.chdir(results_directory)

        try:
            full_data = pd.read_csv("interplane_angles.csv")
            x = full_data["Structure"]
            y = full_data["Interplane Angle"]
            graph = Grapher()
            graph.single_scatter_graph(
                x,
                y,
                "Structure Number",
                "Angle($^\circ$)",
                "Inter-plane Angles",
                "interplane_angles.png",
            )
        except FileNotFoundError:
            logging.error("No interplane_angles.csv file found...")
