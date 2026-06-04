#!/usr/bin/env python3

###################################################################################################
# ------------------------------------CX-ASAP: centroids_pipeline----------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Config, Directory_Browse, Grapher
from post_refinement_analysis.modules.centroids import Centroids
import os
import pandas as pd
import logging

# ----------Class Definition----------#


class Centroids_Pipeline:
    def __init__(self) -> None:
        """Initialises the class

        Sets up the yaml parameters input by the user

        Also defines the location for the system yaml file

        which stores a yaml of code-only parameters accessible throughout

        the software package
        """

        config = Config()

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    def centroid_distance_analysis(
        self,
        working_directory: str,
        atom_list_1: list,
        atom_list_2: list,
        results_directory: str,
        label: str = "Centroid Distance",
    ) -> None:
        """Calculates centroid-to-centroid distances for a series of .lst files
        in separate folders within a common parent folder.

        Produces a scatter graph of centroid distance vs structure number.

        Args:
            working_directory (str): full path to the parent folder containing
                                     folders with .lst files
            atom_list_1 (list): atom labels for the first centroid group
            atom_list_2 (list): atom labels for the second centroid group
            results_directory (str): full path to the output directory
            label (str): column header for the output csv and graph y-axis
        """

        centroid = Centroids()
        tree = Directory_Browse(working_directory)

        for index, item in enumerate(tree.directories):
            tree.enter_directory(item, ".lst")
            centroid.analysis_centroid_distance(
                tree.item_file,
                index + 1,
                results_directory,
                atom_list_1,
                atom_list_2,
                label,
            )
            tree.exit_directory()

        os.chdir(results_directory)
        csv_name = label.lower().replace(" ", "_") + ".csv"

        try:
            full_data = pd.read_csv(csv_name)
            x = full_data["Structure"]
            y = full_data[label]
            graph = Grapher()
            graph.single_scatter_graph(
                x,
                y,
                "Structure Number",
                r"Distance ($\AA$)",
                label,
                csv_name.replace(".csv", ".png"),
            )
        except FileNotFoundError:
            logging.error(f"No {csv_name} file found...")


