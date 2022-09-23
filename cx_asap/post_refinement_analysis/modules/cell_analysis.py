#!/usr/bin/env python3

###################################################################################################
# --------------------------------------CX-ASAP: cell_analysis-------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

import logging
import os
import pathlib

import pandas as pd
from system_files.utils import Cell_Import, Config, Grapher, Nice_YAML_Dumper

# ----------Class Definition----------#


class Cell_Deformation:
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

    def import_data(self, csv: str, ref: str) -> None:

        """Imports unit cell data from a .csv file

        BEST TO RUN CIF_READ.py FIRST

        Also imports a reference .ins file to extract a reference unit cell

        Args:
            csv (str): full path to .csv file with cell data
            ref (str): full path to reference .ins file
        """

        # Imports data from a .csv file

        self.csv = pathlib.Path(csv)
        try:
            self.df = pd.read_csv(self.csv)
        except:
            logging.critical("No data in csv. If used in pipeline, this likely means no datasets successfully refined.")
            print("Check Error Logs")
            exit()

        ref_cell = Cell_Import(self.test_mode)
        ref_cell.cell_import(ref)

        # Reloading sys.yaml file to make sure that it is accepting the updated reference cell

        config = Config(self.test_mode)

        self.sys = config.sys
        self.sys_path = config.sys_path

    def import_dataframe(self, df: "pd.Dataframe") -> None:

        """Option to import data directly from a dataframe

        This is more useful in larger pipelines so that the user

        does not have to define the data path and the step of

        exporting and reimporting a .csv can be skipped

        Args:
            df ("pd.Dataframe") = dataframe with unit cell data
        """

        # This is the option to directly import a dataframe rather than a .csv file

        self.df = df

    def deformation(self, df: "pd.DataFrame", ref: list) -> "pd.DataFrame":

        """Performs the actual calculations for the deformation of all unit cell parameters in the

        imported dataframes based on the reference unit cell

        Args:
            df (pd.DataFrame) = dataframe with unit cell information
            ref(list) = list of reference unit cell parameters in the order:
                        a, b, c, alpha, beta, gamma, volume

        """

        # Calculating the cell deformations based on the references in the sys.yaml file

        df["a_axis_deformation"] = ((df["_cell_length_a"] / ref[0]) - 1) * 100
        df["b_axis_deformation"] = ((df["_cell_length_b"] / ref[1]) - 1) * 100
        df["c_axis_deformation"] = ((df["_cell_length_c"] / ref[2]) - 1) * 100
        df["alpha_deformation"] = ((df["_cell_angle_alpha"] / ref[3]) - 1) * 100
        df["beta_deformation"] = ((df["_cell_angle_beta"] / ref[4]) - 1) * 100
        df["gamma_deformation"] = ((df["_cell_angle_gamma"] / ref[5]) - 1) * 100
        df["volume_deformation"] = ((df["_cell_volume"] / ref[6]) - 1) * 100

        return df

    def calculate_deformations(self) -> None:

        """Calculates the deformation of all unit cell parameters in the

        imported dataframes based on the reference unit cell

        """

        # testing to make sure the headings are of the correct format

        try:
            test = self.df["_cell_length_a"]
        except:
            logging.critical(
                __name__
                + " : CIF headings for unit cell parameters not found in imported data - please use CIF syntax, ie _cell_length_a"
            )
            print("Check Error Logs")
            exit()

        system_df = [
            self.sys["ref_a"],
            self.sys["ref_b"],
            self.sys["ref_c"],
            self.sys["ref_alpha"],
            self.sys["ref_beta"],
            self.sys["ref_gamma"],
            self.sys["ref_volume"],
        ]

        self.df = self.deformation(self.df, system_df)

        # Writes the values out to the csv file for future reference

        self.df.to_csv(self.csv, index=None)

    def graphical_analysis(
        self,
        x_axis: str,
        x_label: str,
        title: str = "Cell Deformation",
        figure_title: str = "cell_parameters.png",
        figure_title_1: str = "axis_deformation.png",
        figure_title_2: str = "angle_deformation.png",
        df: "pd.Dataframe" = None,
    ) -> None:

        """Makes graphs of the changes in cell axes deformation,

        cell angle deformation, and unit cell parameter deformation

        Args:
            x_axis (str): heading in dataframe for x-data
            x_label (str): label for x-axis
            title (str): graph title
            figure_title (str): output file name for cell parameter graph
            figure_title_1 (str): output file name for axis deformation graph
            figure_title_2 (str): output file name for angle deformation graph
            df (pd.Dataframe): dataframe for plotting
        """

        # Will put the graphs in the same directory as the .csv file for analysis

        os.chdir(self.csv.parent)

        analysis = Grapher(self.test_mode)

        try:
            test = df[x_axis]
        except TypeError:
            df = self.df

        x = df[x_axis]

        y1 = [
            df["a_axis_deformation"],
            df["b_axis_deformation"],
            df["c_axis_deformation"],
            df["volume_deformation"],
        ]
        y2 = [
            df["alpha_deformation"],
            df["beta_deformation"],
            df["gamma_deformation"],
            df["volume_deformation"],
        ]
        y3 = [
            df["_cell_length_a"],
            df["_cell_length_b"],
            df["_cell_length_c"],
            df["_cell_angle_alpha"],
            df["_cell_angle_beta"],
            df["_cell_angle_gamma"],
            df["_cell_volume"],
        ]

        # A graph of the axis deformations will be printed

        analysis.single_scatter_graph(
            x,
            y1,
            x_label,
            "Deformation (%)",
            title,
            figure_title_1,
            [
                "$\epsilon$ a-axis",
                "$\epsilon$ b-axis",
                "$\epsilon$ c-axis",
                "$\epsilon$ volume",
            ],
            ["red", "green", "blue", "black"],
            ["x", "+", "$*$", "."],
            [1, 1, 0.5, 0.75],
            [100, 100, 200, 100],
        )

        # A graph of the angle deformations will also be printed separately

        analysis.single_scatter_graph(
            x,
            y2,
            x_label,
            "Deformation (%)",
            title,
            figure_title_2,
            [
                "$\epsilon$ alpha",
                "$\epsilon$ beta",
                "$\epsilon$ gamma",
                "$\epsilon$ volume",
            ],
            ["red", "green", "blue", "black"],
            ["x", "+", "$*$", "."],
            [1, 1, 0.5, 0.75],
            [100, 100, 200, 100],
        )

        # A graph of just the unit cell parameters will be output

        analysis.multi_scatter_graph(
            x,
            y3,
            ["Cell Axes", "Cell Angles", "Volume"],
            3,
            3,
            [1, 2, 3, 4, 5, 6, 8],
            ["a-axis", "b-axis", "c-axis", "alpha", "beta", "gamma", "volume"],
            x_label,
            [
                r"Distance ($\AA$)",
                r"Distance ($\AA$)",
                r"Distance ($\AA$)",
                "Angle(" + chr(176) + ")",
                "Angle(" + chr(176) + ")",
                "Angle(" + chr(176) + ")",
                "Volume (\u212B\u00B3)",
            ],
            figure_title,
        )

    def quality_analysis(
        self,
        x_axis: str,
        y_dict: dict,
        x_title: str,
        title: str = "Quality Statistics",
        figure_name: str = "Quality_Statistics.png",
        df: "pd.Dataframe" = None,
    ) -> None:

        """Makes graphs of the data quality statistics from dataframe

        Args:
            x_axis (str): heading in dataframe for x-data
            y_dict (dict): headers + data for a range of quality statistics
                            (usually R1, Rint and Completeness)
            x_title (str): label for x-axis
            title (str): graph title
            figure_name (str): output file name for quality statistics graph
            df (pd.Dataframe): dataframe for plotting
        """

        # Finally, a graph of the statistics will be output

        os.chdir(self.csv.parent)

        analysis = Grapher(self.test_mode)

        try:
            test = df[x_axis]
        except TypeError:
            df = self.df

        x = df[x_axis]

        analysis.multi_multi_scatter_graph(x, y_dict, title, 3, 1, [1, 2, 3], x_title, "Statistic", figure_name)
