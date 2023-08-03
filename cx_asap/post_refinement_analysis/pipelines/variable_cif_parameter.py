#!/usr/bin/env python3

###################################################################################################
# -------------------------------CX-ASAP: Variable_Temperature_Analysis----------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config, Grapher
from post_refinement_analysis.modules.cif_read import CIF_Read
from post_refinement_analysis.modules.cell_analysis import Cell_Deformation
from post_refinement_analysis.modules.structural_analysis import Structural_Analysis
from post_refinement_analysis.modules.ADP_analysis import ADP_analysis
import yaml
import logging

# ----------Class Definition----------#


class Variable_Analysis_Pipeline:
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

    def determine_behaviour(self, df: "pd.DataFrame", param: str) -> list:

        """Searches through the data and classifies everything as a

        minima, maxima, increasing or decreasing

        This is done so that graphical output can have different colours

        for increasing vs decreasing

        Makes it easier to understand any hysteresis

        THIS FUNCTION CAN BE CUSTOMISED FOR ANY CHANGING CIF PARAMETER

        Args:
            df ("pd.DataFrame") = dataframe containing all the... data
            param (str) = heading in the dataframe for the variable of interest
                            function analyses how THIS param in the df is changing

        Returns:
            behaviour (list) = a list of the behaviour with the same length
                                as the data in the dataframe
        """

        behaviour = []

        data = df[param]

        for index, i in enumerate(data):
            if index != 0 and index != len(data) - 1:
                if i < data[index - 1] and i < data[index + 1]:
                    behaviour.append("Minima")
                elif i > data[index - 1] and i > data[index + 1]:
                    behaviour.append("Maxima")
                elif i < data[index - 1] and i > data[index + 1]:
                    behaviour.append("Decreasing")
                elif i > data[index - 1] and i < data[index + 1]:
                    behaviour.append("Increasing")
                elif i > data[index - 1] and i == data[index + 1]:
                    behaviour.append("Increasing")
                elif i < data[index - 1] and i == data[index + 1]:
                    behaviour.append("Decreasing")
                elif i == data[index - 1]:
                    behaviour.append("Did Not Change")
                else:
                    behaviour.append("Error")

            # This else statement is needed to classify the first and last points

            else:
                if index == 0:
                    if i < data[index + 1]:
                        behaviour.append("Increasing")
                    elif i == data[index + 1]:
                        behaviour.append("Did Not Change")
                    else:
                        behaviour.append("Decreasing")
                else:
                    if i < data[index - 1]:
                        behaviour.append("Decreasing")
                    elif i == data[index - 1]:
                        behaviour.append("Did Not Change")
                    else:
                        behaviour.append("Increasing")

        df["behaviour"] = behaviour

        return behaviour

    def analyse_data(
        self,
        ref_cell: str,
        location: str,
        cif_parameters: list,
        atoms_for_analysis: list,
        param: str,
        bonds: bool = False,
        angles: bool = False,
        torsions: bool = False,
        adps: bool = False,
    ) -> None:

        """Performs much analysis on CIF files

        Reads them and extracts the desired parameters

        Performs structural anaysis on bonds, angles and torsions as requested

        Performs ADP analysis if requested

        Outputs graphs corresponding to these analyses

        THIS FUNCTION CAN BE CUSTOMISED FOR ANY CHANGING CIF PARAMETER

        Args:
            ref_cell (str): full path to reference .ins file
            location (str): full path to the folder containing all CIFs for analysis
            cif_parameters (list): list of cif parameters to extract
            atoms_for_analysis (list): list of important atoms to separate
            param (str): the CIF parameter that is varying in correct CIF syntax
                        ie "_diffrn_ambient_temperature" for a VT experiment
            bonds (bool): whether or not bond analysis should be run
            angles (bool): whether or not angle analysis should be run
            torsions (bool): whether or not torsion analysis should be run
            adps (bool): whether or not ADP analysis should be run
        """

        CIF_Data = CIF_Read(self.test_mode)
        CIF_Data.configure(cif_parameters)
        CIF_Data.get_data(location, bonds, angles, torsions, adps, param)
        CIF_Data.data_output()

        geometry = Structural_Analysis(self.test_mode)

        if bonds != False:
            bonds = "Bond_Lengths.csv"
        if angles != False:
            angles = "Bond_Angles.csv"
        if torsions != False:
            torsions = "Bond_Torsions.csv"
        if adps != False:
            adps = "ADPs.csv"

        geometry.import_and_analyse(
            bonds,
            angles,
            torsions,
            atoms_for_analysis,
            location,
            varying_parameter=param,
        )

        cell = Cell_Deformation(self.test_mode)
        cell.import_data("CIF_Parameters.csv", ref_cell)
        cell.calculate_deformations()
        cell.quality_analysis(
            param,
            {
                "R1": cell.df["_refine_ls_R_factor_gt"],
                "Rint": cell.df["_diffrn_reflns_av_R_equivalents"],
                "Completeness": cell.df["_diffrn_measured_fraction_theta_full"],
            },
            param,
        )
        cell.graphical_analysis(param, param)

        if adps != False:
            adp_object = ADP_analysis(self.test_mode)
            adp_object.analyse_data(adps, "CIF_Parameters.csv")
            

        graph = Grapher(self.test_mode)
        discrete_behaviour = list(
            dict.fromkeys(self.determine_behaviour(cell.df, param))
        )
        separated_by_behaviour_dfs = []
        for item in discrete_behaviour:
            condition = cell.df["behaviour"] == item
            separated_by_behaviour_dfs.append(cell.df[condition])

        y_headers = [
            "_cell_length_a",
            "_cell_length_b",
            "_cell_length_c",
            "_cell_angle_alpha",
            "_cell_angle_beta",
            "_cell_angle_gamma",
            "_cell_volume",
        ]
