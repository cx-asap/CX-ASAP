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


class VT_Analysis_Pipeline:
    def __init__(self) -> None:

        """Initialises the class

        Sets up the yaml parameters input by the user

        Also defines the location for the system yaml file

        which stores a yaml of code-only parameters accessible throughout

        the software package
        """

        # Setup yaml files and logger

        config = Config()

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    def determine_temp_behaviour(self, df: "pd.DataFrame") -> list:

        """Searches through the data and classifies everything as a

        minima, maxima, increasing or decreasing

        This is done so that graphical output can have different colours

        for increasing vs decreasing

        Makes it easier to understand any hysteresis

        THIS FUNCTION IS SPECIFIC TO VT, SO LOOKS FOR "_diffrn_ambient_temperature"

        Args:
            df ("pd.DataFrame") = dataframe containing all the... data

        Returns:
            behaviour (list) = a list of the behaviour with the same length
                                as the data in the dataframe
        """

        behaviour = []

        data = df["_diffrn_ambient_temperature"]

        # Searches through the data and classifies everything as a temperature minima, maxima, increasing or decreasing

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
        ref_cell: list,
        ref_ins: str,
        location: str,
        cif_parameters: list,
        atoms_for_analysis: list,
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

        THIS FUNCTION IS SPECIFIC TO VT, SO HAS CHANGING PARAM AS "_diffrn_ambient_temperature"

        Args:
            ref_cell (list): reference unit cell 
            ref_ins (str): full path to the reference .ins
            location (str): full path to the folder containing all CIFs for analysis
            cif_parameters (list): list of cif parameters to extract
            atoms_for_analysis (list): list of important atoms to separate
            bonds (bool): whether or not bond analysis should be run
            angles (bool): whether or not angle analysis should be run
            torsions (bool): whether or not torsion analysis should be run
            adps (bool): whether or not ADP analysis should be run
        """

        CIF_Data = CIF_Read()
        CIF_Data.configure(cif_parameters)
        CIF_Data.get_data(location, bonds, angles, torsions, adps)
        CIF_Data.data_output()

        geometry = Structural_Analysis()

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
            varying_parameter="_diffrn_ambient_temperature",
        )

        cell = Cell_Deformation()
        cell.import_data("CIF_Parameters.csv", ref_ins)
        cell.calculate_deformations()
        cell.quality_analysis(
            "_diffrn_ambient_temperature",
            {
                "R1": cell.df["_refine_ls_R_factor_gt"],
                "Rint": cell.df["_diffrn_reflns_av_R_equivalents"],
                "Completeness": cell.df["_diffrn_measured_fraction_theta_full"],
            },
            "Temperature (K)",
        )
        cell.graphical_analysis("_diffrn_ambient_temperature", "Temperature (K)")

        if adps != False:
            adp_object = ADP_analysis()
            adp_object.analyse_data(adps, "CIF_Parameters.csv")

        graph = Grapher()
        discrete_behaviour = list(dict.fromkeys(self.determine_temp_behaviour(cell.df)))
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
