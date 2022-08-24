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
import pandas as pd
import logging

# ----------Class Definition----------#


class VP_Analysis_Pipeline:
    def __init__(self) -> None:

        """Initialises the class

        Sets up the yaml parameters input by the user

        Also defines the location for the system yaml file

        which stores a yaml of code-only parameters accessible throughout

        the software package
        """

        # Setup yaml files and logger

        self.config = Config()

        self.cfg = self.config.cfg
        self.sys = self.config.sys
        self.conf_path = self.config.conf_path
        self.sys_path = self.config.sys_path

        self.extra_parameters = {}

        self.parameters = ["angles", "spot_MC", "min_pixels", "strong_pixels", "sepmin"]

        for item in self.parameters:
            self.extra_parameters[item] = []

    def analyse_data(
        self,
        ref_cell: str,
        location: str,
        cif_parameters: list,
        atoms_for_analysis: list,
        atoms_for_plane: str,
        step_size: int,
        spot_MC: list,
        min_pixels: list,
        strong_pixels: list,
        sepmin: list,
        wedge_angles: list,
        reference_plane: list,
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

        As this is for a flexible mapping experiment, additional information about

        the XDS parameters are included so that the output data is linked to which

        indexation parameters were used for it

        Also analyses the changes in angle between the MPLA command and a reference plane

        Args:
            ref_cell (str): full path to reference .ins file
            location (str): full path to the folder containing all CIFs for analysis
            cif_parameters (list): list of cif parameters to extract
            atoms_for_analysis (list): list of important atoms to separate
            atoms_for_plane (str): atoms used for mean plane analysis
            step_size (int): step size used in mapping experiment
            spot_MC (list): list of SPOT_MAXIMUM-CENTROID values to test (XDS param)
            min_pixels (list): list of MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT values to test (XDS param)
            strong_pixels (list): list of STRONG_PIXEL values to test (XDS param)
            sepmin (list): list of SEPMIN values to test (XDS param)
            wedge_angles (list): list of wedge angles to test
            reference_plane (list): crystallographic plane to complare MPLA with
                                    Eg [1,0,0] to compare MPLA to the (100) plane
            bonds (bool): whether or not bond analysis should be run
            angles (bool): whether or not angle analysis should be run
            torsions (bool): whether or not torsion analysis should be run
            adps (bool): whether or not ADP analysis should be run
        """

        CIF_Data = CIF_Read()
        CIF_Data.configure(cif_parameters)
        CIF_Data.get_data(location, bonds, angles, torsions, adps)
        CIF_Data.data_output()

        cell = Cell_Deformation()
        cell.import_data("CIF_Parameters.csv", ref_cell)

        self.cfg, self.sys = self.config.yaml_reload()

        cell.calculate_deformations()

        # Calculating the movement from the step-size and the position number

        distance_moved = [
            int(x) * int(step_size) for x in self.sys["Successful_Positions"]
        ]
        atoms = [atoms_for_plane for x in self.sys["Successful_Positions"]]
        ref_plane = [reference_plane for x in self.sys["Successful_Positions"]]

        # Duplicates the different parameters tested by the number of structures they apply to so list lengths match when adding to the data frame

        index = 0
        for item1 in spot_MC:
            for item2 in min_pixels:
                for item3 in strong_pixels:
                    for item4 in sepmin:
                        for item5 in wedge_angles:
                            count = 0
                            while count < int(
                                self.sys["Structures_in_each_CIF"][index]
                            ):
                                count += 1
                                self.extra_parameters["angles"].append(item5)
                                self.extra_parameters["spot_MC"].append(item1)
                                self.extra_parameters["min_pixels"].append(item2)
                                self.extra_parameters["strong_pixels"].append(item3)
                                self.extra_parameters["sepmin"].append(item4)
                            index += 1

        # Adding all of the extra data relevant to flexible mapping experiments

        cell.df["Wedge Angle"] = self.extra_parameters["angles"]
        cell.df["Spot Maximum-Centroid"] = self.extra_parameters["spot_MC"]
        cell.df["Minimum Number of Pixels in a Spot"] = self.extra_parameters[
            "min_pixels"
        ]
        cell.df["Strong Pixel"] = self.extra_parameters["strong_pixels"]
        cell.df["Sepmin"] = self.extra_parameters["sepmin"]
        cell.df["Distance"] = distance_moved
        cell.df["Position"] = self.sys["Successful_Positions"]
        # cell.df['Rotation'] = self.angles_df['Rotation Angle']
        cell.df["Components in Mean Plane"] = atoms
        cell.df["Reference Plane"] = ref_plane

        graph = Grapher()

        # Separate DataFrames by CIF Name

        discrete_cif_names = list(dict.fromkeys(cell.df["CIF_File"]))

        separated_by_cif_dfs = []

        for item in discrete_cif_names:
            condition = cell.df["CIF_File"] == item
            separated_by_cif_dfs.append(cell.df[condition])

        for index, item in enumerate(separated_by_cif_dfs):
            item.to_csv(discrete_cif_names[index] + "_parameters.csv", index=None)

            # Statistics Graph

            cell.quality_analysis(
                "Distance",
                {
                    "R1": item["_refine_ls_R_factor_gt"],
                    "Rint": item["_diffrn_reflns_av_R_equivalents"],
                    "Completeness": item["_diffrn_measured_fraction_theta_full"],
                },
                "Distance($\mu$m)",
                df=item,
            )

            # Deformation Graph

            cell.graphical_analysis(
                "Distance",
                "Distance($\mu$m)",
                "Cell Deformation - " + discrete_cif_names[index],
                "cell_parameters_" + discrete_cif_names[index] + ".png",
                "axis_deformation_" + discrete_cif_names[index] + ".png",
                "angle_deformation_" + discrete_cif_names[index] + ".png",
                item,
            )

        # ADPs

        if adps != False:
            adp_df = pd.read_csv("ADPs.csv")
            discrete_cif_names_bond = list(dict.fromkeys(bond_df["CIF_File"]))
            separated_by_cif_bond = []

            for item in discrete_cif_names_bond:
                condition = bond_df["CIF_File"] == item
                separated_by_cif_bond.append(bond_df[condition])

            for index, item in enumerate(separated_by_cif_bond):
                item.to_csv("ADPs_" + discrete_cif_names_bond[index] + ".csv")

                adp_object = ADP_analysis(self.test_mode)
                adp_object.analyse_data(
                    "ADPs_" + discrete_cif_names_bond[index] + ".csv",
                    discrete_cif_names[index] + "_parameters.csv",
                )

        # Bonds/Angles/Torsions

        geometry = Structural_Analysis()

        if bonds != False:
            bond_df = pd.read_csv("Bond_Lengths.csv")
            discrete_cif_names_bond = list(dict.fromkeys(bond_df["CIF_File"]))
            separated_by_cif_bond = []

            for item in discrete_cif_names_bond:
                condition = bond_df["CIF_File"] == item
                separated_by_cif_bond.append(bond_df[condition])

            for index, item in enumerate(separated_by_cif_bond):
                item.to_csv("Bond_Lengths_" + discrete_cif_names_bond[index] + ".csv")

                geometry.import_and_analyse(
                    "Bond_Lengths_" + discrete_cif_names_bond[index] + ".csv",
                    False,
                    False,
                    atoms_for_analysis,
                    location,
                    str(index + 1),
                    True,
                )

        if angles != False:
            angle_df = pd.read_csv("Bond_Angles.csv")
            discrete_cif_names_angle = list(dict.fromkeys(angle_df["CIF_File"]))
            separated_by_cif_angle = []

            for item in discrete_cif_names_angle:
                condition = angle_df["CIF_File"] == item
                separated_by_cif_angle.append(angle_df[condition])

            for index, item in enumerate(separated_by_cif_angle):
                item.to_csv("Bond_Angles_" + discrete_cif_names_bond[index] + ".csv")

                geometry.import_and_analyse(
                    False,
                    "Bond_Angles_" + discrete_cif_names_bond[index] + ".csv",
                    False,
                    atoms_for_analysis,
                    location,
                    str(index + 1),
                    True,
                )

        if torsions != False:
            torsion_df = pd.read_csv("Bond_Torsions.csv")
            discrete_cif_names_torsion = list(dict.fromkeys(torsion_df["CIF_File"]))
            separated_by_cif_torsion = []

            for item in discrete_cif_names_torsion:
                condition = torsion_df["CIF_File"] == item
                separated_by_cif_torsion.append(torsion_df[condition])

            for index, item in enumerate(separated_by_cif_torsion):
                item.to_csv("Bond_Torsions_" + discrete_cif_names_bond[index] + ".csv")

                geometry.import_and_analyse(
                    False,
                    False,
                    "Bond_Torsions_" + discrete_cif_names_bond[index] + ".csv",
                    atoms_for_analysis,
                    location,
                    str(index + 1),
                    True,
                )
