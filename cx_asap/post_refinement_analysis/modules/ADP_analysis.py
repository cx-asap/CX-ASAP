#!/usr/bin/env python3

###################################################################################################
# --------------------------------------CX-ASAP: ADP_analysis--------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

import logging
import math
import os
import pathlib
from typing import Tuple

import numpy as np
import pandas as pd
from system_files.utils import Config, Nice_YAML_Dumper

# ----------Class Definition----------#


class ADP_analysis:
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

        self.test_mode = test_mode

        # Setup yaml files and logger

        config = Config(self.test_mode)

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    def scale_vector(self, vector: "np.array", G: "np.array") -> "np.array":

        """Calculates the scaling factor and scales a given vector

        based on the metrix matrix G

        Args:
            vector (np.array): vector to be scaled
            G (np.array): metrix matrix of the unit cell

        Returns:
            scaled_vector (np.array): vector scaled by calculated factor
        """

        dot_product = np.matmul(np.matmul(vector, G), np.transpose(vector))
        scaling_factor = math.sqrt(dot_product[0][0])

        scaled_vector = vector / scaling_factor

        return scaled_vector

    def calculate_angle(
        self,
        vector: "np.array",
        cell_axis: "np.array",
        cell_length: float,
        G: "np.array",
    ) -> Tuple[float, float]:

        """Calculates the angle between a vector and a cell axis

        Requires the metric matrix because not cartesian coordinates

        Args:
            vector (np.array): vector to be compared
            cell_axis (np.array): unit cell vector to be compared
            cell_length (float): length of cell axis
            G (np.array): metric matrix of the unit cell

        Returns:
            angle (float): angle between unit cell vector and input vector
            supplementary_angle (float): supplementary angle of var 'angle'
        """

        dot_product = np.matmul(vector, np.matmul(G, cell_axis))
        div = dot_product / cell_length
        angle = math.acos(div) * (180 / math.pi)

        angle = round(angle, 1)

        supplementary_angle = 180 - angle

        return angle, supplementary_angle

    def analyse_data(self, csv_file: str, cell_data: str) -> None:

        """RECOMMENDED THAT CIF_READ.PY IS RUN FIRST FOR CORRECT FORMATTING

        OF INPUT .CSV FILES

        Will take information from ADP data and unit cell data and

        calculate the principle axes of all the atoms

        Also calculates the angle between these principle axes and the

        various unit cell axes

        Args:
            csv_file (str): path to csv containing ADP parameters
            cell_data (str): path to csv containing cell parameters
        """

        adp_df = pd.read_csv(csv_file)

        new_adp_df = pd.DataFrame()

        new_adp_df["Data_Block"] = adp_df["Data_Block"]
        new_adp_df["Atom"] = adp_df["_atom_site_aniso_label"]

        cell_df = pd.read_csv(cell_data)

        adp_by_cell = adp_df.groupby("Data_Block")

        counter = 0

        principle_A = []
        principle_B = []
        principle_C = []
        principle_A2 = []
        principle_B2 = []
        principle_C2 = []
        vector_A = []
        vector_B = []
        vector_C = []
        vector_A_angle_a = []
        vector_A_angle_b = []
        vector_A_angle_c = []
        vector_B_angle_a = []
        vector_B_angle_b = []
        vector_B_angle_c = []
        vector_C_angle_a = []
        vector_C_angle_b = []
        vector_C_angle_c = []
        vector_A_angle_a_sup = []
        vector_A_angle_b_sup = []
        vector_A_angle_c_sup = []
        vector_B_angle_a_sup = []
        vector_B_angle_b_sup = []
        vector_B_angle_c_sup = []
        vector_C_angle_a_sup = []
        vector_C_angle_b_sup = []
        vector_C_angle_c_sup = []

        # Loop through each CIF file

        for name, group in adp_by_cell:

            # Calculate metric matrix for the structure

            a = cell_df.iloc[counter]["_cell_length_a"]
            b = cell_df.iloc[counter]["_cell_length_b"]
            c = cell_df.iloc[counter]["_cell_length_c"]
            alpha = cell_df.iloc[counter]["_cell_angle_alpha"]
            beta = cell_df.iloc[counter]["_cell_angle_beta"]
            gamma = cell_df.iloc[counter]["_cell_angle_gamma"]

            alpha_rad = alpha * (math.pi / 180)
            beta_rad = beta * (math.pi / 180)
            gamma_rad = gamma * (math.pi / 180)

            G = np.array(
                [
                    [
                        (a * a),
                        (a * b * math.cos(gamma_rad)),
                        (a * c * math.cos(beta_rad)),
                    ],
                    [
                        (b * a * math.cos(gamma_rad)),
                        (b * b),
                        (b * c * math.cos(alpha_rad)),
                    ],
                    [
                        (c * a * math.cos(beta_rad)),
                        (c * b * math.cos(alpha_rad)),
                        (c * c),
                    ],
                ]
            )

            G_recip = np.linalg.inv(G)

            a_recip = math.sqrt(G_recip[0][0])
            b_recip = math.sqrt(G_recip[1][1])
            c_recip = math.sqrt(G_recip[2][2])

            # Loop through each atom

            for idx, row in group.iterrows():

                # Calculate the matricies of thermal parameters for each atom

                U11 = row["_atom_site_aniso_U_11"]
                U12 = row["_atom_site_aniso_U_12"]
                U13 = row["_atom_site_aniso_U_13"]
                U22 = row["_atom_site_aniso_U_22"]
                U23 = row["_atom_site_aniso_U_23"]
                U33 = row["_atom_site_aniso_U_33"]

                B11 = 2 * (math.pi**2) * U11 * (a_recip**2)
                B12 = 2 * (math.pi**2) * U12 * (a_recip * b_recip)
                B13 = 2 * (math.pi**2) * U13 * (a_recip * c_recip)
                B22 = 2 * (math.pi**2) * U22 * (b_recip**2)
                B23 = 2 * (math.pi**2) * U23 * (b_recip * c_recip)
                B33 = 2 * (math.pi**2) * U33 * (c_recip**2)

                B = np.array([[B11, B12, B13], [B12, B22, B23], [B13, B23, B33]])

                BG = np.matmul(B, G)

                values_unsorted, vectors = np.linalg.eig(BG)

                test = values_unsorted.tolist()

                values = sorted(test)

                values.reverse()

                positions = {}

                for i_1, i in enumerate(values_unsorted):
                    for i_2, j in enumerate(values):
                        if i == j:
                            positions[i_1] = i_2

                vector_1 = np.array(
                    [
                        [
                            vectors[0][positions[0]],
                            vectors[1][positions[0]],
                            vectors[2][positions[0]],
                        ]
                    ]
                )

                vector_2 = np.array(
                    [
                        [
                            vectors[0][positions[1]],
                            vectors[1][positions[1]],
                            vectors[2][positions[1]],
                        ]
                    ]
                )

                vector_3 = np.array(
                    [
                        [
                            vectors[0][positions[2]],
                            vectors[1][positions[2]],
                            vectors[2][positions[2]],
                        ]
                    ]
                )

                try:
                    principle_A.append(math.sqrt(values[0] / (2 * (math.pi**2))))
                except ValueError:
                    principle_A.append("NPD")
                try:
                    principle_B.append(math.sqrt(values[1] / (2 * (math.pi**2))))
                except ValueError:
                    principle_B.append("NPD")
                try:
                    principle_C.append(math.sqrt(values[2] / (2 * (math.pi**2))))
                except ValueError:
                    principle_C.append("NPD")

                principle_A2.append(values[0] / (2 * (math.pi**2)))
                principle_B2.append(values[1] / (2 * (math.pi**2)))
                principle_C2.append(values[2] / (2 * (math.pi**2)))

                scaled_vector_A = self.scale_vector(vector_1, G)
                scaled_vector_B = self.scale_vector(vector_2, G)
                scaled_vector_C = self.scale_vector(vector_3, G)

                vector_A.append(scaled_vector_A)
                vector_B.append(scaled_vector_B)
                vector_C.append(scaled_vector_C)

                a = np.array([[1], [0], [0]])
                b = np.array([[0], [1], [0]])
                c = np.array([[0], [0], [1]])

                v1_a, v1_a_sup = self.calculate_angle(scaled_vector_A, a, math.sqrt(G[0][0]), G)
                v1_b, v1_b_sup = self.calculate_angle(scaled_vector_A, b, math.sqrt(G[1][1]), G)
                v1_c, v1_c_sup = self.calculate_angle(scaled_vector_A, c, math.sqrt(G[2][2]), G)
                v2_a, v2_a_sup = self.calculate_angle(scaled_vector_B, a, math.sqrt(G[0][0]), G)
                v2_b, v2_b_sup = self.calculate_angle(scaled_vector_B, b, math.sqrt(G[1][1]), G)
                v2_c, v2_c_sup = self.calculate_angle(scaled_vector_B, c, math.sqrt(G[2][2]), G)
                v3_a, v3_a_sup = self.calculate_angle(scaled_vector_C, a, math.sqrt(G[0][0]), G)
                v3_b, v3_b_sup = self.calculate_angle(scaled_vector_C, b, math.sqrt(G[1][1]), G)
                v3_c, v3_c_sup = self.calculate_angle(scaled_vector_C, c, math.sqrt(G[2][2]), G)

                vector_A_angle_a.append(v1_a)
                vector_A_angle_b.append(v1_b)
                vector_A_angle_c.append(v1_c)
                vector_B_angle_a.append(v2_a)
                vector_B_angle_b.append(v2_b)
                vector_B_angle_c.append(v2_c)
                vector_C_angle_a.append(v3_a)
                vector_C_angle_b.append(v3_b)
                vector_C_angle_c.append(v3_c)

                vector_A_angle_a_sup.append(v1_a_sup)
                vector_A_angle_b_sup.append(v1_b_sup)
                vector_A_angle_c_sup.append(v1_c_sup)
                vector_B_angle_a_sup.append(v2_a_sup)
                vector_B_angle_b_sup.append(v2_b_sup)
                vector_B_angle_c_sup.append(v2_c_sup)
                vector_C_angle_a_sup.append(v3_a_sup)
                vector_C_angle_b_sup.append(v3_b_sup)
                vector_C_angle_c_sup.append(v3_c_sup)

            counter += 1

        new_adp_df["Root_Mean_Square_Principle_Axis_Displacement_1"] = principle_A
        new_adp_df["Root_Mean_Square_Principle_Axis_Displacement_2"] = principle_B
        new_adp_df["Root_Mean_Square_Principle_Axis_Displacement_3"] = principle_C
        new_adp_df["Mean_Square_Principle_Axis_Displacement_1"] = principle_A2
        new_adp_df["Mean_Square_Principle_Axis_Displacement_2"] = principle_B2
        new_adp_df["Mean_Square_Principle_Axis_Displacement_3"] = principle_C2
        new_adp_df["Principle_Vector_1"] = vector_A
        new_adp_df["Principle_Vector_2"] = vector_B
        new_adp_df["Principle_Vector_3"] = vector_C
        new_adp_df["Vector_1_Angle_to_a_axis"] = vector_A_angle_a
        new_adp_df["Supplementary_Angle_V1_a"] = vector_A_angle_a_sup
        new_adp_df["Vector_1_Angle_to_b_axis"] = vector_A_angle_b
        new_adp_df["Supplementary_Angle_V1_b"] = vector_A_angle_b_sup
        new_adp_df["Vector_1_Angle_to_c_axis"] = vector_A_angle_c
        new_adp_df["Supplementary_Angle_V1_c"] = vector_A_angle_c_sup
        new_adp_df["Vector_2_Angle_to_a_axis"] = vector_B_angle_a
        new_adp_df["Supplementary_Angle_V2_a"] = vector_B_angle_a_sup
        new_adp_df["Vector_2_Angle_to_b_axis"] = vector_B_angle_b
        new_adp_df["Supplementary_Angle_V2_b"] = vector_B_angle_b_sup
        new_adp_df["Vector_2_Angle_to_c_axis"] = vector_B_angle_c
        new_adp_df["Supplementary_Angle_V2_c"] = vector_B_angle_c_sup
        new_adp_df["Vector_3_Angle_to_a_axis"] = vector_C_angle_a
        new_adp_df["Supplementary_Angle_V3_a"] = vector_C_angle_a_sup
        new_adp_df["Vector_3_Angle_to_b_axis"] = vector_C_angle_b
        new_adp_df["Supplementary_Angle_V3_b"] = vector_C_angle_b_sup
        new_adp_df["Vector_3_Angle_to_c_axis"] = vector_C_angle_c
        new_adp_df["Supplementary_Angle_V3_c"] = vector_C_angle_c_sup

        adp_by_atom = new_adp_df.groupby("Atom")

        os.chdir(pathlib.Path(csv_file).parent)

        new_folder = "Individual_Atomic_ADP_Analysis"

        try:
            os.mkdir(new_folder)
        except FileExistsError:
            pass

        for i, j in adp_by_atom:

            csv_location = pathlib.Path(csv_file).parent / new_folder / (str(i) + "_ADPs.csv")

            j.to_csv(csv_location, index=False)
