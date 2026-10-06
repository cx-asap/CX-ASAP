#!/usr/bin/env python3

###################################################################################################
# ----------------------------------------CX-ASAP: TEMPLATE----------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config, File_Sorter
from data_reduction.modules.xprep_intensity_compare import Intensity_Compare
import os
import pandas as pd
import logging

# ----------Class Definition----------#


class Intensity_Pipeline:
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

    def multiple_intensity(
        self,
        location: str,
        h_condition_1: str,
        k_condition_1: str,
        l_condition_1: str,
        h_condition_2: str,
        k_condition_2: str,
        l_condition_2: str,
        hk_condition_1: str,
        hk_condition_2: str,
        hl_condition_1: str,
        hl_condition_2: str,
        kl_condition_1: str,
        kl_condition_2: str,
        hkl_condition_1: str,
        hkl_condition_2: str,
        a_axis: str,
        b_axis: str,
        c_axis: str,
        alpha: str,
        beta: str,
        gamma: str,
        h1_0: bool,
        k1_0: bool,
        l1_0: bool,
        h2_0: bool,
        k2_0: bool,
        l2_0: bool,
    ) -> None:

        """Goes through each .hkl file in a specified folder and runs the

        xprep_intensity_module on each

        This allows for examination of the differences between two criteria

        over a set of files to see if additional changes are occurring

        Two combined data frames are output:

        - full

        - simplified

        Args:
            location (str): Path to the folder containing a range of reflection files
            h_condition_1 (str): Group 1 condition for h
            k_condition_1 (str): Group 1 condition for k
            l_condition_1 (str): Group 1 condition for l
            h_condition_2 (str): Group 2 condition for h
            k_condition_2 (str): Group 2 condition for k
            l_condition_2 (str): Group 2 condition for l
            hk_condition_1 (str): Group 1 condition for h+k
            hk_condition_2 (str): Group 2 condition for h+k
            hl_condition_1 (str): Group 1 condition for h+l
            hl_condition_2 (str): Group 2 condition for h+l
            kl_condition_1 (str): Group 1 condition for k+l
            kl_condition_2 (str): Group 2 condition for k+l
            hkl_condition_1 (str): Group 1 condition for h+k+l
            hkl_condition_2 (str): Group 2 condition for h+k+l
            file_name (str): full path to the reflection file for analysis
            a_axis (str): a_axis of the unit cell
            b_axis (str): b_axis of the unit cell
            c_axis (str): c_axis of the unit cell
            alpha (str): alpha angle of the unit cell
            beta (str): beta angle of the unit cell
            gamma (str): gamma angle of the unit cell
            h1_0 (bool): Group 1 - include reflections with h=0? T or F
            k1_0 (bool): Group 1 - include reflections with k=0? T or F
            l1_0 (bool): Group 1 - include reflections with l=0? T or F
            h2_0 (bool): Group 2 - include reflections with h=0? T or F
            k2_0 (bool): Group 2 - include reflections with k=0? T or F
            l2_0 (bool): Group 2 - include reflections with l=0? T or F


        """

        os.chdir(location)

        sorter = File_Sorter()

        files = sorter.sorted_properly(os.listdir())

        hkl_files = []

        for item in files:
            if item.endswith(".hkl"):
                hkl_files.append(item)

        reflections = Intensity_Compare(self.test_mode)

        smallest_d_spacing, largest_d_spacing = reflections.analyse_reflections(
            h_condition_1,
            k_condition_1,
            l_condition_1,
            h_condition_2,
            k_condition_2,
            l_condition_2,
            hk_condition_1,
            hk_condition_2,
            hl_condition_1,
            hl_condition_2,
            kl_condition_1,
            kl_condition_2,
            hkl_condition_1,
            hkl_condition_2,
            hkl_files[0],
            a_axis,
            b_axis,
            c_axis,
            alpha,
            beta,
            gamma,
            h1_0,
            k1_0,
            l1_0,
            h2_0,
            k2_0,
            l2_0,
        )

        for item in hkl_files:

            print("Analysing " + str(item))

            reflections.analyse_reflections(
                h_condition_1,
                k_condition_1,
                l_condition_1,
                h_condition_2,
                k_condition_2,
                l_condition_2,
                hk_condition_1,
                hk_condition_2,
                hl_condition_1,
                hl_condition_2,
                kl_condition_1,
                kl_condition_2,
                hkl_condition_1,
                hkl_condition_2,
                item,
                a_axis,
                b_axis,
                c_axis,
                alpha,
                beta,
                gamma,
                h1_0,
                k1_0,
                l1_0,
                h2_0,
                k2_0,
                l2_0,
                smallest_d_spacing,
                largest_d_spacing,
            )

        combined = pd.DataFrame()

        I1 = []
        R1 = []
        IS1 = []
        I2 = []
        R2 = []
        IS2 = []
        file_names = []

        if os.path.exists("Combined.csv"):
            os.remove("Combined.csv")
        if os.path.exists("Simplified_Combined.csv"):
            os.remove("Simplified_Combined.csv")

        for item in sorter.sorted_properly(os.listdir()):
            if item.endswith(".csv"):
                temp = pd.read_csv(item)
                try:
                    combined = pd.merge(combined, temp, how="outer", on=["d-spacing"])
                except KeyError:
                    combined = temp

                all_columns = temp.columns.values.tolist()
                needed_columns = [
                    "Overall Average Intensity Group 1",
                    "Overall Rsigma Group 1",
                    "I/sigma Group 1",
                    "Overall Average Intensity Group 2",
                    "Overall Rsigma Group 2",
                    "I/sigma Group 2",
                ]
                extracted_columns = []

                for i in all_columns:
                    for j in needed_columns:
                        if j in i:
                            extracted_columns.append(i)

                I1.append(temp[extracted_columns[0]][0])
                R1.append(temp[extracted_columns[1]][0])
                IS1.append(temp[extracted_columns[2]][0])
                I2.append(temp[extracted_columns[3]][0])
                R2.append(temp[extracted_columns[4]][0])
                IS2.append(temp[extracted_columns[5]][0])
                file_names.append(item)

        data = {
            "HKL File": file_names,
            "Group 1 Average Intensities": I1,
            "Group 1 Average Rsigma": R1,
            "Group 1 I/sigma": IS1,
            "Group 2 Average Intensities": I2,
            "Group 2 Average Rsigma": R2,
            "Group 2 I/sigma": IS2,
        }

        combined_easy = pd.DataFrame(data)

        data = combined.to_csv("Combined.csv", index=False)
        data2 = combined_easy.to_csv("Simplified_Combined.csv", index=False)
