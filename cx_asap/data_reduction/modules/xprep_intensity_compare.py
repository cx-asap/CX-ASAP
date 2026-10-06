#!/usr/bin/env python3

###################################################################################################
# --------------------------------CX-ASAP: xprep intensity comparison------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config
import os
import subprocess
from numpy import mean
import numpy as np
import math
import pandas as pd
import pathlib
import logging
from typing import Tuple

# ----------Class Definition----------#


class Intensity_Compare:
    def __init__(self, test_mode=False) -> None:

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

    def analyse_condition(self, condition: str) -> Tuple[int, int, int, str]:

        """This function analyses the reflection condition input by the user

        and separates it into a form the code can understand

        Args:
            conditon (str): The reflection condition input by the user
                            Example: "2n+1"

        Returns:
            number_to_be_multiplied_by_n (int): as the name suggests (ie "2")
            lhs (int): left hand side of the operator without "n" (ie "2")
            rhs (int): right hand side of the operator without "n" (ie "1")
            operator (str): either "+" or "-"
        """

        l_or_r = -1

        if condition == "n":
            number_to_be_multiplied_by_n = False
            lhs = False
            rhs = False
            operator = False
        else:

            # for subtracting two terms where one is 'n'

            if "-" in str(condition):
                test = condition.split("-")
                for index, item in enumerate(test):
                    if "n" in item:
                        item.strip("n")
                        l_or_r = index

                lhs = int(test[0].strip("n "))
                rhs = int(test[1].strip("n "))
                operator = "-"

                if l_or_r == 1:
                    number_to_be_multiplied_by_n = rhs
                elif l_or_r == 0:
                    number_to_be_multiplied_by_n = lhs
                else:
                    number_to_be_multiplied_by_n = 0

            # for adding two terms where one is 'n'

            elif "+" in str(condition):
                test = condition.split("+")
                for index, item in enumerate(test):
                    if "n" in item:
                        item.strip("n")
                        l_or_r = index

                lhs = int(test[0].strip("n "))
                rhs = int(test[1].strip("n "))
                operator = "+"

                if l_or_r == 1:
                    number_to_be_multiplied_by_n = rhs
                elif l_or_r == 0:
                    number_to_be_multiplied_by_n = lhs
                else:
                    number_to_be_multiplied_by_n = 0

            # for single term conditions (both with and without n)

            else:
                rhs = False
                operator = False
                if "n" in str(condition):
                    number_to_be_multiplied_by_n = int(condition.strip("n"))
                    lhs = False

                else:
                    number_to_be_multiplied_by_n = False
                    lhs = int(condition)

        return number_to_be_multiplied_by_n, lhs, rhs, operator

    def math_doer(
        self, number_to_be_multiplied_by_n: int, lhs: int, rhs: int, operator: str
    ) -> list:

        """This function returns a set of values from -40 -> 40 that

        match the reflection conditions input by the user

        Args:
            number_to_be_multiplied_by_n (int): as the name suggests (ie "2")
            lhs (int): left hand side of the operator without "n" (ie "2")
            rhs (int): right hand side of the operator without "n" (ie "1")
            operator (str): either "+" or "-"

        Returns:
            numbers_that_match_condition (list): a list of numbers that
                                                match the condition

        """

        n_values = range(-40, 40, 1)
        numbers_that_match_condition = []
        for n in n_values:
            if number_to_be_multiplied_by_n == lhs and operator == "+":
                value = (number_to_be_multiplied_by_n * n) + rhs
            elif number_to_be_multiplied_by_n == lhs and operator == "-":
                value = (number_to_be_multiplied_by_n * n) - rhs
            elif number_to_be_multiplied_by_n == rhs and operator == "+":
                value = lhs + (number_to_be_multiplied_by_n * n)
            elif number_to_be_multiplied_by_n == rhs and operator == "-":
                value = lhs - (number_to_be_multiplied_by_n * n)
            elif number_to_be_multiplied_by_n == 0 and operator == "+":
                value = lhs + rhs
            elif number_to_be_multiplied_by_n == 0 and operator == "-":
                value = lhs - rhs
            elif (
                operator == False
                and lhs == False
                and number_to_be_multiplied_by_n != False
            ):
                value = number_to_be_multiplied_by_n * n
            elif operator == False and lhs != False:
                value = lhs
            elif operator == False and lhs == False and rhs == False:
                value = n
            numbers_that_match_condition.append(str(value))
        return numbers_that_match_condition

    def metric_matrix(self, unit_cell: list) -> "np.array":

        """This function calculates the metric matrix for a given unit cell

        Args:
            unit_cell (list): the 6 unit cell parameters in a list
                                (ie [a, b, c, alpha, beta, gamma])

        Returns:
            metric_unit_cell (np.array): the metrix matrix of the unit cell

        """

        metric_unit_cell = np.array(
            [
                [
                    unit_cell[0] * unit_cell[0],
                    unit_cell[0]
                    * unit_cell[1]
                    * math.cos(unit_cell[5] * (math.pi / 180)),
                    unit_cell[0]
                    * unit_cell[2]
                    * math.cos(unit_cell[4] * (math.pi / 180)),
                ],
                [
                    unit_cell[1]
                    * unit_cell[0]
                    * math.cos(unit_cell[5] * (math.pi / 180)),
                    unit_cell[1] * unit_cell[1],
                    unit_cell[1]
                    * unit_cell[2]
                    * math.cos(unit_cell[3] * (math.pi / 180)),
                ],
                [
                    unit_cell[2]
                    * unit_cell[0]
                    * math.cos(unit_cell[4] * (math.pi / 180)),
                    unit_cell[2]
                    * unit_cell[1]
                    * math.cos(unit_cell[3] * (math.pi / 180)),
                    unit_cell[2] * unit_cell[2],
                ],
            ]
        )

        return metric_unit_cell

    def d_spacing(self, h: int, k: int, l: int, reciprocal_matrix: "np.array") -> float:

        """Calculates the d_spacing of a reflection

        using the reciprocoal of the metric matrix

        Args:
            h (int): 'h' of the hkl reflection
            k (int): 'k' of the hkl reflection
            l (int): 'l' of the hkl reflection
            reciprocal_matrix (np.array): reciprocal of the metrix matrix

        Returns:
            d_spacing (float): the calculated d_spacing of the reflection

        """

        try:
            d_spacing = 1 / math.sqrt(
                h * h * reciprocal_matrix[0][0]
                + h * k * reciprocal_matrix[0][1]
                + h * l * reciprocal_matrix[0][2]
                + k * h * reciprocal_matrix[1][0]
                + k * k * reciprocal_matrix[1][1]
                + k * l * reciprocal_matrix[1][2]
                + l * h * reciprocal_matrix[2][0]
                + l * k * reciprocal_matrix[2][1]
                + l * l * reciprocal_matrix[2][2]
            )
        except ZeroDivisionError:
            d_spacing = 0

        return d_spacing

    def break_into_resolution_shell(
        self,
        a_axis: str,
        b_axis: str,
        c_axis: str,
        alpha: str,
        beta: str,
        gamma: str,
        reflections_1: dict,
        reflections_2: dict,
        file_name: str,
        smallest_d_spacing: float,
        largest_d_spacing: float,
    ) -> Tuple[float, float]:

        """Sorts the reflections into resolution shells and

        averages the intensities of each bin.

        What the bins are is calculated based on the smallest

        and largest d_spacing present in the reflection file

        Args:
            a_axis (str): a_axis of the unit cell
            b_axis (str): b_axis of the unit cell
            c_axis (str): c_axis of the unit cell
            alpha (str): alpha angle of the unit cell
            beta (str): beta angle of the unit cell
            gamma (str): gamma angle of the unit cell
            reflections_1 (dict): reflections from the first condition to compare
            reflections_2 (dict): reflections from the second condition to compare
            file_name (str): name of the reflection file
            smallest_d_spacing (float): smallest d_spacing of any reflection
            largest_d_spacing (float): largest d_spacing of any reflection

        Returns:
            smallest_d_spacing (float): the smallest calculated d_spacing of the reflection
            largest_d_spacing (float): the largest calculated d_spacing of the reflection

        """

        overall_sum_1_intensities = 0

        overall_sum_2_intensities = 0

        overall_sum_1_errors = 0

        overall_sum_2_errors = 0

        intensities1 = pd.DataFrame(columns=["d_spacing", "intensities", "sigma"])
        intensities2 = pd.DataFrame(columns=["d_spacing", "intensities", "sigma"])

        unit_cell = [
            float(a_axis),
            float(b_axis),
            float(c_axis),
            float(alpha),
            float(beta),
            float(gamma),
        ]

        metric_cell = self.metric_matrix(unit_cell)

        reciprocal_cell = np.linalg.inv(metric_cell)

        if reflections_1 != False:

            for item in reflections_1:
                d_spacing = self.d_spacing(
                    float(reflections_1[item][0]),
                    float(reflections_1[item][1]),
                    float(reflections_1[item][2]),
                    reciprocal_cell,
                )
                if d_spacing != 0:
                    new_row = {
                        "d_spacing": d_spacing,
                        "intensities": float(reflections_1[item][3]),
                        "sigma": float(reflections_1[item][4]),
                    }
                    intensities1 = intensities1.append(new_row, ignore_index=True)

                    overall_sum_1_intensities += float(reflections_1[item][3])
                    overall_sum_1_errors += float(reflections_1[item][4])

                    final_number_of_reflns_1 = int(item)

            df_smallest1 = intensities1.nsmallest(1, ["d_spacing"])
            df_largest1 = intensities1.nlargest(1, ["d_spacing"])
            smallest_d_spacing1 = df_smallest1.iloc[0]["d_spacing"]
            largest_d_spacing1 = df_largest1.iloc[0]["d_spacing"]

        if reflections_2 != False:

            for item in reflections_2:
                d_spacing = self.d_spacing(
                    float(reflections_2[item][0]),
                    float(reflections_2[item][1]),
                    float(reflections_2[item][2]),
                    reciprocal_cell,
                )
                if d_spacing != 0:
                    new_row = {
                        "d_spacing": d_spacing,
                        "intensities": float(reflections_2[item][3]),
                        "sigma": float(reflections_2[item][4]),
                    }
                    intensities2 = intensities2.append(new_row, ignore_index=True)

                    overall_sum_2_intensities += float(reflections_2[item][3])
                    overall_sum_2_errors += float(reflections_2[item][4])

                    final_number_of_reflns_2 = int(item)

            df_smallest2 = intensities2.nsmallest(1, ["d_spacing"])
            df_largest2 = intensities2.nlargest(1, ["d_spacing"])

            smallest_d_spacing2 = df_smallest2.iloc[0]["d_spacing"]
            largest_d_spacing2 = df_largest2.iloc[0]["d_spacing"]

        average_intensity_1 = overall_sum_1_intensities / final_number_of_reflns_1
        average_intensity_2 = overall_sum_2_intensities / final_number_of_reflns_2

        overall_rsigma_1 = overall_sum_1_errors / overall_sum_1_intensities
        overall_rsigma_2 = overall_sum_2_errors / overall_sum_2_intensities

        i_on_sigma_1 = 1 / overall_rsigma_1
        i_on_sigma_2 = 1 / overall_rsigma_2

        if smallest_d_spacing == False or largest_d_spacing == False:

            if reflections_1 != False and reflections_2 != False:
                if smallest_d_spacing2 >= smallest_d_spacing1:
                    smallest_d_spacing = smallest_d_spacing1
                else:
                    smallest_d_spacing = smallest_d_spacing2

                if largest_d_spacing2 >= largest_d_spacing1:
                    largest_d_spacing = largest_d_spacing2
                else:
                    largest_d_spacing = largest_d_spacing1

            elif reflections_1 == False:
                largest_d_spacing = largest_d_spacing2
                smallest_d_spacing = smallest_d_spacing2

            elif reflections_2 == False:
                largest_d_spacing = largest_d_spacing1
                smallest_d_spacing = smallest_d_spacing1

        def custom_output(x) -> "pd.Series":

            """Firstly I'm sorry this is a function within a function :(

            This was a pandas thing that I was doing to try and make everything neater.

            Am now confused as to what it does, but it was important for making

            the getting the right data in the output dataframe

            Args:
                x(???): ???

            """

            d = {}
            d["intensities"] = x["intensities"].mean()
            try:
                d["Rsigma"] = sum(x["sigma"]) / sum(x["intensities"])
            except ZeroDivisionError:
                d["Rsigma"] = 0

            d["number of reflections"] = len(x["d_spacing"])

            return pd.Series(
                d, index=["number of reflections", "intensities", "Rsigma"]
            )

        if reflections_1 != False:

            intensities1 = intensities1.groupby(
                pd.cut(
                    intensities1["d_spacing"],
                    np.arange(smallest_d_spacing, largest_d_spacing, 0.1),
                ),
                as_index=False,
            ).apply(custom_output)

        if reflections_2 != False:

            intensities2 = intensities2.groupby(
                pd.cut(
                    intensities2["d_spacing"],
                    np.arange(smallest_d_spacing, largest_d_spacing, 0.1),
                ),
                as_index=False,
            ).apply(custom_output)

        test = np.arange(smallest_d_spacing, largest_d_spacing, 0.1)
        test = test.round(decimals=2)
        test2 = np.arange(smallest_d_spacing + 0.1, largest_d_spacing + 0.1, 0.1)
        test2 = test2.round(decimals=2)

        ranges = []

        for index, item in enumerate(test):
            ranges.append(str(test[index]) + " - " + str(test2[index]))

        intensities = pd.DataFrame()

        intensities["d-spacing"] = ranges

        if reflections_1 != False:
            intensities[
                "Average Intensity Group 1 " + pathlib.Path(file_name).stem
            ] = intensities1["intensities"]

            intensities[
                "Rsigma Group 1 " + pathlib.Path(file_name).stem
            ] = intensities1["Rsigma"]

            intensities[
                "Number of reflections Group 1 " + pathlib.Path(file_name).stem
            ] = intensities1["number of reflections"]

        if reflections_2 != False:
            intensities[
                "Average Intensity Group 2 " + pathlib.Path(file_name).stem
            ] = intensities2["intensities"]

            intensities[
                "Rsigma Group 2 " + pathlib.Path(file_name).stem
            ] = intensities2["Rsigma"]

            intensities[
                "Number of reflections Group 2 " + pathlib.Path(file_name).stem
            ] = intensities2["number of reflections"]

        if reflections_1 != False and reflections_2 != False:

            ratios = []
            biggest_intensity = 0
            smallest_intensity = 0

            for index, item in enumerate(ranges):
                if (
                    intensities.iloc[index][
                        "Average Intensity Group 1 " + pathlib.Path(file_name).stem
                    ]
                    >= intensities.iloc[index][
                        "Average Intensity Group 2 " + pathlib.Path(file_name).stem
                    ]
                ):
                    biggest_intensity = intensities.iloc[index][
                        "Average Intensity Group 1 " + pathlib.Path(file_name).stem
                    ]
                    smallest_intensity = intensities.iloc[index][
                        "Average Intensity Group 2 " + pathlib.Path(file_name).stem
                    ]
                    ratios.append(
                        str(round(biggest_intensity / smallest_intensity, 2))
                        + " : "
                        + str(round(smallest_intensity / smallest_intensity, 2))
                    )
                else:
                    biggest_intensity = intensities.iloc[index][
                        "Average Intensity Group 2 " + pathlib.Path(file_name).stem
                    ]
                    smallest_intensity = intensities.iloc[index][
                        "Average Intensity Group 1 " + pathlib.Path(file_name).stem
                    ]
                    ratios.append(
                        str(round(smallest_intensity / smallest_intensity, 2))
                        + " : "
                        + str(round(biggest_intensity / smallest_intensity, 2))
                    )

            intensities["Intensity Ratio " + pathlib.Path(file_name).stem] = ratios

        n = len(ranges) - 1

        aaa = [average_intensity_1]
        aaa += n * [" "]
        bbb = [overall_rsigma_1]
        bbb += n * [" "]
        ccc = [i_on_sigma_1]
        ccc += n * [" "]
        ddd = [average_intensity_2]
        ddd += n * [" "]
        eee = [overall_rsigma_2]
        eee += n * [" "]
        fff = [i_on_sigma_2]
        fff += n * [" "]

        intensities["Overall Average Intensity Group 1"] = aaa
        intensities["Overall Rsigma Group 1"] = bbb
        intensities["I/sigma Group 1"] = ccc
        intensities["Overall Average Intensity Group 2"] = ddd
        intensities["Overall Rsigma Group 2"] = eee
        intensities["I/sigma Group 2"] = fff

        os.chdir(pathlib.Path(file_name).parent)

        data = intensities.to_csv(pathlib.Path(file_name).stem + ".csv", index=False)

        return smallest_d_spacing, largest_d_spacing

    def analyse_reflections(
        self,
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
        file_name: str,
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
        small_d: bool = False,
        big_d: bool = False,
    ):

        """Compares the intensities of two criteria of reflections in a file

        For example, it analyses the intensity differences between

        even ("2n") and odd ("2n+1") reflections in a single file

        It does this by using the above functions to do a range of calculations

        Args:
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
            small_d(bool):=False : Iterative to have smallest_d_spacing (can also be str)
            big_d:bool=False : Iterative to have largest_d_spacing (can also be str)

        Returns:
            smallest_d_spacing (float): the smallest calculated d_spacing of the reflection
            largest_d_spacing (float): the largest calculated d_spacing of the reflection

        """

        reflections_1 = {}
        reflections_2 = {}

        with open(file_name, "rt") as f:
            flag = False
            for index, line in enumerate(f):

                refln = []

                for item in line.split(" "):
                    try:
                        float(item.strip("\n"))
                    except:
                        pass
                    else:
                        refln.append(item.strip("\n"))

                if len(refln) != 0 and flag == False:

                    reflections_1[str(index)] = refln
                    reflections_2[str(index)] = refln

                try:
                    test = refln[2]
                except IndexError:
                    flag = True
                else:
                    if refln[0] == "0" and refln[1] == "0" and refln[2] == "0":
                        flag = True

        a, b, c, d = self.analyse_condition(h_condition_1)

        e = self.math_doer(a, b, c, d)

        f, g, h, i = self.analyse_condition(h_condition_2)

        j = self.math_doer(f, g, h, i)

        k, l, m, n = self.analyse_condition(k_condition_1)

        o = self.math_doer(k, l, m, n)

        p, q, r, s = self.analyse_condition(k_condition_2)

        t = self.math_doer(p, q, r, s)

        u, v, w, x = self.analyse_condition(l_condition_1)

        y = self.math_doer(u, v, w, x)

        z, aa, bb, cc = self.analyse_condition(l_condition_2)

        dd = self.math_doer(z, aa, bb, cc)

        ### Combination conditions ###

        ee, ff, gg, hh = self.analyse_condition(hk_condition_1)

        ii = self.math_doer(ee, ff, gg, hh)

        jj, kk, ll, mm = self.analyse_condition(hk_condition_2)

        nn = self.math_doer(jj, kk, ll, mm)

        oo, pp, qq, rr = self.analyse_condition(hl_condition_1)

        ss = self.math_doer(oo, pp, qq, rr)

        tt, uu, vv, ww = self.analyse_condition(hl_condition_2)

        xx = self.math_doer(tt, uu, vv, ww)

        yy, zz, aaa, bbb = self.analyse_condition(kl_condition_1)

        ccc = self.math_doer(yy, zz, aaa, bbb)

        ddd, eee, fff, ggg = self.analyse_condition(kl_condition_2)

        hhh = self.math_doer(ddd, eee, fff, ggg)

        iii, jjj, kkk, lll = self.analyse_condition(hkl_condition_1)

        mmm = self.math_doer(iii, jjj, kkk, lll)

        nnn, ooo, ppp, qqq = self.analyse_condition(hkl_condition_2)

        rrr = self.math_doer(nnn, ooo, ppp, qqq)

        ones_to_delete_1 = []

        ones_to_delete_2 = []

        for item in reflections_1:

            # Find values that do not satisfy the conditions for group 1

            if (
                reflections_1[item][0] not in e
                or reflections_1[item][1] not in o
                or reflections_1[item][2] not in y
                or (str(int(reflections_1[item][0]) + int(reflections_1[item][1])))
                not in ii
                or (str(int(reflections_1[item][0]) + int(reflections_1[item][2])))
                not in ss
                or (str(int(reflections_1[item][1]) + int(reflections_1[item][2])))
                not in ccc
                or (
                    str(
                        int(reflections_1[item][0])
                        + int(reflections_1[item][1])
                        + int(reflections_1[item][2])
                    )
                )
                not in mmm
            ):
                ones_to_delete_1.append(item)

            # Find values that do not satisfy the conditions for group 2

            if (
                reflections_2[item][0] not in j
                or reflections_2[item][1] not in t
                or reflections_2[item][2] not in dd
                or (str(int(reflections_2[item][0]) + int(reflections_2[item][1])))
                not in nn
                or (str(int(reflections_2[item][0]) + int(reflections_2[item][2])))
                not in xx
                or (str(int(reflections_2[item][1]) + int(reflections_2[item][2])))
                not in hhh
                or (
                    str(
                        int(reflections_2[item][0])
                        + int(reflections_2[item][1])
                        + int(reflections_2[item][2])
                    )
                )
                not in rrr
            ):
                ones_to_delete_2.append(item)

            # See if any reflections containing '0' need to be removed

            if h1_0 == False and reflections_1[item][0] == "0":
                ones_to_delete_1.append(item)

            if h2_0 == False and reflections_2[item][0] == "0":
                ones_to_delete_2.append(item)

            if k1_0 == False and reflections_1[item][1] == "0":
                ones_to_delete_1.append(item)

            if k2_0 == False and reflections_2[item][1] == "0":
                ones_to_delete_2.append(item)

            if l1_0 == False and reflections_1[item][2] == "0":
                ones_to_delete_1.append(item)

            if l2_0 == False and reflections_2[item][2] == "0":
                ones_to_delete_2.append(item)

        # Remove duplicates (ie from reflections with multiple 0's)

        ones_to_delete_1 = list(dict.fromkeys(ones_to_delete_1))
        ones_to_delete_2 = list(dict.fromkeys(ones_to_delete_2))

        for item in ones_to_delete_1:
            del reflections_1[item]

        for item in ones_to_delete_2:
            del reflections_2[item]

        if len(reflections_2) != 0 and len(reflections_1) != 0:

            ###Break up by resolution shells

            smallest_d_spacing, largest_d_spacing = self.break_into_resolution_shell(
                a_axis,
                b_axis,
                c_axis,
                alpha,
                beta,
                gamma,
                reflections_1,
                reflections_2,
                file_name,
                small_d,
                big_d,
            )

        elif len(reflections_1) == 0 and len(reflections_2) != 0:
            print("Criteria 1 is fully absent")
            smallest_d_spacing, largest_d_spacing = self.break_into_resolution_shell(
                a_axis,
                b_axis,
                c_axis,
                alpha,
                beta,
                gamma,
                False,
                reflections_2,
                file_name,
                small_d,
                big_d,
            )
        elif len(reflections_2) == 0 and len(reflections_1) != 0:
            print("Criteria 2 is fully absent")
            smallest_d_spacing, largest_d_spacing = self.break_into_resolution_shell(
                a_axis,
                b_axis,
                c_axis,
                alpha,
                beta,
                gamma,
                reflections_1,
                False,
                file_name,
                small_d,
                big_d,
            )
        else:
            print("Both criteria are fully absent - no analysis performed")

        return smallest_d_spacing, largest_d_spacing
