#!/usr/bin/env python3

###################################################################################################
# ---------------------------------------CX-ASAP: refinement---------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import (
    Nice_YAML_Dumper,
    Config,
    Directory_Browse,
    File_Sorter,
    Grapher,
)

import os
import re
import statistics
import pandas as pd
import pathlib
import subprocess
import shutil
import logging
from typing import Tuple

# ----------Class Definition----------#


class Structure_Refinement:
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

    def merge_data(self, structure: str, cell: str) -> str:
        """Adds reference structure data into a new .ins file

        From new file, take TITL->LATT

        From reference file, take LATT->END

        Args:
            structure (str): data from new structure file
            cell (str): data from reference file

        Returns:
            complete_file(str): combined file
        """

        ref_x = re.search("LATT", structure)
        ref_y = re.search("END", structure)

        new_x = re.search("TITL", cell)
        new_y = re.search("LATT", cell)

        complete_file = ""

        if new_x is None or new_y is None or ref_x is None or ref_y is None:
            if os.path.exists(self.reference) == True:
                pass
            else:
                logging.critical(__name__ + " : .ins files of incorrect format")
                print("Error in creating .ins files - check error log")
                exit()
        else:
            complete_file = (
                cell[new_x.start() : new_y.start()]
                + structure[ref_x.start() : ref_y.end()]
            )

        return complete_file

    def import_refinement(self, file_name: str, ref_struct: str) -> None:

        """Imports the reference and the new .ins files

        Outputs the combined file

        Args:
            file_name (str): name of the new structure file
            ref_struct (str): full path to the reference structure
        """

        self.reference = ref_struct

        with open(self.reference, "rt") as reference:
            structure = reference.read()

        try:
            with open(file_name, "rt") as new_file:
                cell = new_file.read()
        except FileNotFoundError:
            logging.info(__name__ + " : Missing .ins file")
            print("Error - See Error Log for more info")
        else:
            with open(file_name, "rt") as new_file:
                cell = new_file.read()

            merged = self.merge_data(structure, cell)

            with open(file_name, "w") as combined:
                for line in merged:
                    combined.write(line)

    def converge(
        self,
        convergence: bool,
        refinement_failed: bool,
        shift: list,
        refinements: int,
        tolerance: float,
        weights_old: list,
        weights_new: list,
    ) -> Tuple[bool, bool]:

        """Check if the structure has converged or not

        For the structure to have converged, it requires:

        1) That the weighting scheme has remained unchanged between refinements

        2) The maximum number of refinements allowed has occurred

        3) The average shift of the last "X" refinements is smaller than "Y"

                "X" = refinements
                "Y" = tolerance

        Args:
            convergence (bool): whether or not the refinement is converged
            refinement_failed (bool): whether or not the refinement has failed
            shift (list): a list of previous shifts to be appended to
            refinements (int): number of refinements to check for shift convergence
            tolerance (float): target average shift
            weights_old (list): weight values from previous refinement
            weights_new (list): weight values from current refinement


        Returns:
            convergence (bool): whether or not the refinement has converged
            refinement_failed (bool): whether or not the refinement failed
        """

        try:
            statistics.mean(shift)
        except statistics.StatisticsError:
            logging.info(__name__ + " : No atoms for refinement")
            convergence = True
        else:

            if abs(statistics.mean(shift[-int(refinements) :])) <= float(tolerance):
                if str(weights_old[0]).strip("0.") == str(weights_new[0]).strip(
                    "0."
                ) and str(weights_old[1]).strip("0.") == str(weights_new[1]).strip(
                    "0."
                ):
                    convergence = True
                    logging.info(__name__ + " : Refinement has converged")
            else:
                convergence = False
                logging.info(__name__ + " : Refinement has not converged")

        if len(weights_new) != 2:
            logging.info(__name__ + " : Refinement has failed")
            convergence = True
            refinement_failed = True
        else:
            logging.info(__name__ + " : Stats check!")
            logging.info(
                __name__
                + " : mean shift is: "
                + str(statistics.mean(shift[-int(refinements) :]))
            )
            logging.info(
                __name__
                + " : old weight is: "
                + str(weights_old[0])
                + " "
                + str(weights_old[1])
            )
            logging.info(
                __name__
                + " : new weight is: "
                + str(weights_new[0])
                + " "
                + str(weights_new[1])
            )

        return convergence, refinement_failed

    def convergence_check(
        self,
        lines: str,
        shift: list,
        old_weight: str,
        new_weight: str,
        refinements: int,
        tolerance: float,
        refinement_failed: bool,
    ) -> Tuple[bool, list, bool]:

        """Function sets up values to check the convergence of the refinements

        The actual convergence check is in a separate function (above).

        Args:
            input_file (str): name of the file to be checked
            shift (list): a list of previous shifts to be appended to
            old_weight (str): weighting scheme of the previous refinement
            new_weight (str): weighting scheme of the current refinement
            refinements (int): number of refinements to check for shift convergence
            tolerance (float): target average shift
            refinement_failed (bool): whether or not the refinement has failed

        Returns:
            convergence (bool): whether or not the refinement has converged
            shift (list): updated list of shift parameters
            refinement_failed (bool): whether or not the refinement failed
        """

        convergence = False

        # Investigate shifts

        shift_param = []

        for line in lines:
            if "Mean shift" in line:
                shift_param.append(line)

        for item in shift_param:
            try:
                float(item.split(" ")[6])
            except ValueError:
                logging.info(__name__ + " : Structure likely exploded")
                shift.append(float(99))
            else:
                shift.append(float(item.split(" ")[6]))

        # Investigate weights

        weights_old = []
        weights_new = []

        for item in old_weight.split(" "):
            if item != "WGHT" and item != "":
                try:
                    float(item.strip("\n"))
                except:
                    pass
                else:
                    weights_old.append(item.strip("\n"))
        for item in new_weight.split(" "):
            if item != "WGHT" and item != "":
                try:
                    float(item.strip("\n"))
                except:
                    pass
                else:
                    weights_new.append(item.strip("\n"))

        if len(weights_old) == 1:
            weights_old.append(0)
        if len(weights_new) == 1:
            weights_new.append(0)

        convergence, refinement_failed = self.converge(
            convergence,
            refinement_failed,
            shift,
            refinements,
            tolerance,
            weights_old,
            weights_new,
        )

        return convergence, shift, refinement_failed

    def run_shelxl(
        self,
        ins_file: str,
        reference: str,
        refinements: int,
        tolerance: float,
        max_cycles: int,
    ) -> bool:

        """Runs SHELXL on a single structure that has had a reference model

        imported into it

        Will continue to run SHELXL until the structure has either converged

        or until a maximum number of cycles has been reached

        Checks if structure has successfully refined

        Outputs statistics graphs (R1, Weight, Shift)

        Args:
            ins_file (str): full path to the new .ins file
            reference (str): full path to the reference .ins file
            refinements (int): number of refinements to check for shift convergence
            tolerance (float): target shift value
            max_cycles (int): maximum cycles SHELXL can run before stopping

        Returns:
            worked_flag (bool): whether or not the structure refined successfully
        """

        failure = False

        new_structure = pathlib.Path(ins_file)

        res_file = str(new_structure.stem) + ".res"

        os.chdir(new_structure.parent)

        self.import_refinement(new_structure.name, reference)

        df_weights = pd.DataFrame()
        df_shifts = pd.DataFrame()
        df_r_factor = pd.DataFrame()

        weight_list_1 = []
        weight_list_2 = []
        refinement_shifts = []
        r_factor_list = []

        convergence = False

        refine_count = 0

        while convergence == False and refine_count < max_cycles:
            refine_count += 1
            weight = ""
            new_weight = ""
            shelxl = subprocess.call(["shelxl", new_structure.stem])

            with open(res_file, "rt") as refinement:
                lines = refinement.readlines()

            worked_flag = False

            # The data has worked if the .res hasn't emptied itself AND if there is a CIF file present in the folder

            for line in lines:
                if "TITL" in line:
                    for item in os.listdir():
                        if new_structure.stem + ".cif" == item:
                            file_size = os.path.getsize(item)
                            if file_size > 0:
                                worked_flag = True

            if worked_flag == False:

                refine_count == max_cycles

            else:

                try:
                    shutil.copy(res_file, new_structure.name)
                except FileNotFoundError:
                    logging.info(__name__ + " : Refinement Failed")
                else:
                    with open(res_file, "rt") as refinement:
                        lines = refinement.readlines()
                        end_flag = False
                        second_weight_flag = False
                        for line in lines:
                            if end_flag == False and "WGHT" in line:
                                flag = 1
                                for item in line.split():
                                    if flag == 1:
                                        if item != "WGHT" and item != "":
                                            weight_list_1.append(float(item))
                                            flag += 1
                                    else:
                                        if item != "WGHT" and item != "":
                                            weight_list_2.append(float(item))
                                            second_weight_flag = True

                            elif end_flag == True and "WGHT" in line:
                                new_weight = line

                            elif "END" in line:
                                end_flag = True

                        if second_weight_flag == False:
                            weight_list_2.append(0)

                    # Copying the res to the ins file and changing the weight in the instruction section of the file

                    with open(new_structure.name, "rt") as initial:
                        lines = initial.readlines()

                        ACTA_flag = False
                        End_Flag = False

                        for line in lines:
                            if "ACTA" in line:
                                ACTA_flag = True

                        with open(new_structure.name, "w") as initial:
                            for line in lines:
                                if End_Flag == False:

                                    if "WGHT" in line and ACTA_flag == False:
                                        initial.write("ACTA \n")
                                        ACTA_flag = True
                                        initial.write(new_weight)
                                        old_weight = line

                                    elif "WGHT" in line and ACTA_flag == True:
                                        initial.write(new_weight)
                                        old_weight = line

                                    elif "END" in line:
                                        initial.write(line)
                                        End_Flag = True

                                    else:
                                        initial.write(line)
                                else:
                                    initial.write(line)

                        try:
                            logging.info(__name__ + " : " + str(old_weight))
                        except UnboundLocalError:
                            logging.info(__name__ + " : Structure died straight away")
                        else:
                            logging.info(__name__ + " : " + str(new_structure.name))
                            logging.info(__name__ + " : " + str(new_weight))

                            with open(
                                str(new_structure.stem) + ".lst", "rt"
                            ) as lst_file:
                                lst_lines = lst_file.readlines()
                            for line in lst_lines:
                                if "R1" and "Fo > 4sig(Fo)" in line:
                                    if line.split(" ")[3] != "":
                                        r_factor_list.append(float(line.split(" ")[3]))
                                    elif line.split(" ")[4] != "":
                                        r_factor_list.append(float(line.split(" ")[4]))
                                    else:
                                        print(
                                            "Error with reading R1 from shelxl. Plz investigate refinement.py line ~300"
                                        )
                                        exit()

                            if os.path.exists(str(new_structure.stem) + ".lst"):

                                with open(
                                    new_structure.stem + ".lst", "rt"
                                ) as refinement:
                                    lines = refinement.readlines()

                                (
                                    convergence,
                                    refinement_shifts,
                                    failure,
                                ) = self.convergence_check(
                                    lines,
                                    refinement_shifts,
                                    old_weight,
                                    new_weight,
                                    refinements,
                                    tolerance,
                                    failure,
                                )

                            else:
                                continue

                if failure == False:

                    # Because SHELXL will run multiple refinement cycles per execution, this will make sure that the r1, and weights have the same length as the shift stats for better graph communication

                    # Not having this doesn't break the code, it's just to make it nicer

                    if len(r_factor_list) == len(weight_list_1) == len(weight_list_2):

                        refinement_cycles = len(refinement_shifts) / refine_count
                        weight_list_1 = weight_list_1 + (
                            [weight_list_1[-1]] * (int(refinement_cycles) - 1)
                        )
                        weight_list_2 = weight_list_2 + (
                            [weight_list_2[-1]] * (int(refinement_cycles) - 1)
                        )
                        r_factor_list = r_factor_list + (
                            [r_factor_list[-1]] * (int(refinement_cycles) - 1)
                        )

        if failure == False and worked_flag == True:

            graphs = Grapher(self.test_mode)

            # Makes graphs of the changing weights, shifts and R-factors over time to graphically check for convergence

            if (
                len(refinement_shifts)
                == len(weight_list_1)
                == len(weight_list_2)
                == len(r_factor_list)
            ):

                self.figure_name = (
                    "Refinement_Statistics_" + str(new_structure.stem) + ".png"
                )
                x1 = list(range(1, len(weight_list_1) + 1))
                x2 = list(range(1, len(refinement_shifts) + 1))
                x3 = list(range(1, len(r_factor_list) + 1))
                x = [x1, x1, x2, x3]
                y = [weight_list_1, weight_list_2, refinement_shifts, r_factor_list]
                x_title = [
                    "Refinement Cycle",
                    "Refinement Cycle",
                    "Refinement Cycle",
                    "Refinement Cycle",
                ]
                y_title = ["Weight 1", "Weight 2", "Shift", "R1"]
                full_title = "Refinement Statistics - " + str(new_structure.name)
                mini_titles = [
                    "Weighting Convergence - " + str(new_structure.name),
                    "Weighting Convergence - " + str(new_structure.name),
                    "Shift Convergence - " + str(new_structure.name),
                    "R-Factor - " + str(new_structure.name),
                ]

                graphs.four_line_graph(
                    self.figure_name, x, y, x_title, y_title, full_title, mini_titles
                )

            else:
                self.figure_name = "No_Graph_Refinement_Failed"
                logging.info(
                    __name__
                    + " : Refinement fell over (ie **UNIT MISSING** or **REFINEMENT UNSTABLE** error)"
                )

        return worked_flag
