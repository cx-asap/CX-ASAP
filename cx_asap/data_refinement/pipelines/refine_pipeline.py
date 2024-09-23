#!/usr/bin/env python3

###################################################################################################
# -------------------------------------CX-ASAP: refine_pipeline------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config, Directory_Browse
from data_refinement.modules.refinement import Structure_Refinement
import shutil
import os
import logging

# ----------Class Definition----------#


class Refinement_Pipeline:
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

    def multiple_refinement(
        self,
        location: str,
        reference: str,
        graph_output_location: str,
        refinements_to_check: int,
        tolerance: float,
        max_cycles: int,
    ) -> None:
        """Runs SHELXL on a series of structures based on one reference

        All of the new structures should be folders within a common folder

        A report detailing which structures did and did not refine successfully

        is output as both a file and to the terminal

        Args:
            location (str): full path to the folder containing folders of .ins files
            reference (str): full path to the reference .ins/.res file
            graph_output_location (str): full path to the location of output files
            refinements_to_check (int): number of refinements to check for shift convergence
            tolerance (float): target shift value
            max_cycles (int): maximum cycles SHELXL can run before stopping
        """

        successful_structures = []

        failed_structures = []

        self.tree = Directory_Browse(location, self.test_mode)
        self.shelxl = Structure_Refinement(self.test_mode)
        for item in self.tree.directories:

            shelxl_run_flag = False

            self.tree.enter_directory(item, ".ins")
            if self.tree.item_file != "":
                outcome = self.shelxl.run_shelxl(
                    self.tree.item_file,
                    reference,
                    refinements_to_check,
                    tolerance,
                    max_cycles,
                )
                self.tree.check_file_contents()

                shelxl_run_flag = True

                # Copies the statistics graphs to the outer folder for easier comparison

                # Ie it saves the user having to dig through each individual folder to look at them

                try:
                    shutil.copy(self.shelxl.figure_name, graph_output_location)
                except:
                    logging.info(__name__ + " : Refinement failed so no graph :( ")

            else:
                if shelxl_run_flag == True:
                    logging.info(__name__ + " : Failed to get .ins file")
                    outcome = False
                else:
                    logging.info(__name__ + " : No .ins file in folder " + str(item))
                    outcome = False

            if outcome == True and shelxl_run_flag == True:
                successful_structures.append(self.tree.item_file)
            elif outcome == False and shelxl_run_flag == True:
                failed_structures.append(self.tree.item_file)
                for i in os.listdir(os.getcwd()):
                    if i == self.tree.item_name + ".cif":
                        os.rename(i, i + "_old")

            self.tree.exit_directory()

        a = "------------------------------"
        b = "------Refinement Summary------"
        c = "------------------------------"
        d = (
            "Successful refinements: "
            + str(len(successful_structures))
            + " of "
            + str(len(successful_structures) + len(failed_structures))
        )
        e = "Possible reasons for refinement failure = change of symmetry or cell setting"
        i = "If some structures failed to refine, you should check the data reduction"
        g = "Successful refinements:"
        h = "Failed refinements:"

        os.chdir(graph_output_location)

        with open("refinement_summary.txt", "w") as f:
            f.write(str(a) + "\n")
            f.write(str(b) + "\n")
            f.write(str(c) + "\n")
            f.write(str(d) + "\n")
            f.write(str(e) + "\n")
            f.write(str(i) + "\n")
            f.write(str(g) + "\n")

            for item in successful_structures:
                try:
                    f.write(item.name + "\n")
                except:
                    f.write("ERROR" + "\n")
                    logging.info(__name__ + " : Folder without .ins analysed")
            f.write(str(h) + "\n")

            for item in failed_structures:
                try:
                    f.write(item.name + "\n")
                except:
                    f.write("ERROR" + "\n")
                    logging.info(__name__ + " : Folder without .ins analysed")
            f.write("Reference location:" + str(reference) + "\n" + "\N{rocket}")

        print(a)
        print(b)
        print(c)
        print(d)
        print(e)
        print(i)
        print(g)
        for item in successful_structures:
            try:
                print(item.name)
            except:
                print("ERROR")
        print(h)

        if len(failed_structures) != 0:
            for item in failed_structures:
                try:
                    print(item.name)
                except:
                    print("ERROR")
        else:
            print("None")
        print(a)
