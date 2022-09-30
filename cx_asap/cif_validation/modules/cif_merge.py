#!/usr/bin/env python3

###################################################################################################
# ---------------------------------------CX-ASAP: refinement---------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from CifFile import ReadCif
import subprocess
import pathlib
import os
import logging
import shutil
from system_files.utils import Nice_YAML_Dumper, Config, Directory_Browse

# ----------Class Definition----------#


class Cif_Merge:
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

        # Set up yaml files and logger

        self.test_mode = test_mode

        config = Config(self.test_mode)

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

        # These are the entires it ignores from the instrument CIF (ie sometimes the synergy likes to have an extra part where everything is in P1

        self.ignored_entries = [
            "_space_group_it_number",
            "_space_group_crystal_system",
            "_space_group_name_h-m_alt",
            "_reflns_odcompleteness_completeness",
            "_reflns_odcompleteness_theta",
            "_reflns_odcompleteness_iscentric",
            "_chemical_oxdiff_formula",
            "_cell_length_a",
            "_cell_length_b",
            "_cell_length_c",
            "_cell_angle_alpha",
            "_cell_angle_beta",
            "_cell_angle_gamma",
            "_cell_volume",
            "_cell_formula_units_z",
            "_diffrn_oxdiff_digest_hkl",
            "_diffrn_oxdiff_digest_frames",
            "_diffrn_reflns_theta_full",
            "_diffrn_measured_fraction_theta_max",
            "_diffrn_measured_fraction_theta_full",
            "_diffrn_reflns_theta_min",
            "_diffrn_reflns_theta_max",
        ]

        self.flag1 = False
        self.flag2 = False

    def import_CIFs(self, instrument: str, new: str) -> None:

        """This function imports both a defined instrument cif and a structure cif

        It sets up the class parameters self.cif1 and self.cif2, which saves this data

        to the class. The PyCIFRW library is used to import the CIFs

        Args:
            instrument (str): Full path to the instrument cif file
            new (str): Full path to the file that the instrument details get merged into
        """

        instrument = pathlib.Path(instrument)
        new = pathlib.Path(new)

        self.flag1 = False
        self.flag2 = False

        os.chdir(instrument.parent)

        try:
            self.cif1 = ReadCif(instrument.name)
        except:
            self.flag1 = True
            logging.info(__name__ + " : Could not read instrument cif")

        os.chdir(new.parent)

        try:
            self.cif2 = ReadCif(new.name)
        except:
            self.flag2 = True
            logging.info(__name__ + " : Could not read structure cif")

        tree = Directory_Browse(os.getcwd(), self.test_mode)
        tree.check_file_contents(new)

    def synchrotron_cif_edit(self, instrument_cif: str) -> None:

        """This function edits the synchrotron autoprocess.cif so that it

        can be read in by PyCifRW

        This was required due to some weird formatting, and is largely unused for the

        main pipelines dedicated to post-reduction from homesource diffractometers

        Args:
            instrument (str): Full path to the instrument cif file
        """

        try:
            test = ReadCif(instrument_cif)
        except:

            with open(instrument_cif, "r") as f:
                lines = f.readlines()

            os.rename(instrument_cif, "old_autoprocess.txt")

            with open("temp.cif", "w") as f:
                if lines[0] != "data_autoprocess\n":
                    f.writelines("data_autoprocess\n")
                for line in lines:
                    if line.count("'") + line.count('"') > 2:
                        line_test = line.split(" ")
                        separator = " "
                        new_line = (
                            line_test[0] + "\n;" + separator.join(line_test[1:]) + ";\n"
                        )
                        f.writelines(new_line)
                    elif line.count("'") < 2 and "temperature" in line:
                        space_flag = False
                        letter_flag = False
                        changed_flag = False
                        for index, item in enumerate(line):
                            if (
                                item == " "
                                and space_flag == False
                                and letter_flag == False
                            ):
                                space_flag = True
                            elif (
                                item != " "
                                and space_flag == True
                                and letter_flag == False
                            ):
                                letter_flag = True
                            elif (
                                item == " "
                                and space_flag == True
                                and letter_flag == True
                            ):
                                new_line = line[0:index] + line[index + 1 :]
                                changed_flag = True
                        if changed_flag == False:
                            new_line = line
                        f.writelines(new_line)
                    elif line.count("'") == 0 and "temperature" not in line:
                        space_flag = False
                        letter_flag = False
                        changed_flag = False
                        for index, item in enumerate(line):
                            if (
                                item == " "
                                and space_flag == False
                                and letter_flag == False
                            ):
                                space_flag = True
                            elif (
                                item != " "
                                and space_flag == True
                                and letter_flag == False
                            ):
                                letter_flag = True
                                start_of_value = index
                                new_line = (
                                    line[0:index]
                                    + "'"
                                    + line[index:].strip("\n")
                                    + "'\n"
                                )
                                space_flag = False
                        if space_flag == False and letter_flag == True:
                            f.writelines(new_line)
                        else:
                            f.writelines(line)

                    else:
                        f.writelines(line)

            try:
                os.rename("temp.cif", "autoprocess.cif")
            except FileNotFoundError:
                pass

    def merge_CIFs(self) -> None:

        """This function merges the instrument cif and the

        refined cif using the PyCIF module.

        It adds a reference to CX-ASAP

        It also makes some small additional edits for checkCIF completeness :)

        """

        try:
            test = self.cif2.first_block()
        except:
            self.flag2 = True

        if self.flag2 == False and self.flag1 == True:

            self.data_block2 = self.cif2.first_block()

        elif self.flag2 == False and self.flag1 == False:

            self.data_block1 = self.cif1.first_block()
            self.data_block2 = self.cif2.first_block()

            data_to_change = self.data_block1.GetItemOrder()

            for item in data_to_change:
                if item not in self.ignored_entries and type(item) == str:
                    self.data_block2[item] = self.data_block1[item]

            if self.data_block2["_cell_measurement_reflns_used"] == "?":
                self.data_block2["_cell_measurement_reflns_used"] = self.data_block2[
                    "_diffrn_reflns_number"
                ]
                self.data_block2["_cell_measurement_theta_min"] = self.data_block2[
                    "_diffrn_reflns_theta_min"
                ]
                self.data_block2["_cell_measurement_theta_max"] = self.data_block2[
                    "_diffrn_reflns_theta_max"
                ]
            self.data_block2["_cell_measurement_temperature"] = self.data_block2[
                "_diffrn_ambient_temperature"
            ]

            self.data_block2["_refine_special_details"] = "CX-ASAP (Thompson, 2021)"

    def user_edits(
        self,
        solution: str,
        formula: str,
        colour: str,
        habit: str,
        cryst_max: str,
        cryst_mid: str,
        cryst_min: str,
    ) -> None:

        """This function puts all the additional user-defined parameters in

        The user defines these parameters in the yaml file

        Args:
            solution(str): the structure solution program used
            formula(str): the chemical formula
            habit(str): the crystal habit
            cryst_max(str): the largest dimension of the crystal
            cryst_mid(str): the middle dimension of the crystal
            cryst_min(str): the smallest dimension of the crystal
        """

        if self.flag2 == False and self.flag1 == False:

            self.data_block2["_computing_structure_solution"] = solution
            self.data_block2["_chemical_formula_moiety"] = formula
            self.data_block2["_exptl_crystal_colour"] = colour
            self.data_block2["_exptl_crystal_description"] = habit
            self.data_block2["_exptl_crystal_size_max"] = cryst_max
            self.data_block2["_exptl_crystal_size_mid"] = cryst_mid
            self.data_block2["_exptl_crystal_size_min"] = cryst_min

    def single_write_out(self, file_name: str) -> None:

        """This function will write out the edited cif stored in the class

        Args:
            file_name(str): the name of the output file (or full path to)
        """

        if self.flag2 == False:

            with open(file_name, "w") as f:

                f.write(self.cif2.WriteOut())

    def write_out(
        self, location: str, cif_name: str, validation_name: str, single_cif_name: str
    ) -> None:

        """This function writes out the edited cif stored in the class

        It differes from the 'single_write_out' function by

        also running checkCIF, and has the ability to append to

        a pre-existing CIF. Ie, this is where the CIFs get MERGED.

        It also prints out the checkCIF report.

        If run as a stand-alone module, files should be output to the same

        location as the original cif

        If run as a part of the pipeline, the files should be output where

        the results files go

        Args:
            location(str): where the files are written out to
            cif_name(str): name of the output cif
            validation_name(str): name of the output check_CIF
            single_cif_name(str): the name of the original file
        """

        if self.flag2 == False:

            current_folder = os.getcwd()
            os.chdir(location)

            bad_flag = False

            with open(single_cif_name, "rt") as f:
                lines = f.readlines()

            for line in lines:
                if "data_xcalibur" in line:
                    bad_flag = True

            if bad_flag == False:

                logging.info(
                    __name__
                    + " : CIF being added to combined is: "
                    + str(single_cif_name)
                )

                with open(cif_name, "a") as f:
                    f.write(self.cif2.WriteOut())

                try:
                    test = self.validation
                except:
                    pass
                else:

                    with open(validation_name, "a") as f:
                        for line in self.validation:
                            f.write(line)

            os.chdir(current_folder)

    def validate_CIFs(self, file_name: str) -> None:

        """Run checkCIF on a specified file

        then saves the data to the class as self.validation

        Args:
            location(str): the name of the cif file for platon checkCIF
        """

        # This function runs each cif through platon for validation, and adds the report to a combined check cif file
        # This is variable that isn't used yet:
        # windows_wait = 15

        if self.flag2 == False:
            checkCIF = subprocess.Popen(
                ["platon", "-u", file_name], stdin=subprocess.PIPE, encoding="utf8"
            )
            try:
                checkCIF.wait(45)
            except subprocess.TimeoutExpired:
                checkCIF.kill()

            try:
                shutil.move(pathlib.Path(file_name).stem + ".chk", "check_CIF.chk")
            except FileNotFoundError:
                self.validation = ""
            else:
                for item in os.listdir():
                    if item.endswith(".chk"):
                        with open(item, "rt") as f:
                            self.validation = f.readlines()
