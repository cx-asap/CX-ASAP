#!/usr/bin/env python3

###################################################################################################
# ----------------------------------------CX-ASAP: utilities---------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

import pathlib
import yaml
import os
import shutil
import logging
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math
import re
import fileinput
from typing import Tuple

# ----------Class Definition----------#


class Configure_Flexible:
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

        self.XDS = XDS_File_Edit()
        self.cell = Cell_Import()
        self.home_path = self.sys["home_path"]
        self.bkg_files = [
            "GAIN.cbf",
            "BKGINIT.cbf",
            "BLANK.cbf",
            "X-CORRECTIONS.cbf",
            "Y-CORRECTIONS.cbf",
        ]

        # This is a variable kept for naming the files to keep track of experiments - ie if testing various angles, each test needs a unique identifier, and this is it! So it needs to be reset every time the pipeline is run!

        self.sys["process_counter"] = 0

    def file_movement(
        self,
        background_files_reference_path: str,
        instrument_parameters_path: str,
        instrument_cif: str,
    ) -> None:

        """Sets up the directory tree for a specialised flexible mapping experiment

        Basically its just moving around files and making things pretty

        Args:
            background_files_reference_path (str): full path to reference XDS background files
                                                    BKGINIT.cbf
                                                    BLANK.cbf
                                                    GAIN.cbf
                                                    X-CORRECTIONS.cbf
                                                    Y-CORRECTIONS.cbf
            instrument_parameters_path (str): full path to the reference GXPARM.XDS file
            instrument_cif (str): full path to the instrument CIF

        """

        os.chdir(self.sys["analysis_path"])

        # Copies instrument cif into ref folder if not done previously

        self.instrument_cif_path = (
            pathlib.Path(self.sys["ref_path"]) / pathlib.Path(instrument_cif).name
        )

        if not os.path.exists(self.instrument_cif_path):
            try:
                shutil.copy(pathlib.Path(instrument_cif), self.sys["ref_path"])
            except FileNotFoundError as error:
                logging.critical(
                    __name__
                    + " : "
                    + f"Failed to find instrument cif with error {error}"
                )
                print("Error! see logs")
                exit()
            else:
                self.instrument_cif_path = (
                    pathlib.Path(self.sys["ref_path"])
                    / pathlib.Path(instrument_cif).name
                )

        # Moves background files into reference folder AND copies into each run folder if not done previously

        for folder in os.listdir():
            for item in self.bkg_files:
                if not os.path.exists(pathlib.Path(self.sys["ref_path"]) / item):
                    shutil.copy(
                        pathlib.Path(background_files_reference_path) / item,
                        pathlib.Path(self.sys["ref_path"]) / item,
                    )
                shutil.copy(pathlib.Path(self.sys["ref_path"]) / item, folder)
            shutil.copy(self.instrument_cif_path, folder)

        os.chdir(self.home_path)

        # Copies GXPARM into ref folder if not done previously

        if not os.path.exists(pathlib.Path(self.sys["ref_path"], "GXPARM.XDS")):
            try:
                shutil.copy(instrument_parameters_path, self.sys["ref_path"])
            except FileNotFoundError as error:
                logging.critical(
                    __name__
                    + " : "
                    + f"Failed to find instrument parameters with error {error}"
                )
                print("Error! See Logs")
                exit()

    def flexible_setup(
        self,
        max_processors: int,
        neggia_library: str,
        space_group_number: int,
        MPLA_atoms: str,
        total_angle: int,
    ) -> None:

        """Primarily sets up the XDS.INP file for a flexible crystal experiment

        Also adds in the MPLA command to the reference structure

        Args:
            max_processors (int): maximum number of processors for XDS to use
            neggia_library (str): full path to dectris-neggia.so
            space_group_number (int): space group number to go into XDS.inp
            MPLA_atoms (str): atoms for MPLA command to be written into reference
                                .ins/.res file
            total_angle (int): total wedge angle measured in the experiment

        """

        # Finds the starting and total angle and saves them (allows for later calculations of wedge angles and also allows for them to be put back into the file at the end

        self.sys["start_angle"] = self.XDS.start_angle(self.sys["XDS_inp_organised"])

        self.sys["total_angle"] = total_angle

        self.sys["total_frames"] = self.XDS.get_value(
            self.sys["XDS_inp_organised"], "DATA_RANGE"
        )[1:]

        with open(self.sys_path, "w") as f:
            yaml.dump(
                self.sys,
                f,
                default_flow_style=False,
                Dumper=Nice_YAML_Dumper,
                sort_keys=False,
            )

        self.cfg, self.sys = self.config.yaml_reload()

        self.XDS.new_line_rewrite(self.sys["XDS_inp_organised"])

        with open(
            pathlib.Path(self.sys["ref_path"]) / "GXPARM.XDS", "rt"
        ) as instrument:

            # Locates the instrument parameters needed to copy into the XDS.INP File

            for index, line in enumerate(instrument):
                if index == 8:
                    data = line.strip("\n").split(" ")
        data2 = []
        for element in data:
            if "." in element:
                data2.append(element)

        orgx = data2[0]
        orgy = data2[1]
        dist = data2[2]

        self.XDS.change(self.sys["XDS_inp_organised"], "ORGX", orgx)
        self.XDS.change(self.sys["XDS_inp_organised"], "ORGY", orgy)
        self.XDS.change(self.sys["XDS_inp_organised"], "DETECTOR_DISTANCE", dist)

        # Sets the initial XDS.INP Parameters

        self.XDS.change(
            self.sys["XDS_inp_organised"],
            "MAXIMUM_NUMBER_OF_PROCESSORS",
            max_processors,
        )
        self.XDS.change(self.sys["XDS_inp_organised"], "LIB", neggia_library)
        self.XDS.change(
            self.sys["XDS_inp_organised"],
            "JOB",
            "COLSPOT IDXREF DEFPIX INTEGRATE CORRECT",
        )
        self.XDS.change(
            self.sys["XDS_inp_organised"], "NAME_TEMPLATE_OF_DATA_FRAMES", "img"
        )
        self.XDS.change(self.sys["XDS_inp_organised"], "SPOT_MAXIMUM-CENTROID", 2)
        self.XDS.change(
            self.sys["XDS_inp_organised"], "MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT", 2
        )
        self.XDS.change(self.sys["XDS_inp_organised"], "STRONG_PIXEL", 2)

        with open(self.sys["XDS_inp_organised"], "rt") as in_file:
            flag1 = 0
            flag2 = 0
            flag3 = 0
            flag4 = 0
            flag5 = 0
            flag6 = 0
            for line in in_file:
                if "MINIMUM_FRACTION_OF_INDEXED_SPOTS" in line:
                    self.XDS.change(
                        self.sys["XDS_inp_organised"],
                        "MINIMUM_FRACTION_OF_INDEXED_SPOTS",
                        0.2,
                    )
                    flag1 += 1
                elif "SEPMIN" in line:
                    self.XDS.change(self.sys["XDS_inp_organised"], "SEPMIN", 7)
                    flag2 += 1
                elif "REFINE(INTEGRATE)" in line:
                    self.XDS.change(
                        self.sys["XDS_inp_organised"],
                        "REFINE(INTEGRATE)",
                        "ORIENTATION CELL AXIS !BEAM POSITION",
                    )
                    flag3 += 1
                elif "REFINE(IDXREF)" in line:
                    self.XDS.change(
                        self.sys["XDS_inp_organised"],
                        "REFINE(IDXREF)",
                        "ORIENTATION CELL AXIS !BEAM POSITION",
                    )
                    flag4 += 1
                elif "REFINE(CORRECT)" in line:
                    self.XDS.change(
                        self.sys["XDS_inp_organised"],
                        "REFINE(CORRECT)",
                        "ORIENTATION CELL AXIS !BEAM POSITION",
                    )
                    flag5 += 1
                elif "CLUSTER_RADIUS" in line:
                    self.XDS.change(
                        self.sys["XDS_inp_organised"], "CLUSTER_RADIUS", 3.5
                    )
                    flag6 += 1

        with open(self.sys["XDS_inp_organised"], "a") as in_file:
            if flag1 == 0:
                in_file.write(" MINIMUM_FRACTION_OF_INDEXED_SPOTS= 0.2\n")
            if flag2 == 0:
                in_file.write(" SEPMIN= 7\n")
            if flag3 == 0:
                in_file.write(
                    " REFINE(INTEGRATE)= ORIENTATION CELL AXIS !BEAM POSITION\n"
                )
            if flag4 == 0:
                in_file.write(" REFINE(IDXREF)= ORIENTATION CELL AXIS !BEAM POSITION\n")
            if flag5 == 0:
                in_file.write(
                    " REFINE(CORRECT)= ORIENTATION CELL AXIS !BEAM POSITION\n"
                )
            if flag6 == 0:
                in_file.write(" CLUSTER_RADIUS= 3.5\n")

        self.cfg["wavelength"] = self.XDS.get_value(
            self.sys["XDS_inp_organised"], "X-RAY_WAVELENGTH"
        )

        with open(self.sys_path, "w") as f:
            yaml.dump(
                self.sys,
                f,
                default_flow_style=False,
                Dumper=Nice_YAML_Dumper,
                sort_keys=False,
            )

        self.cell.cell_import(self.sys["ref_path_organised"])

        self.cfg, self.sys = self.config.yaml_reload()

        # IF CANNOT FIND START_ANGLE AND TOTAL_ANGLE, COMMENT OUT CELL_IMPORT ABOVE, RUN CODE, THEN MAGICALLY FIXED??? should only happen after a fresh install of cxasap.....

        # Puts the reference unit cell into the XDS.INP file

        self.XDS.change(
            self.sys["XDS_inp_organised"],
            "UNIT_CELL_CONSTANTS",
            str(self.sys["ref_a"])
            + " "
            + str(self.sys["ref_b"])
            + " "
            + str(self.sys["ref_c"])
            + " "
            + str(self.sys["ref_alpha"])
            + " "
            + str(self.sys["ref_beta"])
            + " "
            + str(self.sys["ref_gamma"]),
        )

        # XDS Accepts the space group NUMBER not the name - currently this has to be entered manually as the reference INS file expresses the space group as the NAME not the number ( future need to implement a dictionary to translate between NAME and NUMBER

        self.XDS.change(
            self.sys["XDS_inp_organised"], "SPACE_GROUP_NUMBER", space_group_number
        )

        self.cell.ref_edit(self.sys["ref_path_organised"], MPLA_atoms)


# ----------Class Definition----------#


class Directory_Browse:
    def __init__(self, location, error_mode: bool = False) -> None:

        """Initialises the class

        Sets up the yaml parameters input by the user

        Also defines the location for the system yaml file

        which stores a yaml of code-only parameters accessible throughout

        the software package

        Args:
            error_mode (bool): whether or not error mode is run (relevant to
                                yaml config)

        """

        # Set up yaml files and logger

        # for some reason changed test_mode to error_mode??? Same thing?

        config = Config(error_mode)
        self.error_mode = error_mode

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

        # Finds all of the directories in the desired location

        self.home_directory = location
        current_directory = os.getcwd()
        os.chdir(self.home_directory)
        self.directories = []
        self.sort = File_Sorter()
        for item in self.sort.sorted_properly(os.listdir(self.home_directory)):
            if os.path.isdir(item):
                self.directories.append(pathlib.Path(item).absolute())
        os.chdir(current_directory)

    def enter_directory(
        self,
        folder: str,
        file_suffix: str,
        ignored_files: list = False,
        ignored_folders: list = [],
        ignore_check: bool = False,
    ) -> None:

        """Enter a directory and search for a specified file

        The idea here is that upon creating the class, it will have

        the property object.directories

        The pipeline will then loop through object.directories

        and within this loop use functions such as 'enter_directory' etc

        folder is used as a variable, so (for example) if there is not

        actually a series of folders to scan,

        then you can still use this function to enter just one

        directory of your chosing (without using self.directories)

        Imports a unit cell from a reference .ins/.res file

        saves it in the sys.yaml file

        ALSO runs a check to see if there are any duplicate files in the folder

        that would break the code

        Args:
            folder (str): folder name for entering
            file_suffix (str): suffix to check for duplicates
            ignored_files (list): any files to ignore
            ignored_folders (list): any folders to ignore
            ignore_check (bool): if true, ignores the file check,
                                if false, checks for  duplicates
        """

        if folder not in ignored_folders:
            os.chdir(folder)
            logging.info(__name__ + " : Performing tasks in folder: " + folder.name)
            self.item_file = ""
            self.item_name = ""

            if ignore_check == False:

                check = File_Check(self.error_mode)

                file_list = check.duplicate_check(file_suffix, ignored_files)

                if len(file_list) >= 2:
                    logging.info(
                        __name__
                        + " : Too many files present with suffix: "
                        + file_suffix
                    )
                    logging.info(
                        __name__
                        + " : Please make sure only one of this file type is present in the folder"
                    )

                    print(
                        "There are too many files present with suffix: " + file_suffix
                    )
                    print(file_list)
                    print(
                        "Please make sure only one of this file type is present in the folder."
                    )

                    exit()

            for item in os.listdir():
                if item.endswith(file_suffix) and item != ignored_files:
                    self.item_file = pathlib.Path(pathlib.Path(os.getcwd()) / item)
                    self.item_name = self.item_file.stem
                    logging.info(
                        __name__ + " : File name for analysis: " + str(self.item_file)
                    )

    def enter_directory_multiple(
        self, folder: str, file_suffix: str, ignored_files: list = False
    ) -> None:

        """Similar to above function, but it is for when there are multiple

        of the desired files in the folder

        Ie a bunch of CIFs that all get analysed in a single folder

        Args:
            folder (str): folder name for entering
            file_suffix (str): suffix to check for duplicates
            ignored_files (list): any files to ignore
        """

        self.item_files = []
        self.item_names = []

        os.chdir(folder)
        logging.info(__name__ + " : Performing tasks in folder: " + folder.name)
        for item in self.sort.sorted_properly(os.listdir()):
            if item.endswith(file_suffix) and item != ignored_files:
                self.item_files.append(pathlib.Path(item))
                self.item_names.append(pathlib.Path(item).stem)
                logging.info(__name__ + " : File name for analysis: " + str(item))

    def exit_directory(self) -> None:

        """Exits back to the home directory"""

        os.chdir(self.home_directory)

    def check_file_contents(self, file_check: str = False) -> None:

        """Checks the contents of a specified file to see if it is empty

        Some crystallographic programs will write out an empty file if

        the process has failed

        This checks for that and deletes empty files

        Args:
            file_check (str): the file to be checked
        """

        if file_check == False:
            file_check = self.item_file

        # Ability to check file contents to see if anything in them

        # Will remove any empty files

        try:
            file_size = os.stat(file_check)
        except FileNotFoundError:
            logging.info(__name__ + " : file not found: " + str(file_check))
        else:
            try:
                if file_size.st_size < 1:
                    logging.info(__name__ + " : Empty file: " + str(file_check))
                    os.remove(file_check)
            except PermissionError:
                logging.info(
                    __name__ + " : windows permission error: " + str(file_check)
                )
                print(
                    "WINDOWS CREATED AN EMPTY STUPID FILE CALLED "
                    + str(file_check)
                    + " in "
                    + str(os.getcwd())
                )
            pass


# ----------Class Definition----------#


class File_Check:
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

    def duplicate_check(self, file_suffix: str, ignored_files: list = False) -> list:

        """Checks if there are two files with the same suffix in a folder

        Ie two .ins files means in a folder would break the code

        Args:
            file_suffix (str): file ending to check for duplicates
            ignored_files (list): list of any ignored files

        Returns:
            files_with_suffix (list): list of files with a common ending
        """

        files_with_suffix = []

        for item in os.listdir():
            if item.endswith(file_suffix) and item != ignored_files:
                files_with_suffix.append(item)

        return files_with_suffix


# ----------Class Definition----------#


class Autoprocess_Setup:
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

    def config(self, location: str, experiment_name: str, experiment_type: str) -> None:

        """Sets up the directory tree structure and updates the paths in

        the sys.yaml file

        Args:
            location (str): full path to the location of folders of autoprocessed data
            experiment_name (str): common name between data sets
            experiment_type (str): type of experiment (only for naming new folder)
        """

        os.chdir(location)
        self.experiment_name = experiment_name
        self.experiment_type = experiment_type

        self.home_path = pathlib.Path(location) / (
            experiment_name + "_" + experiment_type
        )
        self.tree_structure = ["analysis", "ref", "results", "failed_autoprocessing"]

        self.analysis_path = pathlib.Path(self.home_path) / self.tree_structure[0]
        self.ref_path = pathlib.Path(self.home_path) / self.tree_structure[1]
        self.results_path = pathlib.Path(self.home_path) / self.tree_structure[2]
        self.failed_path = pathlib.Path(self.home_path) / self.tree_structure[3]

        self.sys["home_path"] = str(self.home_path.absolute())
        self.sys["analysis_path"] = str(self.analysis_path.absolute())
        self.sys["ref_path"] = str(self.ref_path.absolute())
        self.sys["results_path"] = str(self.results_path.absolute())
        self.sys["failed_path"] = str(self.failed_path.absolute())

        with open(self.sys_path, "w") as f:
            yaml.dump(
                self.sys,
                f,
                default_flow_style=False,
                Dumper=Nice_YAML_Dumper,
                sort_keys=False,
            )

    def File_Rename_AS(self) -> None:

        """This next part adds numbers so folders are correctly ordered

        Also allows for the fact that if the screenings auto-process

        they will be present too

        """

        os.chdir(self.analysis_path)

        list_split = []
        data_dict = {}
        data = []

        for item in os.listdir():
            if os.path.isdir(item):
                list_split.append(item.split("_"))
                data.append(item)
        index = 0

        for item in list_split:
            data_dict[item[-3]] = []
        for item in list_split:
            data_dict[item[-3]].append(index)
            index += 1

        keys = sorted(data_dict.keys())

        counter = 1

        for key in keys:
            for value in data_dict[key]:
                os.rename(data[value], str(counter) + "_" + data[value])
                counter += 1

        os.chdir(self.home_path)

    def Organise_Directory_Tree(
        self, reference_location: str, Synchrotron: bool = False
    ) -> None:

        """Organises the directory tree, moves a lot of files around, and

        makes sure that any that failed autoprocessing are separted from

        those that have worked

        Args:
            reference_location (str): full path to the reference .ins/.res file
            Synchrotron (bool): if from the Australian Synchrotron, performs the
                                above file renaming function
        """

        original_reference_file = pathlib.Path(reference_location)

        # Folders will be made and files organised if not already done in a previous run

        if not os.path.exists(self.home_path):
            os.mkdir(self.home_path.name)
            os.chdir(self.home_path)
            for item in self.tree_structure:
                os.mkdir(item)

            # Move reference into ref folder

            try:
                shutil.copy(reference_location, self.ref_path)
            except FileNotFoundError as error:
                logging.critical(
                    __name__
                    + " : "
                    + f"Failed to find reference file with error {error}"
                )
                print("Error - Please check logs")
                exit()

            os.chdir("..")

            # Move experiments into analysis folder based on their name

            for item in os.listdir():
                if self.experiment_name in item:
                    if self.experiment_type not in item:

                        # makes sure that the folder that everything is being moved into is not considered

                        shutil.move(item, pathlib.Path(self.analysis_path) / item)

            os.chdir(self.analysis_path)

            # Finds out which runs failed to autoprocess and moves them to the failed folder

            for run in os.listdir():
                os.chdir(run)
                print(run)
                if "autoprocess.cif" not in os.listdir():
                    shutil.move(
                        pathlib.Path(self.analysis_path) / run,
                        pathlib.Path(self.failed_path) / run,
                    )
                    print(run + "failed")
                os.chdir(self.analysis_path)

            if Synchrotron == True:
                self.File_Rename_AS()

        os.chdir(self.results_path)

        # Makes a new folder in the results for each time the code is run - keeps tests separated

        tmp = os.listdir()
        tmp2 = []
        for item in tmp:
            if os.path.isdir(item):
                tmp2.append(int(item))
        tmp2.sort()

        try:
            self.new_folder = int(tmp2[-1]) + 1
        except IndexError:
            self.new_folder = 1

        os.mkdir(pathlib.Path(self.results_path) / str(self.new_folder))

        self.sys["current_results_path"] = str(
            pathlib.Path(
                pathlib.Path(self.results_path) / str(self.new_folder)
            ).absolute()
        )
        self.sys["ref_path_organised"] = str(
            pathlib.Path(
                pathlib.Path(self.ref_path) / original_reference_file.name
            ).absolute()
        )

        with open(self.sys_path, "w") as f:
            yaml.dump(
                self.sys,
                f,
                default_flow_style=False,
                Dumper=Nice_YAML_Dumper,
                sort_keys=False,
            )

        # Gets rid of all previously made .ins files from the analysis folders to make sure that the code runs ok later

        os.chdir(self.analysis_path)

        for run in os.listdir():
            if os.path.isdir(run):
                os.chdir(run)
                for item in os.listdir():
                    if ".ins" in item:
                        os.remove(item)
                os.chdir("..")

        os.chdir(self.home_path)


# ----------Class Definition----------#


class Reprocess_Setup:
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

    def config(self, location: str, experiment_name: str, experiment_type: str) -> None:

        """Sets up the directory tree structure and updates the paths in

        the sys.yaml file

        Args:
            location (str): full path to the location of folders of autoprocessed data
            experiment_name (str): common name between data sets
            experiment_type (str): type of experiment (only for naming new folder)
        """

        os.chdir(location)
        self.experiment_name = experiment_name
        self.experiment_type = experiment_type

        self.home_path = pathlib.Path(location) / (
            experiment_name + "_" + experiment_type
        )
        self.tree_structure = ["analysis", "ref", "results", "frames"]

        self.analysis_path = pathlib.Path(self.home_path) / self.tree_structure[0]
        self.ref_path = pathlib.Path(self.home_path) / self.tree_structure[1]
        self.results_path = pathlib.Path(self.home_path) / self.tree_structure[2]
        self.frames_path = pathlib.Path(self.home_path) / self.tree_structure[3]

        self.sys["home_path"] = str(self.home_path.absolute())
        self.sys["analysis_path"] = str(self.analysis_path.absolute())
        self.sys["ref_path"] = str(self.ref_path.absolute())
        self.sys["results_path"] = str(self.results_path.absolute())
        self.sys["frames_path"] = str(self.frames_path.absolute())

        with open(self.sys_path, "w") as f:
            yaml.dump(
                self.sys,
                f,
                default_flow_style=False,
                Dumper=Nice_YAML_Dumper,
                sort_keys=False,
            )

    def Organise_Directory_Tree(
        self, reference_location: str, xds_inp_location: str
    ) -> None:

        """Moves a lot of files around and sets everything up for automatic analysis

        Args:
            reference_location (str): full path to the reference .ins/.res file
            xds_inp_location (str): full path to the reference XDS.INP file
        """

        original_reference_file = pathlib.Path(reference_location)
        original_xds_file = pathlib.Path(xds_inp_location)

        # Folders will be made and files organised if not already done in a previous run

        if not os.path.exists(self.home_path):
            os.mkdir(self.home_path.name)
            os.chdir(self.home_path)
            for item in self.tree_structure:
                os.mkdir(item)

            # Move references into ref folder

            try:
                shutil.copy(reference_location, self.ref_path)
            except FileNotFoundError as error:
                logging.critical(
                    __name__
                    + " : "
                    + f"Failed to find reference file with error {error}"
                )
                print("Error - Please check logs")

                exit()

            try:
                shutil.copy(xds_inp_location, self.ref_path)
            except FileNotFoundError as error:
                # self.logger.critical(f'Failed to find reference file with error {error}')
                logging.critical(
                    __name__
                    + " : "
                    + f"Failed to find reference file with error {error}"
                )
                print("Error - Please check logs")
                exit()

            os.chdir("..")

            # For old detector images:

            folders_made = []

            for item in os.listdir():
                if self.experiment_name in item:
                    if "analysis" not in item:
                        shutil.move(item, pathlib.Path(self.frames_path) / item)

                        if item.endswith("master.h5"):
                            b = item.replace("_master.h5", "")
                            os.mkdir(pathlib.Path(self.analysis_path) / b)

                        # For old detector images:

                        if item.endswith(".img"):
                            b = "_".join(item.split("_")[:-1])
                            if b not in folders_made:
                                os.mkdir(pathlib.Path(self.analysis_path) / b)
                                folders_made.append(b)

            os.chdir(self.analysis_path)

            for folder in os.listdir():
                os.symlink(self.frames_path, os.path.join(folder, "img"))

        os.chdir(self.results_path)

        # Makes a new folder in the results for each time the code is run - keeps tests separated

        tmp = os.listdir()
        tmp2 = []
        for item in tmp:
            if os.path.isdir(item):
                tmp2.append(int(item))
        tmp2.sort()

        try:
            self.new_folder = int(tmp2[-1]) + 1
        except IndexError:
            self.new_folder = 1

        os.mkdir(pathlib.Path(self.results_path) / str(self.new_folder))

        self.sys["current_results_path"] = str(
            pathlib.Path(
                pathlib.Path(self.results_path) / str(self.new_folder)
            ).absolute()
        )
        self.sys["ref_path_organised"] = str(
            pathlib.Path(
                pathlib.Path(self.ref_path) / original_reference_file.name
            ).absolute()
        )
        self.sys["XDS_inp_organised"] = str(
            pathlib.Path(
                pathlib.Path(self.ref_path) / original_xds_file.name
            ).absolute()
        )

        with open(self.sys_path, "w") as f:
            yaml.dump(
                self.sys,
                f,
                default_flow_style=False,
                Dumper=Nice_YAML_Dumper,
                sort_keys=False,
            )

        # Gets rid of all previously made .ins files from the analysis folders to make sure that the code runs ok later

        os.chdir(self.home_path)


# ----------Class Definition----------#


class Grapher:
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

        config = Config(test_mode)

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    def single_scatter_graph(
        self,
        x: list,
        y: list,
        x_title: str,
        y_title: str,
        graph_title: str,
        figure_name: str,
        y_series_title: list = None,
        colour: list = None,
        marker: list = None,
        linewidth: list = None,
        s: list = None,
    ) -> None:

        """This function is used to print out a single scatter graph

        Options include having multiple series on the same graph

        Args:
            x (list): x-data
            y (list): y-data
            x_title (str): label for x-axis
            y_title (str): label for y-axis
            graph_title (str): title for graph
            figure_name (str): name for the output file
            y_series_title (list): names of y-series if multiple plots to go in legend
            colour (list): colours for multiple series
            marker (list): markers for multiple series
            linewidth (list): line widths for multiple series
            s (list): marker sizes for multiple series
        """

        for index, item in enumerate(y):

            if type(item) != float:

                if len(x) != len(item):
                    logging.info(
                        __name__
                        + " : Possible error with plotting structural changes. Check the structures in the output CIF for unreasonable structures."
                    )

                    # This will manually make x the correct length by repeating value

                    # Often just an artefact of weird things going on

                    # Also stops the pipeline-variable-analysis throwing error in test

                    # Was because when this pipeline is run individually on test data,

                    # all the temperatures hadn't been edited yet

                    to_repeat = x[0]

                    x = [to_repeat] * len(item)

                if colour == None and y_series_title != None:
                    plt.scatter(x, item, label=y_series_title[index])
                elif y_series_title == None:
                    plt.scatter(x, y)
                else:
                    plt.scatter(
                        x,
                        item,
                        c=colour[index],
                        marker=marker[index],
                        linewidth=linewidth[index],
                        s=s[index],
                        label=y_series_title[index],
                    )

            else:
                plt.scatter(x, y)

        plt.xlabel(x_title, fontsize=12)
        plt.ylabel(y_title, fontsize=12)

        plt.title(graph_title)
        if y_series_title != None:
            plt.legend(fontsize=12)

        figure = plt.gcf()
        figure.set_size_inches(16, 12)

        plt.savefig(figure_name, bbox_inches="tight", dpi=100)
        plt.clf()
        plt.close("all")

    def multi_scatter_graph(
        self,
        x: list,
        y: list,
        title: list,
        rows: int,
        columns: int,
        subplots: list,
        y_series_title: str,
        x_title: str,
        y_title: str,
        figure_name: str,
    ) -> None:

        """This function is for doing multiple separate scatter graphs on the

        same figure (ie 8 different scatter graphs)

        Args:
            x (list): x-data
            y (list): y-data
            title (str): title for graph
            rows (int): number of rows of graphs in figure
            columns (int): number of columns of graphs in figure
            subplots (list): list of graphs ie for 4 graphs it should be [1,2,3,4]
            y_series_title (list): names of y-series if multiple plots to go in legend
            x_title (str): label for x-axis
            y_title (str): label for y-axis
            figure_name (str): name for the output file
        """

        for index, item in enumerate(subplots):

            plt.rc("font", size=25)

            ax = plt.subplot(rows, columns, item)
            ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
            ax.yaxis.set_major_locator(ticker.MaxNLocator(10))
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))

            if index == 0 or index == 3:
                plt.scatter(
                    x, y[index], label=y_series_title[index], s=100, color="red"
                )
            elif index == 1 or index == 4:
                plt.scatter(
                    x, y[index], label=y_series_title[index], s=100, color="green"
                )
            elif index == 2 or index == 5:
                plt.scatter(
                    x, y[index], label=y_series_title[index], s=100, color="blue"
                )
            else:
                plt.scatter(
                    x, y[index], label=y_series_title[index], s=100, color="black"
                )

            if index >= 0 and index <= 2:
                ax.set_facecolor("wheat")
                ax.patch.set_alpha(0.3)
            elif index >= 3 and index <= 5:
                ax.set_facecolor("lavender")
                ax.patch.set_alpha(0.3)
            elif index == 6:
                ax.set_facecolor("silver")
                ax.patch.set_alpha(0.3)

            if index == 1:
                plt.title(title[0], fontsize=30, fontweight="bold")
            elif index == 4:
                plt.title(title[1], fontsize=30, fontweight="bold")
            elif index == 6:
                plt.title(title[2], fontsize=30, fontweight="bold")

            plt.xlabel(x_title)
            plt.ylabel(y_title[index])

            plt.legend(fontsize=30)

        figure = plt.gcf()
        figure.set_size_inches(35, 35)

        plt.savefig(figure_name, bbox_inches="tight", dpi=100)
        plt.clf()
        plt.close("all")

    def multi_multi_scatter_graph(
        self,
        x: list,
        y: dict,
        title: str,
        rows: int,
        columns: int,
        subplots: list,
        x_title: str,
        y_title: str,
        figure_name: str,
    ) -> None:

        """This function is for doing multiple separate scatter graphs,

        each with multiple series on them

        This works by having y as a dictionary of values

        (where each key is the series name)

        Example y is below:

        EXAMPLE y = {'Rint': [[dataset1], [dataset2]], 'R1': [[dataset1],[dataset2]], 'completeness': [[dataset1], [dataset2]]}

        Args:
            x (list): x-data
            y (dict): y-data
            title (str): title for graph
            rows (int): number of rows of graphs in figure
            columns (int): number of columns of graphs in figure
            subplots (list): list of graphs ie for 4 graphs it should be [1,2,3,4]
            y_series_title (list): names of y-series if multiple plots to go in legend
            x_title (str): label for x-axis
            y_title (str): label for y-axis
            figure_name (str): name for the output file
        """

        for index, i in enumerate(y):
            plt.rc("font", size=12)
            ax = plt.subplot(rows, columns, subplots[index])
            ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
            ax.yaxis.set_major_locator(ticker.MaxNLocator(10))
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3f"))
            for j in y[i]:
                if type(j) == list:
                    plt.scatter(x, j, label=i)
                else:
                    plt.scatter(x, y[i], label=i)
                    break
            plt.xlabel(x_title)
            if type(y_title) == str:
                plt.ylabel(y_title)
            else:
                plt.ylabel(y_title[index])
            plt.title(title)
            plt.legend()

        figure = plt.gcf()
        figure.set_size_inches(12, 20)

        plt.savefig(figure_name, bbox_inches="tight", dpi=100)
        plt.clf()
        plt.close("all")

    def four_line_graph(
        self,
        figure_name,
        x: list,
        y: list,
        x_title: list,
        y_title: list,
        full_title: str,
        mini_titles: list,
    ) -> None:

        """This is for doing a 2 x 2 grid of line graphs

        Args:
            figure_name (str): name for the output file
            x (list): x-data (list of list data - ie lists within a list, one for each subplot)
            y (list): y-data (list of list data - ie lists within a list, one for each subplot)
            x_title (list): list of labels for x-axis
            y_title (list: list of labels for y-axis
            full_title (str): title for graph
            mini_titles (list): titles for each of the 4 graphs

        """

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        plt.rc("font", size=20)
        ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3f"))
        ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3f"))
        ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3f"))
        ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3f"))
        fig.suptitle(full_title)
        if len(x[0]) == 1:
            ax1.scatter(x[0], y[0], color="r")
        else:
            ax1.plot(x[0], y[0], "r")
        ax1.set_xlabel(x_title[0])
        ax1.set_ylabel(y_title[0])
        ax1.set_title(mini_titles[0])

        if len(x[1]) == 1:
            ax2.scatter(x[1], y[1], color="b")
        else:
            ax2.plot(x[1], y[1], "b")
        ax2.set_xlabel(x_title[1])
        ax2.set_ylabel(y_title[1])
        ax2.set_title(mini_titles[1])

        if len(x[2]) == 1:
            ax3.scatter(x[2], y[2], color="g")
        else:
            ax3.plot(x[2], y[2], "g")
        ax3.set_xlabel(x_title[2])
        ax3.set_ylabel(y_title[2])
        ax3.set_title(mini_titles[2])

        if len(x[3]) == 1:
            ax4.scatter(x[3], y[3], color="tab:orange")
        else:
            ax4.plot(x[3], y[3], "tab:orange")
        ax4.set_xlabel(x_title[3])
        ax4.set_ylabel(y_title[3])
        ax4.set_title(mini_titles[3])

        figure = plt.gcf()
        figure.set_size_inches(25, 17)
        plt.savefig(figure_name, bbox_inches="tight", dpi=100)
        plt.clf()
        plt.close("all")

    def graph_multi_series(
        self, separated_dfs: list, x_title: str, y_title: list, figure_name: str
    ) -> None:

        """Multi-series graph where y is a list of dataframes separated by a parameter

        Args:
            separated_dfs (list): list of separated dataframes, where each item in list is a
                                    pd.DataFrame
            x_title (str): dataframe heading of x-axis
            y_title (list): list of dataframe labels for y-axis
            figure_name (str): name for the output file

        """

        for index, item in enumerate(separated_dfs):
            x = item[x_title]
            for index, param in enumerate(y_title):
                y = item[param]
                self.mini_graph(x, y, param, (index + 1), item["behaviour"][0])

        figure = plt.gcf()
        figure.set_size_inches(16, 12)
        plt.savefig(figure_name, bbox_inches="tight", dpi=100)
        plt.clf()
        plt.close("all")

    def mini_graph(
        self, x: list, y: list, title: str, subplot: int, behaviour: str
    ) -> None:

        """Pretty standard graphing function with nothing special

        Args:
            x (list): x-data
            y (list): y-data
            title (str): title for graph
            subplot (int): for defining subplot
            behaviour (str): label for scatter plot

        """

        ax = plt.subplot(2, 4, subplot)
        ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(10))
        plt.scatter(x, y, label=behaviour)
        plt.xlabel("Temperature(K)")
        if "alpha" or "beta" or "gamma" in title.lower():
            plt.ylabel("Angle(" + chr(176) + ")")
        elif "vol" in title.lower():
            plt.ylabel("Volume (\u212B\u00B3)")
        else:
            plt.ylabel(r"Distance ($\AA$)")
        plt.title(title)
        plt.legend()
        plt.close("all")


# ----------Class Definition----------#


class Cell_Import:
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

        self.config = Config(self.test_mode)

        self.cfg = self.config.cfg
        self.sys = self.config.sys
        self.conf_path = self.config.conf_path
        self.sys_path = self.config.sys_path

    def cell_import(self, ins: str) -> None:

        """Imports a unit cell from a reference .ins/.res file

        saves it in the sys.yaml file

        Args:
            ins (str): full path to the .ins/.res file
        """

        self.cfg, self.sys = self.config.yaml_reload(self.test_mode)

        # This function takes a .ins file and imports the unit cell parameters and saves it to the sys.yaml file

        with open(ins, "rt") as ins_file:
            split_line = []
            for line in ins_file:
                if "CELL" in line:
                    split_line = line.split()
                elif "TITL" in line and "REM" not in line:
                    self.sys["space_group"] = line.split()[3]

        ref_parameters = [
            "ref_a",
            "ref_b",
            "ref_c",
            "ref_alpha",
            "ref_beta",
            "ref_gamma",
        ]

        for index, item in enumerate(ref_parameters):
            self.sys[item] = float(split_line[index + 2])

        # Calculates the volume because it's not actually written anywhere

        self.sys["ref_volume"] = (
            self.sys["ref_a"]
            * self.sys["ref_b"]
            * self.sys["ref_c"]
            * math.sqrt(
                1
                - (math.cos(self.sys["ref_alpha"] * (math.pi / 180)) ** 2)
                - (math.cos(self.sys["ref_beta"] * (math.pi / 180)) ** 2)
                - (math.cos(self.sys["ref_gamma"] * (math.pi / 180)) ** 2)
                + (
                    2
                    * math.cos(self.sys["ref_alpha"] * (math.pi / 180))
                    * math.cos(self.sys["ref_beta"] * (math.pi / 180))
                    * math.cos(self.sys["ref_gamma"] * (math.pi / 180))
                )
            )
        )

        with open(self.sys_path, "w") as f:
            yaml.dump(
                self.sys,
                f,
                default_flow_style=False,
                Dumper=Nice_YAML_Dumper,
                sort_keys=False,
            )

        self.cfg, self.sys = self.config.yaml_reload(self.test_mode)

    def ref_edit(self, ins: str, MPLA_atoms: str) -> None:

        """Edits the reference .ins/.res file to put the MPLA

        command in with the user defined atoms

        Args:
            ins (str): full path to the .ins/.res file
            MPLA_atoms (str): list of atoms for MPLA command
        """

        with open(ins, "rt") as ins_file:
            content = ins_file.readlines()

        flag = False

        with open(ins, "w") as ins_file:
            for line in content:
                if "MPLA" in line:
                    flag = True
            if flag == False:
                for line in content:
                    if "PLAN" in line:
                        ins_file.write(line)
                        ins_file.write("MPLA " + MPLA_atoms + "\n")
                        ins_file.write("CONF\n")
                    else:
                        ins_file.write(line)
            else:
                for line in content:
                    ins_file.write(line)


# ----------Class Definition----------#


class File_Sorter:
    def sorted_properly(self, data: list) -> list:

        """A function to sort numbers like 1 and 10 properly

        Args:
            data (list): list to be sorted properly

        Returns:
            data (list): sorted list
        """

        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]

        return sorted(data, key=alphanum_key)


# ----------Class Definition----------#


class XDS_File_Edit:
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

    def change(self, file_path: str, parameter: str, new_value: str) -> None:

        """Changes a parameter in an XDS.INP file

        Args:
            file_path (str): full path to the XDS.INP file for editing
            parameter (str): the XDS.INP parameter that needs changing
            new_value (str): the new value of the XDS.INP parameter
        """

        # Honestly I forget why this is so complicated, but I remember having a lot of problems so it is the way that it is *shrugs in code*

        flag = 0
        editing = ""
        with open(file_path, "rt") as in_file:
            for line in in_file:
                if parameter in line and flag == 0:
                    to_edit = line.split()[1:]
                    flag += 1

        try:
            test = to_edit
        except UnboundLocalError:
            with open(file_path, "a") as f:
                f.write(" " + parameter + "= " + new_value)
        else:
            for element in to_edit:
                editing += " " + str(element)
            edit = editing.strip(" ")
            flag = 0
            for line in fileinput.input(file_path, inplace=True):
                if parameter in line and flag == 0:
                    line = line.rstrip("\r\n")
                    print(line.replace(edit, str(new_value)))
                    flag += 1
                else:
                    line = line.rstrip("\r\n")
                    print(line)

    def get_value(self, file_path: str, parameter: str) -> str:

        """Gets the value of a parameter in an XDS.INP

        Args:
            file_path (str): full path to the XDS.INP file for checking
            parameter (str): the parameter of interest

        Returns:
            value (str): the value of the specified parameter in the XDS.INP file
        """

        # Gets a value from the XDS INP file and saves it as a parameter

        self.value = ""
        with open(file_path, "rt") as in_file:
            for line in in_file:
                if parameter in line:
                    val_list = line.split()[1:]
        for element in val_list:
            self.value += "" + str(element)
        return self.value

    def start_angle(self, file_path: str) -> float:

        """Gets start and total angle from XDS.INP

        required later during editing to 'remember' what originally was

        I think it's a separate function because there's not just one value

        to the parameter (ie instead of frames = 100 its frames = 1 100)

        NOPE BECAUSE DEPENDENT ON THE EXPOSURE TIME

        Args:
            file_path (str): full path to the XDS.INP file

        Returns:
            angle (float): starting angle of the experiment
        """

        with open(file_path, "rt") as in_file:
            for line in in_file:
                if "STARTING_ANGLE= " in line:
                    angle = float(line.split()[1])

        return angle

    def new_line_rewrite(self, file_path: str) -> None:

        """To make editing the file easier later,

        the XDS.INP file is rewritten upon initialisation

        to put each keyword on a new line

        Args:
            file_path (str): full path to the XDS.INP file
        """

        with open(file_path, "r+") as in_file:
            lines = in_file.readlines()
        with open("temp.INP", "w") as out_file:
            for line in lines:

                # Multiple '=' signs indicate more than one keyword, EXCEPT for the headers used by XDS

                if line.count("=") > 1 and "!" not in line:
                    index_current = 0
                    index_old = 0
                    new_line = ""
                    for character in line:

                        # Splits the line at the spaces between different keywords/parameters (elif statement to make sure the last keyword isn't deleted due to a space NOT being present at the end of the line)

                        if (
                            character == " "
                            and line[index_current - 1] != " "
                            and line[index_current + 1] == " "
                        ):
                            new_line += (
                                " " + line[index_old : index_current + 1].strip() + "\n"
                            )
                            index_old = index_current
                        elif line.endswith(character) and len(new_line) != 0:
                            new_line += (
                                " " + line[index_old : index_current + 1].strip() + "\n"
                            )
                        index_current += 1
                    out_file.write(new_line)
                else:
                    out_file.write(line)
        shutil.move("temp.INP", file_path)


# ----------Class Definition----------#


class Nice_YAML_Dumper(yaml.SafeDumper):
    def write_line_break(self, data: dict = None) -> None:

        """Makes the yaml have better formatting when edited

        Unsure of the specifics, got this code from StackOverflow
        """

        super().write_line_break(data)

        if len(self.indents) == 1:
            super().write_line_break()


# ----------Class Definition----------#


class Config:
    def __init__(self, test_mode: bool = False) -> None:

        """Initialises the class

        Sets up the logbook and config file

        Args:
            test_mode (bool): Automatically false, if true it will

                            make the functions compatible with the testing script
        """

        self.conf_path = (
            pathlib.Path(os.path.abspath(__file__)).parent.parent / "conf.yaml"
        )
        self.sys_path = pathlib.Path(os.path.abspath(__file__)).parent / "sys.yaml"

        if test_mode == False:

            with open(self.conf_path, "r") as f:
                try:
                    self.cfg = yaml.load(f, yaml.FullLoader)
                except:
                    logging.critical(__name__ + " : Failed to open config file")
                    print("Error - See Log")
                    exit()

        if test_mode == True:
            self.cfg = "no conf.yaml file loaded"

        with open(self.sys_path, "r") as f:
            try:
                self.sys = yaml.load(f, yaml.FullLoader)
            except:
                logging.critical(__name__ + " : Failed to open system file")
                print("Error - See Log")
                exit()

    def yaml_reload(self, test_mode=False) -> Tuple[dict, dict]:

        """Reloads the yaml files

        Args:
            test_mode (bool): Automatically false, if true it will

                            make the functions compatible with the testing script

        Returns:
            cfg (dict): reloaded conf.yaml file
            sys (dict): reloaded sys.yaml file
        """

        if test_mode == False:
            with open(self.conf_path, "r") as f:
                self.cfg = yaml.load(f, yaml.FullLoader)

        with open(self.sys_path, "r") as f:
            self.sys = yaml.load(f, yaml.FullLoader)

        return self.cfg, self.sys


# ----------Class Definition----------#


class Generate:
    def __init__(self) -> None:

        """Initialises the class

        Used to generate the conf.yaml specific to each module/pipeline

        This basically just reads iin the parameter.yaml file and stores

        the data for access by external functions

        """
        ##Set up logbook and config file

        self.parameter_conf = (
            pathlib.Path(os.path.abspath(__file__)).parent / "parameter.yaml"
        )
        self.template_conf = (
            pathlib.Path(os.path.abspath(__file__)).parent / "blank_c.yaml"
        )

        with open(self.parameter_conf, "r") as f:
            try:
                self.param = yaml.load(f, yaml.FullLoader)
            except:
                logging.critical(__name__ + " : Failed to open parameter dictionary")
                print("Error - See Log")
                exit()


if __name__ == "__main__":
    parameters = set(Generate().param)
    required = set(
        {
            "module-platon-squeeze": ["file_name"],
            "pipeline-platon-squeeze": ["experiment_location"],
            "module-adp-analysis": ["csv_path", "cell_path"],
            "module-refinement": [
                "structure_location",
                "reference_path",
                "refinements_to_check",
                "tolerance",
                "maximum_cycles",
            ],
            "pipeline-refinement": [
                "experiment_location",
                "reference_path",
                "refinements_to_check",
                "tolerance",
                "maximum_cycles",
            ],
            "pipeline-variable-position": [
                "location_of_frames",
                "experiment_name",
                "reference_structure_location",
                "reference_XDS_INP_location",
                "reference_background_files_location",
                "instrument_parameters_path",
                "instrument_cif_path",
                "maximum_processors",
                "neggia_library",
                "space_group_number",
                "atoms_for_rotation_analysis",
                "reference_plane",
                "chemical_formula",
                "crystal_colour",
                "crystal_habit",
                "max_crystal_dimension",
                "middle_crystal_dimension",
                "min_crystal_dimension",
                "refinements_to_check",
                "tolerance",
                "maximum_cycles",
                "cif_parameters",
                "structural_analysis_bonds",
                "structural_analysis_angles",
                "structural_analysis_torsions",
                "ADP_analysis",
                "atoms_for_analysis",
                "wedge_angles",
                "min_pixels",
                "spot_maximum_centroid",
                "strong_pixels",
                "sepmin",
                "mapping_step_size",
                "total_angle",
            ],
            "pipeline-general": [
                "experiment_location",
                "reference_cif_location",
                "refinements_to_check",
                "tolerance",
                "maximum_cycles",
                "cif_parameters",
                "structural_analysis_bonds",
                "structural_analysis_angles",
                "structural_analysis_torsions",
                "ADP_analysis",
                "atoms_for_analysis",
                "varying_cif_parameter",
                "varying_parameter_values",
            ],
            "pipeline-general-extra": [
                "reference_location",
                "experiment_location",
                "_chemical_formula_moiety",
                "_chemical_absolute_configuration",
                "_exptl_absorpt_correction_T_max",
                "_exptl_absorpt_correction_T_min",
                "_exptl_absorpt_correction_type",
                "_exptl_crystal_colour",
                "_exptl_crystal_description",
                "_exptl_crystal_size_max",
                "_exptl_crystal_size_mid",
                "_exptl_crystal_size_min",
                "_diffrn_detector",
                "_diffrn_detector_area_resol_mean",
                "_diffrn_detector_type",
                "_diffrn_measurement_device",
                "_diffrn_measurement_device_type",
                "_diffrn_measurement_method",
                "_diffrn_radiation_monochromator",
                "_diffrn_radiation_probe",
                "_diffrn_radiation_type",
                "_diffrn_radiation_wavelength",
                "_diffrn_source",
                "_diffrn_source_type",
                "_computing_cell_refinement",
                "_computing_data_collection",
                "_computing_data_reduction",
                "_computing_molecular_graphics",
                "_computing_publication_material",
                "_computing_structure_refinement",
                "_computing_structure_solution",
                "refinements_to_check",
                "tolerance",
                "maximum_cycles",
                "cif_parameters",
                "structural_analysis_bonds",
                "structural_analysis_angles",
                "structural_analysis_torsions",
                "ADP_analysis",
                "atoms_for_analysis",
                "varying_cif_parameter",
                "varying_parameter_values",
            ],
            "module-intensity-compare": [
                "xprep_file_name",
                "h_condition_1",
                "h_condition_2",
                "k_condition_1",
                "k_condition_2",
                "l_condition_1",
                "l_condition_2",
                "a_axis",
                "b_axis",
                "c_axis",
                "alpha",
                "beta",
                "gamma",
                "h+k_condition_1",
                "h+k_condition_2",
                "h+l_condition_1",
                "h+l_condition_2",
                "k+l_condition_1",
                "k+l_condition_2",
                "h+k+l_condition_1",
                "h+k+l_condition_2",
                "include_h_0_condition_1",
                "include_h_0_condition_2",
                "include_k_0_condition_1",
                "include_k_0_condition_2",
                "include_l_0_condition_1",
                "include_l_0_condition_2",
            ],
            "pipeline-intensity-compare": [
                "data_location",
                "h_condition_1",
                "h_condition_2",
                "k_condition_1",
                "k_condition_2",
                "l_condition_1",
                "l_condition_2",
                "h+k_condition_1",
                "h+k_condition_2",
                "h+l_condition_1",
                "h+l_condition_2",
                "k+l_condition_1",
                "k+l_condition_2",
                "h+k+l_condition_1",
                "h+k+l_condition_2",
                "a_axis",
                "b_axis",
                "c_axis",
                "alpha",
                "beta",
                "gamma",
                "include_h_0_condition_1",
                "include_h_0_condition_2",
                "include_k_0_condition_1",
                "include_k_0_condition_2",
                "include_l_0_condition_1",
                "include_l_0_condition_2",
            ],
            "module-cif-merge": ["instrument_cif", "new_cif"],
            "module-make-instrument-cif": ["reference_cif"],
            "pipeline-cif-combine": ["location_of_cifs"],
            "pipeline-cif": [
                "experiment_location",
                "instrument_ending",
                "instrument_file",
                "structure_solution",
                "chemical_formula",
                "crystal_habit",
                "crystal_colour",
                "max_crystal_dimension",
                "middle_crystal_dimension",
                "min_crystal_dimension",
            ],
            "pipeline-rigaku-vt": [
                "experiment_location",
                "reference_location",
                "chemical_formula",
                "crystal_colour",
                "crystal_habit",
                "max_crystal_dimension",
                "middle_crystal_dimension",
                "min_crystal_dimension",
                "refinements_to_check",
                "tolerance",
                "maximum_cycles",
                "cif_parameters",
                "structural_analysis_bonds",
                "structural_analysis_angles",
                "structural_analysis_torsions",
                "ADP_analysis",
                "atoms_for_analysis",
            ],
            "pipeline-aus-synch-vt": [
                "location_of_autoprocess_folders",
                "experiment_name",
                "reference_location",
                "atoms_for_rotation_analysis",
                "reference_plane",
                "chemical_formula",
                "crystal_colour",
                "crystal_habit",
                "max_crystal_dimension",
                "middle_crystal_dimension",
                "min_crystal_dimension",
                "refinements_to_check",
                "tolerance",
                "maximum_cycles",
                "cif_parameters",
                "structural_analysis_bonds",
                "structural_analysis_angles",
                "structural_analysis_torsions",
                "ADP_analysis",
                "atoms_for_analysis",
            ],
            "module-xds-cell-transform": ["XDS_ASCII_File_1", "XDS_ASCII_File_2"],
            "module-xds-reprocess": [
                "xds_template_name",
                "XDS_INP_path",
                "experiment_location",
            ],
            "module-xprep": [
                "xprep_file_name",
                "transformation_matrix",
                "space_group",
                "chemical_formula",
                "experiment_location",
            ],
            "pipeline-xds-reprocess": ["XDS_INP_path", "experiment_location"],
            "pipeline-xprep": [
                "xprep_file_name",
                "transformation_matrix",
                "space_group",
                "chemical_formula",
                "experiment_location",
            ],
            "module-cell-analysis": [
                "csv_location",
                "reference_unit_cell",
                "x_axis_header",
            ],
            "module-cif-read": [
                "cif_parameters",
                "structural_analysis_bonds",
                "structural_analysis_angles",
                "structural_analysis_torsions",
                "ADP_analysis",
                "folder_containing_cifs",
            ],
            "module-rotation-planes": ["reference_plane", "lst_file_location"],
            "module-structural-analysis": [
                "atoms_for_analysis",
                "bond_data",
                "angle_data",
                "torsion_data",
            ],
            "pipeline-rotation-planes": ["reference_plane", "experiment_location"],
            "pipeline-variable-analysis": [
                "cif_parameters",
                "atoms_for_analysis",
                "structural_analysis_bonds",
                "structural_analysis_angles",
                "structural_analysis_torsions",
                "ADP_analysis",
                "varying_cif_parameter",
                "experiment_location",
                "reference_unit_cell",
            ],
            "pipeline-position-analysis": [
                "cif_parameters",
                "atoms_for_analysis",
                "reference_plane",
                "structural_analysis_bonds",
                "structural_analysis_angles",
                "structural_analysis_torsions",
                "ADP_analysis",
                "experiment_location",
                "reference_unit_cell",
                "atoms_for_plane",
                "wedge_angles",
                "min_pixels",
                "spot_maximum_centroid",
                "strong_pixels",
                "sepmin",
                "mapping_step_size",
            ],
            "pipeline-temperature-analysis": [
                "cif_parameters",
                "atoms_for_analysis",
                "structural_analysis_bonds",
                "structural_analysis_angles",
                "structural_analysis_torsions",
                "ADP_analysis",
                "experiment_location",
                "reference_unit_cell",
            ],
            "pipeline-AS-Brute": ["experiment_location"],
            "module-molecule-reconstruction": [
                "reference_path",
                "atom_list",
                "bond_list",
                "starting_atom",
                "starting_coordinates",
                "a_gradient",
                "b_gradient",
                "c_gradient",
                "alpha_gradient",
                "beta_gradient",
                "gamma_gradient",
                "max_position",
                "min_position",
                "position_step_size",
            ],
        }
    )
    setA = set(parameters)
    setB = set(required)
    setC = setB.intersection(setA)
    present = setC
    print(len(setA))
    print(len(setB))
    print(len(setC))
