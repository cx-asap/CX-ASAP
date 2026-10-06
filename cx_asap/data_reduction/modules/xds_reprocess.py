#!/usr/bin/env python3

###################################################################################################
# --------------------------------------CX-ASAP: XDS_reprocess-------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config, XDS_File_Edit
import shutil
import os
import pathlib
import logging

# ----------Class Definition----------#


class XDS_reprocess:
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

        self.XDS_edit = XDS_File_Edit()

    def run_xds(
        self,
        template_name: str,
        XDS_INP_path: str,
        location: str,
        structure_number: int,
    ) -> None:

        """Sets the correct frame name in the XDS_INP file

        Removes key XDS files from previous tests

        This avoids problems with mixing tests

        Runs XDS in the correct location

        Args:
            template_name (str): main name of the frames
            XDS_INP_path (str): full path to the XDS.INP file
            location (str): location where XDS will be run
            structure_number (int): keeps a tally of structures processed
        """

        # Makes a backup file, changes the name of the frames and copies it into the run folder

        shutil.copyfile(
            XDS_INP_path, pathlib.Path(pathlib.Path(XDS_INP_path).parent / "temp.INP")
        )

        os.chdir(location / "img")

        frame_type_flag = "h5"

        for item in os.listdir():
            if item.endswith(".img"):
                frame_type_flag = "img"

        if frame_type_flag == "h5":
            template_ending = "_??????.h5"
        elif frame_type_flag == "img":
            template_ending = "_???.img"

        os.chdir(location)

        self.XDS_edit.change(
            XDS_INP_path,
            "NAME_TEMPLATE_OF_DATA_FRAMES",
            "img/" + template_name + template_ending,
        )
        shutil.copy(XDS_INP_path, os.getcwd())
        os.remove(XDS_INP_path)
        os.rename(
            pathlib.Path(pathlib.Path(XDS_INP_path).parent / "temp.INP"), XDS_INP_path
        )

        # Removes any previous HKL files and indexation files from the system to make sure that if XDS fails, then old data isn't being run through XPREP

        items_to_remove = [
            "XPARM.XDS",
            "INTEGRATE.HKL",
            "SPOT.XDS",
            "BKGPIX.cbf",
            "ABS.cbf",
            "FRAME.cbf",
        ]

        for item in os.listdir():
            if ".HKL" in item:
                os.remove(item)
            elif item in items_to_remove:
                os.remove(item)

        # Runs XDS

        os.system("xds_par")

        # Updates error log with what happened :)

        for item in os.listdir():
            if item == "XDS_ASCII.HKL":
                logging.info(
                    __name__
                    + " : Data successfully indexed and integrated - structure "
                    + structure_number
                )
