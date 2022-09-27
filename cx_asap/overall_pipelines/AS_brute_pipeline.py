#!/usr/bin/env python3

###################################################################################################
# ----------------------------------------CX-ASAP: BRUTE PROCESSEING----------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config
from data_refinement.pipelines.refine_pipeline import Refinement_Pipeline
from tools.pipelines.shelx_t_pipeline import SHELXT_Pipeline
from data_reduction.pipelines.xprep_pipeline import XPREP_Pipeline
import os
import pathlib
import shutil
import logging

# ----------Class Definition----------#


class AS_Brute:
    def __init__(self) -> None:

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

        config = Config()

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    def initialise(self, experiment_location: str) -> None:

        """Checks if this pipeline has already been run

        Also moves into the correct directory

        Args:
            experiment_location (str): full path to the folder containing folders of AS data
        """

        os.chdir(experiment_location)
        if os.path.exists("Brute_Results") == True:
            shutil.rmtree("Brute_Results")

    def xprepreduce(self, location: str) -> None:

        """Will run XPREP with default settings on a series of data sets

        All data sets should be in individual folders in a common location

        Args:
            location (str): full path to the folder containing folders of AS data
        """

        xprep = XPREP_Pipeline()
        xprep.multiple_asdefaults(location)

    def solve(self, location: str) -> None:

        """Will run SHELXT with default settings on a series of data sets

        All data sets should be in individual folders in a common location

        Args:
            location (str): full path to the folder containing folders of AS data
        """

        shelxt = SHELXT_Pipeline()
        shelxt.multi_shelxt(location, "cxasap.ins")

    def report(self, location: str) -> None:

        """Checks if XDS/XPREP/SHELXT has worked for each data set analysed

        Outputs a report to the terminal and to a file

        Args:
            location (str): full path to the folder containing folders of AS data
        """

        if os.path.exists("Brute_Results") != True:
            os.mkdir("Brute_Results")

        XDS_fail = []
        XPREP_fail = []
        SHELXT_fail = []
        XDS_Notsad = []
        XPREP_Notsad = []
        SHELXT_Notsad = []

        os.chdir(location)
        for item in os.listdir():
            XDS_sad = False
            XPREP_sad = False
            SHELXT_sad = False
            if os.path.isdir(item) and item != "Brute_Results":
                os.chdir(item)
                for file in os.listdir():
                    if "XDS_ASCII.HKL_p1" == file:
                        XDS_sad = True
                        XDS_Notsad.append(item)
                    elif "cxasap.ins" == file:
                        XPREP_sad = True
                        XPREP_Notsad.append(item)
                    elif "cxasap_a.res" == file:
                        SHELXT_sad = True
                        SHELXT_Notsad.append(item)
                if XDS_sad == False:
                    XDS_fail.append(item)
                if XPREP_sad == False:
                    XPREP_fail.append(item)
                if SHELXT_sad == False:
                    SHELXT_fail.append(item)
                os.chdir(location)

        a = "------------------------------"
        b = "------AS_Brute Summary--------"
        c = "------------------------------"
        d = (
            "Successful Brutes: "
            + str(len(SHELXT_Notsad))
            + " of "
            + str(len(SHELXT_fail) + len(SHELXT_Notsad))
        )
        e = "Possible reasons for failure = XDS is Sad or XPREP is sad, or Jack put in the wrong formula"
        i = "Check the data reduction and the the actual frames!"
        g = "Successful Brutes:"
        h = "SHELXT Failed (aka Failed Brutes):"
        j = "XPREP Fail:"
        k = "XDS Fail:"
        nl = "\n"

        os.chdir(location + "/" + "Brute_Results")

        with open("Brute_summary.txt", "w") as f:
            f.write(str(a) + "\n")
            f.write(str(b) + "\n")
            f.write(str(c) + "\n")
            f.write(str(d) + "\n")
            f.write(str(e) + "\n")
            f.write(str(i) + "\n")
            f.write("\n")
            f.write(str(g) + "\n")

            for item in SHELXT_Notsad:
                try:
                    f.write(item + "\n")
                except:
                    f.write("ERROR" + "\n")
                    logging.info(__name__ + " : Folder without cxasap_a.res analysed")
            f.write("\n")
            f.write(str(h) + "\n")
            for item in SHELXT_fail:
                try:
                    f.write(item + "\n")
                except:
                    f.write("ERROR" + "\n")
                    logging.info(__name__ + " : Folder without cxasap_a.res analysed")
            f.write("\n")
            f.write(str(j) + "\n")
            
            for item in XPREP_fail:
                try:
                    f.write(item + "\n")
                except:
                    f.write("ERROR" + "\n")
                    logging.info(__name__ + " : Folder without cxasap.ins analysed")
            f.write("\n")
            f.write(str(k) + "\n")
            
            for item in XDS_fail:
                try:
                    f.write(item + "\n")
                except:
                    f.write("ERROR" + "\n")
                    logging.info(
                        __name__ + " : Folder without XDS_ASCII.HKL_p1 analysed"
                    )

        print(a)
        print(b)
        print(c)
        print(nl)
        print(d)
        print(e)
        print(i)
        print(nl)
        print(g)

        for item in SHELXT_Notsad:
            try:
                print(item)
            except:
                print("ERROR")
        print(nl)
        print(h)

        if len(SHELXT_fail) != 0:
            for item in SHELXT_fail:
                try:
                    print(item)
                except:
                    print("ERROR")
        else:
            print("None")
        print(a)

        print("CXASAP - saving you time since 1999")
