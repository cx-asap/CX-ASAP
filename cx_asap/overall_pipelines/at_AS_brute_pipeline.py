#!/usr/bin/env python3

####################################################################################################
# -------------------------------------CX-ASAP: BRUTE PROCESSEING----------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
####################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config
from tools.modules.shelx_t import SHELXT
from data_reduction.modules.xprep_module import XPREP
import os
import pathlib
import shutil
import logging

# ----------Class Definition----------#


class AS_Brute_Single:
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
        path = pathlib.Path(experiment_location)
        os.chdir(path.parent)
        if os.path.exists("Brute_Results") == False:
            os.mkdir("Brute_Results")

    def xprepreduce(self, location: str, formula: str = "C40H30N6FeCl2O8") -> None:

        """Will run XPREP with default settings on a series of data sets

        All data sets should be in individual folders in a common location

        Args:
            location (str): full path to the AS dataset for automatic xprepping
            formula (str): chemical formula
        """

        xprep = XPREP()
        xprep.asdefaults(location, formula)

    def solve(self, location: str) -> None:

        """Will run SHELXT with default settings on a series of data sets

        All data sets should be in individual folders in a common location

        Args:
            location (str): full path to the CX-ASAP_Brute folder containing xprepped data
        """

        shelxt = SHELXT()
        shelxt.run_shelxt(location, "cxasap")

        os.chdir(location)
        if os.path.exists("CX-ASAP_Brute") == False:
            os.mkdir("CX-ASAP_Brute")

        # Edit this for all shelxt files once get it actually running....

        files_to_copy = ["cxasap.hkl", "cxasap.ins", "cxasap.lxt", "cxasap.pcf", "cxasap_a.hkl", "cxasap_a.res"]

        for item in files_to_copy:

            try:
                shutil.move(item, "CX-ASAP_Brute/" + item)
            except FileNotFoundError:
                pass

    def report(self, location: str) -> None:

        """Checks if XDS/XPREP/SHELXT has worked for each data set analysed

        Outputs a report to the terminal and to a file

        Args:
            location (str): full path to the folder containing folders of AS data
        """
        
        XDS_fail = False
        XPREP_fail = False
        SHELXT_fail = False
        XDS_Notsad = False
        XPREP_Notsad = False
        SHELXT_Notsad = False
        
        os.chdir(pathlib.Path(location))
        
        if "XDS_ASCII.HKL_p1" in os.listdir():
            XDS_Notsad = True
        else:
            XDS_fail = True
        
        os.chdir(pathlib.Path(location)/"CX-ASAP_Brute")
        
        if "cxasap.hkl" in os.listdir():
            XPREP_Notsad = True
        else:
            XPREP_fail = True
            
        if "cxasap_a.res" in os.listdir():
            SHELXT_Notsad = True
        else:
            SHELXT_fail = True
                
        a = "------------------------------"
        b = "------AS_Brute Summary--------"
        c = "------------------------------"
        d = "Successful Brutes:"
        
        e = "Failed Brutes - The stage the analysis failed at is in brackets next to each dataset:"
        
        nl = "\n"
        
        os.chdir(pathlib.Path(location).parent/"Brute_Results")
        
        if os.path.exists("Successful_Brutes.txt") == False:
            with open ("Successful_Brutes.txt", "w") as f:
                f.write(a + nl)
                f.write(b + nl)
                f.write(c + nl)
                f.write(d + nl)
                
        if os.path.exists("Failed_Brutes.txt") == False:
            with open ("Failed_Brutes.txt", "w") as f:
                f.write(a + nl)
                f.write(b + nl)
                f.write(c + nl)
                f.write(e + nl)
                
        if SHELXT_Notsad == True: 
            with open ("Successful_Brutes.txt", "a") as f:
                f.write(pathlib.Path(location).name + nl)
                
        if SHELXT_fail == True and XDS_fail == False and XPREP_fail == False:
            with open ("Failed_Brutes.txt", "a") as f:
                f.write(pathlib.Path(location).name + "(SHELXT)" + nl)
                
        if XPREP_fail == True and XDS_fail == False:
            with open ("Failed_Brutes.txt", "a") as f:
                f.write(pathlib.Path(location).name + "(XPREP)" + nl)
                
        if XDS_fail == True:
            with open ("Failed_Brutes.txt", "a") as f:
                f.write(pathlib.Path(location).name + "(XDS)" + nl)
