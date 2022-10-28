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
import datetime

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

    def xprepreduce(self, location: str, file_name: str, formula: str = "C40H30N6FeCl2O8") -> None:

        """Will run XPREP with default settings on a series of data sets

        All data sets should be in individual folders in a common location

        Args:
            location (str): full path to the AS dataset for automatic xprepping
            file_name (str): name for xprep output
            formula (str): chemical formula
        """

        xprep = XPREP()
        if 'sadabs' in file_name:
            xprep.asdefaults_sadabs(location, formula, file_name)
        else:
            xprep.asdefaults(location, formula, file_name)

    def solve(self, location: str, file_name: str) -> None:

        """Will run SHELXT with default settings on a series of data sets

        All data sets should be in individual folders in a common location

        Args:
            location (str): full path to the CX-ASAP_Brute folder containing xprepped data
            file_name (str): file name for shelxt
        """

        shelxt = SHELXT()
        shelxt.run_shelxt(location, file_name)
            
    def move_files(self, parent_location: str, file_location: str, file_prefix: str) -> None:
        
        """Will move files into the main results folder
        
        Args:
            parent_location (str): full path to the raw dataset folder
            file_location (str): full path to the files analysed (ie might be 
                                in a sub folder post-sadabs
            file_prefix (str): name of files for copying into results folder
        """
        os.chdir(parent_location)
        if os.path.exists("CX-ASAP_Brute") == False:
            os.mkdir("CX-ASAP_Brute")
        
        os.chdir(file_location)
        
        files_to_copy = [file_prefix + ".hkl", file_prefix + ".ins", file_prefix + ".lxt", file_prefix + ".pcf", file_prefix + "_a.hkl", file_prefix + "_a.res"]

        for item in files_to_copy:
            try:
                shutil.move(item, parent_location + "/CX-ASAP_Brute/" + item)
            except FileNotFoundError:
                pass

    def report(self, location: str, suffix: str) -> None:

        """Checks if XDS/XPREP/SHELXT has worked for each data set analysed

        Outputs a report to the terminal and to a file

        Args:
            location (str): full path to the folder containing folders of AS data
            suffix (str): additional information to categorise results (ie to distinguish sadabs) 
        """
        
        XDS_fail = False
        XPREP_fail = False
        SHELXT_fail = False
        XDS_Notsad = False
        XPREP_Notsad = False
        SHELXT_Notsad = False
        timestamp = datetime.datetime.now().strftime("%c")
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
                f.write(timestamp + "_" + pathlib.Path(location).name + "_" + suffix + nl)
                
        if SHELXT_fail == True and XDS_fail == False and XPREP_fail == False:
            with open ("Failed_Brutes.txt", "a") as f:
                f.write(timestamp +  "_" + pathlib.Path(location).name + "_" + suffix + "(SHELXT)" + nl)
                
        if XPREP_fail == True and XDS_fail == False:
            with open ("Failed_Brutes.txt", "a") as f:
                f.write(timestamp +  "_" + pathlib.Path(location).name + "_" + suffix + "(XPREP)" + nl)
                
        if XDS_fail == True:
            with open ("Failed_Brutes.txt", "a") as f:
                f.write(timestamp + "_" + pathlib.Path(location).name + "_" + suffix  + "(XDS)" + nl)
