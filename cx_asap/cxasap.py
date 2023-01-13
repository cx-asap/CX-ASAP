#!/usr/bin/env python3

##################################################################################################
# --------------------------------CX-ASAP: Commandline Interface ---------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson Jack K. Clegg & Jason R. Price---#
# ---------------------------------Python Implementation by AJT-----------------------------------#
# ----------------------------------Project Design by JRP & JKC-----------------------------------#
# ------------------------------Valuable Coding Support by KMS & DJE------------------------------#
##################################################################################################

####################
# Overview of how the code works:

# There are a series of modules and pipelines available to run
# Modules perform individual tasks
# Pipelines perform multiple tasks
# There are generally two kinds of pipelines: job-specific and overall pipelines
# Job specific pipelines run a module over a series of structures, performing one job
# Overall pipelines string job specific pipelines together to perform a series of tasks
# The user sets variables using the conf.yaml file, which needs to be configured
# each time the user chooses a module/pipeline to run
# The code generates the conf.yaml depending on what code is chosen
# This avoids having one large file with many many parameters and makes it more user friendly
# Sometimes the code saves variables into the sys.yaml file - this simply allows it to
# to pass variables between modules
# The utils.py file contains a range of functions that are used by multiple modules
# This is simply done to avoid repetition where possible in the code
# If you are going to be developing modules, it is a good idea you understand this file
# so you know what options are available to you

####################

# How to develop this commandline interface:
# 1. import a new code that will be accessible from the commandline pipeline
# 2. add into parameter.yaml which params belong to this
# 3. add a new click command - note that need short help as well as a longer description
# 4. add options to that command
# 5. add the command to the group

#####################

###Template Click Command###

# """DESCRIBE MODULE/PIPELINE HERE
#
#    Most of these click functions are specifying text output to commandline
#    The main coding functions are checking the input of the yaml file
#    and setting up the corresponding class and calling its functions
#    Args:
#        The user will enter one of the four arguments as a flag
#        This will set that parameter as 'TRUE', while the others are 'FALSE'
#        This will define the value in the 'if/elif' statements
#        dependencies (bool): will check for dependencies
#        files (bool): will show the user what files are required
#        configure (bool): will set up the yaml for the user to fill out
#        run (bool): will execute the chosen module/pipeline
#    """

# @click.command('NAME HERE ', short_help='BRIEF DESCRIPTION HERE')
# @click.option('--dependencies', is_flag=True, help='view the software dependencies')
# @click.option('--files', is_flag=True, help='view the required input files')
# @click.option('--configure', is_flag=True, help='generate your conf.yaml file')
# @click.option('--run', is_flag=True, help='run the code!')
# def FUNCTION_NAME_HERE(dependencies, files, configure, run):
# """DESCRIBE MODULE/PIPELINE HERE
# """
# if dependencies:
###click.echo('\nYou require the below software in your path:')
###click.echo('- xxxxxxx')
# elif files:
###click.echo('\nYou require the below files:')
###click.echo(' - xxxxxxxxx')
###click.echo('\nThis folder can be located anywhere ')
# elif configure:
###click.echo('\nWriting a file called conf.yaml in the cx_asap folder...\n')
###click.echo('You will need to fill out the parameters.')
###click.echo('Descriptions are listed below:')
###click.echo(' - nnnnnnn: xxxx')

###fields = yaml_extraction('MODULE-NAME-HERE')
# yaml_creation(fields)

# elif run:
###click.echo('\nChecking to see if experiment configured....\n')

###check, cfg = configuration_check('MODULE-NAME-HERE')

# if check == False:
###click.echo('Make sure you fill in the configuration file!')
###click.echo('If you last ran a different code, make sure you reconfigure for the new script!')
###click.echo('Re-run configuration for description of each parameter\n')
# else:
###click.echo('READY TO RUN SCRIPT!\n')
# reset_logs()
# CREATE CLASS OBJECT AND RUN FUNCTIONS HERE
# copy_logs(OUTPUT_PATH_HERE)
# output_message()

# else:
###click.echo('Please select an option. To view options, add --help')


##############################################

import click
import yaml
import pathlib
import os
import subprocess
import logging
import platform
import time
import shutil
from typing import Union, Tuple

from system_files.utils import Generate, File_Sorter
from system_files.test_installation import Test
from data_reduction.modules.xprep_intensity_compare import Intensity_Compare
from data_reduction.modules.XDS_cell_transformation import XDS_Cell_Transformation
from data_reduction.modules.xds_reprocess import XDS_reprocess
from data_reduction.modules.xprep_module import XPREP
from data_reduction.pipelines.xprep_intensity_pipeline import Intensity_Pipeline
from data_reduction.pipelines.xds_pipeline import XDS_Pipeline
from data_reduction.pipelines.xprep_pipeline import XPREP_Pipeline
from data_refinement.modules.refinement import Structure_Refinement
from data_refinement.pipelines.refine_pipeline import Refinement_Pipeline
from cif_validation.modules.cif_merge import Cif_Merge
from cif_validation.modules.instrument_cif_generation import Instrument_CIF
from cif_validation.pipelines.cif_combine import Cif_Combine
from cif_validation.pipelines.cif_pipeline import CIF_Compile_Pipeline
from post_refinement_analysis.modules.cell_analysis import Cell_Deformation
from post_refinement_analysis.modules.cif_read import CIF_Read
from post_refinement_analysis.modules.rotation_planes import Rotation
from post_refinement_analysis.modules.structural_analysis import Structural_Analysis
from post_refinement_analysis.modules.ADP_analysis import ADP_analysis
from post_refinement_analysis.pipelines.rotation_pipeline import Rotation_Pipeline
from post_refinement_analysis.pipelines.variable_cif_parameter import (
    Variable_Analysis_Pipeline,
)
from post_refinement_analysis.pipelines.variable_position_analysis import (
    VP_Analysis_Pipeline,
)
from post_refinement_analysis.pipelines.variable_temperature_analysis import (
    VT_Analysis_Pipeline,
)
from overall_pipelines.variable_position_pipeline import VP_Pipeline
from overall_pipelines.cxasap_pipeline import General_Pipeline
from overall_pipelines.rigaku_synergy_vt_pipeline import Synergy_VT
from overall_pipelines.variable_temperature_pipeline import VT_Pipeline
from overall_pipelines.AS_brute_pipeline import AS_Brute
from overall_pipelines.at_AS_brute_pipeline import AS_Brute_Single
from tools.modules.platon_squeeze import Platon_Squeeze
from tools.pipelines.platon_squeeze_pipeline import Squeeze_Pipeline
from tools.modules.shelx_t import SHELXT
from tools.pipelines.shelx_t_pipeline import SHELXT_Pipeline
from tools.pipelines.pipeline_shelx_t_auto import SHELXT_Pipeline_auto
from tools.modules.molecule_reconstruction import Molecule_Reconstruction


def reset_logs() -> None:
    """Resets logs at the start of each run, keeps up to the 5 most previous logs - change quickly by editing the variable below max_logs
    Returns:
        None
    """

    log_location = (
        pathlib.Path(os.path.abspath(__file__)).parent / "error_logs/error_output.txt"
    )

    max_logs = 5
    number_list = []

    original_path = os.getcwd()
    os.chdir(log_location.parent)
    for item in os.listdir():
        if "error_output" in item:
            temp = item.split("_")
            try:
                number_list.append(1 + int(temp[2].strip(".txt")))
            except:
                number_list.append(1)

    number_list.sort()
    number_list.reverse()

    for item in number_list:
        if item > max_logs:
            number_list.remove(item)

    sorter = File_Sorter()
    files = sorter.sorted_properly(os.listdir())
    files.reverse()

    for item in files:
        for j in number_list:
            if str(j) in item and str(max_logs) not in item:
                os.rename(item, "error_output_" + str(j + 1) + ".txt")
            elif str(j) in item and str(max_logs) in item:
                os.remove(item)

    try:
        os.rename(log_location.name, "error_output_1.txt")
    except:
        pass

    os.chdir(original_path)
    
    try:
        os.remove(log_location)
    except:
        pass

def copy_logs(destination: str) -> None:
    """Copies logs to the output folder defined.
    If the module / pipeline outputs results, it goes to the results folder
    Otherwise, it goes to the experiment location folder
    Args:
        destination (string): path to the log destination
    """
    log_location = (
        pathlib.Path(os.path.abspath(__file__)).parent / "error_logs/error_output.txt"
    )

    log_location_1 = (
        pathlib.Path(os.path.abspath(__file__)).parent / "error_logs/error_output_1.txt"
    )

    try:
        shutil.copy(log_location, destination)
    except FileNotFoundError:
        shutil.copy(log_location_1, destination)


def yaml_extraction(heading: str) -> dict:
    """Extracts the required yaml parameters for the chosen code.
    Not all yaml parameters are required for every module/pipeline
    The required parameters for each module/pipeline are separated by heading
    This function makes a conf.yaml with only the parameters under this input heading
    This function also serves to set default parameters where possible
    Args:
        heading (str): Corresponds to the module or pipeline name that is being configured
    Returns:
        yaml_dict (dict): A dictionary of the yaml paramters for chosen module/pipeline
    """
    configure = Generate()
    list_of_params = []
    yaml_dict = {}

    for item in configure.param:
        if item == heading:

            # if type(configure.param[item]) == dict:
            # for item2 in configure.param[item]:
            # if type(configure.param[item][item2]) == dict:
            # for item3 in configure.param[item][item2]:
            # list_of_params = (
            # list_of_params + configure.param[item][item2][item3]
            # )
            # else:
            # list_of_params = list_of_params + configure.param[item][item2]

            # else:
            list_of_params = list_of_params + configure.param[item]

    list_params = [
        "wedge_angles",
        "min_pixels",
        "spot_maximum_centroid",
        "strong_pixels",
        "sepmin",
        "atoms_for_analysis",
        "varying_parameter_values",
        "atom_list",
        "bond_list",
    ]

    structure_params = [
        "structural_analysis_angles",
        "structural_analysis_bonds",
        "structural_analysis_torsions",
        "ADP_analysis",
    ]

    for item in list_of_params:
        if item == "cif_parameters":
            yaml_dict[item] = [
                "_cell_length_a",
                "_cell_length_b",
                "_cell_length_c",
                "_cell_angle_alpha",
                "_cell_angle_beta",
                "_cell_angle_gamma",
                "_cell_volume",
                "_diffrn_reflns_av_R_equivalents",
                "_diffrn_measured_fraction_theta_full",
                "_diffrn_ambient_temperature",
                "_refine_ls_R_factor_gt",
            ]
        elif item == "reference_plane" or item == "starting_coordinates":
            yaml_dict[item] = [0, 0, 0]
        elif item in list_params:
            yaml_dict[item] = [0]
        elif item == "refinements_to_check":
            yaml_dict[item] = 8
        elif item == "tolerance":
            yaml_dict[item] = 0.002
        elif item == "transformation_matrix":
            yaml_dict[item] = "1 0 0 0 1 0 0 0 1"
        elif item == "maximum_cycles":
            yaml_dict[item] = 20
        elif item in structure_params:
            yaml_dict[item] = False
        else:
            yaml_dict[item] = 0

    return yaml_dict


def yaml_creation(yaml_dict: dict) -> None:
    """Converts the yaml_dict made from the yaml_extraction function into a yaml file
    This is the conf.yaml that the user will edit to configure the code
    Args:
        yaml_dict (dict): Dictionary of yaml parameters specific to the chosen module/pipeline
    """
    #yaml_path = pathlib.Path(os.path.abspath(__file__)).parent / "conf.yaml"
    yaml_path = pathlib.Path(os.path.join(os.getcwd()), "conf.yaml")

    with open(yaml_path, "w") as f:
        new_yaml = yaml.dump(yaml_dict, f)


def configuration_check(heading: str) -> Tuple[bool, dict]:
    """Checks that there are no code-breaking syntax errors in the conf.yaml
    This is run after the user has entered their parameters
    and also before the requested module/pipeline is called
    This is also used as a check to make sure the user has generated the yaml file
    and also that the yaml file matches the headings required for the chosen module/pipeline
    For example, if the user has changed modules but hasn't reconfigured the yaml
    then this check will throw an error.
    Args:
        heading (str): Corresponds to the module or pipeline name that is being configured

    Returns:
        flag (bool): TRUE if passes checks, FALSE if fails
        cfg (dict): The user's edits in the conf.yaml file imported as a dictionary
    """
    yaml_dict = yaml_extraction(heading)

    #yaml_path = pathlib.Path(os.path.abspath(__file__)).parent / "conf.yaml"
    yaml_path =  pathlib.Path(os.path.join(os.getcwd(),"conf.yaml"))

    zero_exceptions = [
        "a_gradient",
        "alpha_gradient",
        "b_gradient",
        "beta_gradient",
        "c_gradient",
        "gamma_gradient",
    ]

    if heading == "pipeline-AS-Brute-individual":
        zero_exceptions.append("chemical_formula")
    
    if os.path.exists(yaml_path) == False:
        click.echo("No configuration file, please run with --configure option\n")
    else:
        with open(yaml_path, "r") as f:
            try:
                cfg = yaml.load(f, yaml.FullLoader)
            except:
                click.echo("Failed to set up config file. Try reconfiguring")
                exit()

        flag = True

        for item in yaml_dict.keys():
            if item not in cfg.keys():
                flag = False

        if flag == True:

            for item in cfg:
                if cfg[item] == 0 and item not in zero_exceptions:
                    flag = False

                    if cfg[item] == False and type(cfg[item]) == bool:
                        flag = True
                    else:
                        break

                if item == "atoms_for_analysis":
                    for j in cfg[item]:
                        try:
                            len(j.split(" "))
                        except AttributeError:
                            if j == 0:
                                flag = False
                        else:
                            if len(j.split(" ")) > 1:
                                print(
                                    "Error in reading config file. Make sure your atoms_for_analysis are on separate lines in list format. The pre-set dash indicates how to type each atom on a new line"
                                )
                                exit()

    return flag, cfg


@click.group()
def cli():
    """#####################################################################\n
    Welcome to CX-ASAP! Please read the below text to get started.\n
    #####################################################################\n

    How to use CX-ASAP:\n
    1. If you have just installed this software, run the 'test' command \n
    2. Choose a pipeline/module from the commands below \n
    3. View --dependencies and check that you satisfy them \n
    4. View --files and check your directory structure satisfies the requirements \n
    5. Run --configure to write your yaml file \n
    6. Open your newly generated yaml file and fill in the required parameters \n
    7. Run your chosen pipeline/module using --run \n

    #####################################################################\n

    Executing commands:\n
    Practical Usage - 'cxasap COMMAND [ARGS]'\n
    Full Usage - 'cxasap [OPTIONS] COMMAND [ARGS]'\n
    Note that the only option is '--help' and will display this home message\n
    You will mostly be using the 'Practical Usage' of CX-ASAP\n
    The commands are the different modules/pipelines available.\n
    For example, to configure the refinement pipeline, you would type:\n
    'cxasap pipeline-refinement --configure'

    #####################################################################\n

    Choosing a module/pipeline:
    For more information on each available module/pipeline, run the below command:\n
    'cxasap module/pipeline --help' - where module/pipeline is in the below list of commands.\n

    #####################################################################\n

    To access error logs, run the 'errors' command (ie 'cxasap errors')\n
    This will print the log files into the terminal.\n
    Otherwise, they can be found in one of two places:\n
    They will be output with your data processing\n
    And the most recent 5 logs will be found in cx_asap/error_logs

    #####################################################################\n

    Disclaimer:
    What you put in is what you get out!
    If your reference files are incomplete - your output files will be incomplete.
    If your reference structure does not model disorder, none of your output files will model disorder.
    Consider your reference files carefully!!!

    #####################################################################\n

    If you require further help or assistance in setting up required software, consult the documentation \n

    #####################################################################\n

    You are currently running version 1.1.1

    #####################################################################\n

    """


##########-Test Command-##########


@click.command("test", short_help="Test installation of cxasap")
def test() -> None:
    """Executes the test script
    The test script is found in the system files folder
    It is recommended to users that this is run upon installation
    """
    tester = Test()
    tester.test_available_code()


##########-Error Command-##########


@click.command("errors", short_help="Display current error log")
def errors() -> None:
    """Outputs the error logs to the terminal
    Otherwise the user will need to find the .txt file in the documentation folder
    """
    error_path = (
        pathlib.Path(os.path.abspath(__file__)).parent
        / "error_logs"
        / "error_output.txt"
    )

    if os.path.exists(error_path):
        with open(error_path, "rt") as f:
            for line in f:
                click.echo(line)


##########-Output Completion Message-##########


def output_message() -> None:
    """Prints a message upon successful completion of a module/pipeline"""
    message = """
    #################################################################################\n
    CX-ASAP completed!\n
    Please review all statistics and examine output files prior to publication.\n
    If you found this software useful, please consider citing it in your publication.\n
    Thompson, A. J., Smith, K. M. L., Clegg, J. K., Price, J. R. (2022),\n 
    ChemRxiv. Preprint. https://doi.org/10.26434/chemrxiv-2022-7c0cl\n
    Thank you for using this software, and may the Fourier Transforms be with you. \n
    #################################################################################
    """
    
    print(message)

##########-OS Test-##########


def os_test() -> None:
    """This function has not been fully developed
    Aimed to change what functions were available depending on what operating system
    Also a template because xprep functions may be different for Unix vs Windows
    But! The below code had problems with the rest of the click interface
    So this has been left for future testing and development
    """
    os_name = platform.system()
    BadOS = False
    ###print('Testing for operating system...')
    # time.sleep(2)
    ###print('Your OS is: ' + os_name)
    # time.sleep(1)
    if os_name == "Windows":
        BadOS = True
        ###print('Your battlestation is not fully operational...')
        ###print('CX-ASAP loading...')
    # else:
    ###print('Your battlestation is fully operational!')
    ###print('CX-ASAP loading...')

    # time.sleep(2)


##########-Module Refinement-##########

"""This module will refine a single structure to convergence using a 
    single reference file.   
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag       
        This will set that parameter as 'TRUE', while the others are 'FALSE'
        This will define the value in the 'if/elif' statements
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""


@click.command("module-refinement", short_help="Single Structure Refinement")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def module_refinement(dependencies, files, configure, run):
    """This module will refine a single structure to convergence using a
    single reference file.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- SHELXL\n")

    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a reference .ins or .res file")
        click.echo(" - a .hkl file")
        click.echo(" - a .ins file corresponding to the .hkl file")
        click.echo("\nYour .hkl/.ins file should be in the same folder")
        click.echo("Your reference .ins/.res file can be located anywhere\n")

    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(" - maximum_cycles: enter the max number of cycles shelxl can run")
        click.echo(
            " - reference_path: enter the full path to your reference .ins or .res file"
        )
        click.echo(
            " - refinements_to_check: enter the number of refinements you want to check for convergence"
        )
        click.echo(" - structure_location: enter the full path to your new .ins file")
        click.echo(
            " - tolerance: enter the desired mean shift for the number of refinements in refinements_to_check\n"
        )

        fields = yaml_extraction("module-refinement")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("module-refinement")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")

            reset_logs()
            shelxl = Structure_Refinement()
            shelxl.run_shelxl(
                cfg["structure_location"],
                cfg["reference_path"],
                cfg["refinements_to_check"],
                cfg["tolerance"],
                cfg["maximum_cycles"],
            )

            copy_logs(pathlib.Path(cfg["structure_location"]).parent)

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


##########-Pipeline Refinement-##########

"""This module will refine a series of structures to convergence 
    based on a single reference file. Such a dataset might have come from a    
    dynamic experiment, such as variable-temperature, variable-pressure, 
    variable-position, etc
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag 
        This will set that parameter as 'TRUE', while the others are 'FALSE'
        This will define the value in the 'if/elif' statements
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""


@click.command("pipeline-refinement", short_help="Multi Structure Refinement")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_refinement(dependencies, files, configure, run):
    """This module will refine a series of structures to convergence
    based on a single reference file. Such a dataset might have come from a
    dynamic experiment, such as variable-temperature, variable-pressure,
    variable-position, etc
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- SHELXL\n")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a reference .ins or .res file")
        click.echo(" - a series of .hkl files")
        click.echo(" - a series .ins files corresponding to the .hkl files")
        click.echo(
            "\nYour .hkl/.ins files should be in separate folders located in a single parent folder"
        )
        click.echo(
            "Your reference .ins/.res file should be located outside of the individual .hkl/.ins folders\n"
        )
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - experiment_location: enter the full path to the folder containing all dataset folders"
        )
        click.echo(" - maximum_cycles: enter the max number of cycles shelxl can run")
        click.echo(
            " - reference_path: enter the full path to your reference .ins or .res file"
        )
        click.echo(
            " - refinements_to_check: enter the number of refinements you want to check for convergence"
        )
        click.echo(
            " - tolerance: enter the desired mean shift for the number of refinements in refinements_to_check\n"
        )

        fields = yaml_extraction("pipeline-refinement")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-refinement")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            shelxl = Refinement_Pipeline()
            shelxl.multiple_refinement(
                cfg["experiment_location"],
                cfg["reference_path"],
                cfg["experiment_location"],
                cfg["refinements_to_check"],
                cfg["tolerance"],
                cfg["maximum_cycles"],
            )

            copy_logs(cfg["experiment_location"])
        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


##########-Pipeline VP-##########

"""This pipeline is the full analysis for flexible crystal mapping experiments.
    Starts with .h5 diffraction frames and outputs finalised CIF files and analysis.    
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag 
        This will set that parameter as 'TRUE', while the others are 'FALSE'
        This will define the value in the 'if/elif' statements
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""


@click.command(
    "pipeline-variable-position", short_help="Full Variable Position Pipeline"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_vp(dependencies, files, configure, run):
    """This pipeline is the full analysis for flexible crystal mapping experiments.
    Starts with .h5 diffraction frames and outputs finalised CIF files and analysis.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- SHELXL")
        click.echo("- platon")
        click.echo("- XDS (parallel version, with dectris-neggia libraries downloaded)")
        click.echo("- xprep\n")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a reference .ins or .res file")
        click.echo(" - a CIF file with all instrument parameters")
        click.echo(" - BKGINIT.cbf from reference")
        click.echo(" - BLANK.cbf from reference")
        click.echo(" - GAIN.cbf from reference")
        click.echo(" - GXPARM.XDS from reference")
        click.echo(" - X-CORRECTIONS.cbf from reference")
        click.echo(" - Y-CORRECTIONS.cbf from reference")
        click.echo(
            " - XDS.INP from reference - make sure the data/spot ranges are correct"
        )
        click.echo(" - a series of frames with a common root name")
        click.echo("\nYour .h5 frames should be in a common folder")
        click.echo(
            "BKGINIT.cbf, BLANK.cbf, GAIN.cbf, X-CORRECTIONS.cbf and Y-CORRECTIONS.cbf should be in a folder together"
        )
        click.echo("Your reference files can be located anywhere\n")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - ADP_analysis: enter True for ADP analysis, otherwise enter False"
        )
        click.echo(
            " - atoms_for_analysis: enter the atom labels for graphical structural analysis as a list (best suited to small numbers to avoid over cluttering graphs"
        )
        click.echo(
            " - atoms_for_rotation_analysis: enter the labels of the atoms for mean plane analysis"
        )
        click.echo(" - chemical_formula: enter the chemical formula of your crystal")
        click.echo(
            " - cif_parameters: enter the parameters to be extracted from your cifs - the default parameters are usually fine"
        )
        click.echo(" - crystal_colour: enter the colour of your crystal")
        click.echo(" - crystal_habit: enter the habit of your crystal")
        click.echo(
            " - experiment_name: enter the common name between frames (ie AJT-01-018)"
        )
        click.echo(
            " - instrument_cif_path: enter the full path to your instrument cif file"
        )
        click.echo(
            " - instrument_parameters_path: enter the full path to your reference GXPARM.XDS file"
        )
        click.echo(
            " - location_of_frames: enter the full path to the folder containing your .h5 frames"
        )
        click.echo(
            " - mapping_step_size: enter the step-size for the collection (in um)"
        )
        click.echo(
            " - max_crystal_dimension: enter the largest dimension of your crystal (in mm)"
        )
        click.echo(" - maximum_cycles: enter the max number of cycles shelxl can run")
        click.echo(
            " - maximum_processors: enter the maximum number of processors you want XDS to use"
        )
        click.echo(
            " - middle_crystal_dimension: enter the middle diimension of your crystal (in mm)"
        )
        click.echo(
            " - min_crystal_dimension: enter the smallest dimension of your crystal (in mm)"
        )
        click.echo(
            " - min_pixels: enter the values for minimum_pixels_in_a_spot as a list for XDS proessing"
        )
        click.echo(" - neggia_library: enter the full path to dectris_neggia.so")
        click.echo(
            " - reference_background_files_location: enter the full path to the folder containing BKGINIT.cbf, BLANK.cbf, GAIN.cbf, X-CORRECTIONS.cbf and Y-CORRECTIONS.cbf"
        )
        click.echo(
            " - reference_plane: enter the reference plane for mean plane analysis as a list"
        )
        click.echo(
            " - reference_structure_location: enter the full path to your reference .ins/.res file"
        )
        click.echo(
            " - reference_XDS_INP_location: enter the full path to your reference XDS.INP file"
        )
        click.echo(
            " - refinements_to_check: enter the number of refinements you want to check for convergence"
        )
        click.echo(
            " - sepmin: enter the values for sepmin as a list for XDS processing"
        )
        click.echo(
            " - space_group_number: enter the space group number for your compound (ie 14)"
        )
        click.echo(
            " - spot_maximum_centroid: enter the values for spot_maximum_centroid as a list for XDS processing"
        )
        click.echo(
            " - strong_pixels: enter the values for strong_pixels as a list for XDS processing"
        )
        click.echo(
            " - structural_analysis_bonds: enter True for bond length analysis, otherwise enter False"
        )
        click.echo(
            " - structural_analysis_angles: enter True for angle analysis, otherwise enter False"
        )
        click.echo(
            " - structural_analysis_torsions: enter True for torsion analysis, otherwise enter False - note that this will have required the CONF command in your reference .ins/.res file"
        )
        click.echo(
            " - tolerance: enter the desired mean shift for the number of refinements in refinements_to_check"
        )
        click.echo(
            " - wedge_angles: enter the wedge angles as a list for XDS processing"
        )

        fields = yaml_extraction("pipeline-variable-position")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-variable-position")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            full_vp_analysis = VP_Pipeline()
            full_vp_analysis.setup(
                cfg["location_of_frames"],
                cfg["experiment_name"],
                "flexible_mapping_analysis",
                cfg["reference_structure_location"],
                cfg["reference_XDS_INP_location"],
                cfg["reference_background_files_location"],
                cfg["instrument_parameters_path"],
                cfg["maximum_processors"],
                cfg["neggia_library"],
                cfg["space_group_number"],
                cfg["atoms_for_rotation_analysis"],
                cfg["instrument_cif_path"],
                cfg["total_angle"],
            )

            full_vp_analysis.flexible_parameter_loops(
                "1 0 0 0 1 0 0 0 1",
                full_vp_analysis.sys["analysis_path"],
                "XDS_ASCII.HKL",
                full_vp_analysis.sys["space_group"],
                cfg["chemical_formula"],
                full_vp_analysis.sys["ref_path_organised"],
                full_vp_analysis.sys["current_results_path"],
                cfg["refinements_to_check"],
                cfg["tolerance"],
                cfg["maximum_cycles"],
                cfg["reference_plane"],
                full_vp_analysis.sys["XDS_inp_organised"],
                cfg["wedge_angles"],
                cfg["min_pixels"],
                cfg["sepmin"],
                cfg["spot_maximum_centroid"],
                cfg["strong_pixels"],
                "known_structure",
                cfg["crystal_habit"],
                cfg["crystal_colour"],
                cfg["max_crystal_dimension"],
                cfg["middle_crystal_dimension"],
                cfg["min_crystal_dimension"],
                False,
                pathlib.Path(cfg["instrument_cif_path"]).name,
            )

            full_vp_analysis.analyse(
                full_vp_analysis.sys["current_results_path"],
                cfg["structural_analysis_bonds"],
                cfg["structural_analysis_angles"],
                cfg["structural_analysis_torsions"],
                cfg["ADP_analysis"],
                cfg["cif_parameters"],
                cfg["atoms_for_analysis"],
                full_vp_analysis.sys["ref_path_organised"],
                cfg["atoms_for_rotation_analysis"],
                cfg["mapping_step_size"],
                cfg["spot_maximum_centroid"],
                cfg["min_pixels"],
                cfg["strong_pixels"],
                cfg["sepmin"],
                cfg["wedge_angles"],
                cfg["reference_plane"],
            )

            copy_logs(full_vp_analysis.sys["current_results_path"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


##########-Pipeline General-##########

"""This pipeline will refine a series of .hkl/.ins files based on a common  
    reference and provide finalised .cif files, checkCIF reports and structural 
    analysis. Use this only if your files are already in the correct directory 
    structure (check by adding --files to command). Otherwise, consider 
    'pipeline-general-extra'
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag 
        This will set that parameter as 'TRUE', while the others are 'FALSE'
        This will define the value in the 'if/elif' statements
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""


@click.command(
    "pipeline-general", short_help="General refinement->finalisation pipeline"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_general(dependencies, files, configure, run):
    """This pipeline will refine a series of .hkl/.ins files based on a common
    reference and provide finalised .cif files, checkCIF reports and structural
    analysis. Use this only if your files are already in the correct directory
    structure (check by adding --files to command). Otherwise, consider
    'pipeline-general-extra'
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- SHELXL")
        click.echo("- platon")
        click.echo("- shredCIF")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a reference .cif file")
        click.echo(" - a series of .hkl files")
        click.echo(" - a series .ins files corresponding to the .hkl files")
        click.echo(
            "\nYour .hkl/.ins files should be in separate folders located in a single parent folder"
        )
        click.echo(
            "Your reference .cif file should be located anywhere outside of the individual .hkl/.ins folders"
        )
        click.echo(
            "Your reference .cif file should contain all required crystal/instrument parameters\n"
        )
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - ADP_analysis: enter True for ADP analysis, otherwise enter False"
        )
        click.echo(
            " - atoms_for_analysis: enter the atom labels for graphical structural analysis as a list (best suited to small numbers to avoid over cluttering graphs"
        )
        click.echo(
            " - cif_parameters: enter the parameters to be extracted from your cifs - the default parameters are usually fine"
        )
        click.echo(
            " - experiment_location: enter the full path to the folder containing your data sets"
        )
        click.echo(" - maximum_cycles: enter the max number of cycles shelxl can run")
        click.echo(
            " - reference_cif_location: enter the full path to your reference .cif file"
        )
        click.echo(
            " - refinements_to_check: enter the number of refinements you want to check for convergence"
        )
        click.echo(
            " - structural_analysis_bonds: enter True for bond length analysis, otherwise enter False"
        )
        click.echo(
            " - structural_analysis_angles: enter True for angle analysis, otherwise enter False"
        )
        click.echo(
            " - structural_analysis_torsions: enter True for torsion analysis, otherwise enter False - note that this will have required the CONF command in your reference .ins/.res file"
        )
        click.echo(
            " - tolerance: enter the desired mean shift for the number of refinements in refinements_to_check"
        )
        click.echo(
            " - varying_cif_parameter: enter the parameter in your cif files that is varying in proper cif syntax (ie _diffrn_ambient_temperature)"
        )
        click.echo(
            " - varying_parameter_values: enter the values of your varying parameter as a list with errors. Ie if you did a VT experiment at 100(2)K and 200(2)K, enter 100(2) and 200(2) in the provided list"
        )

        fields = yaml_extraction("pipeline-general")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-general")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            full = General_Pipeline()

            full.make_dirs(cfg["experiment_location"])

            full.reference_extract(
                cfg["reference_cif_location"],
                cfg["experiment_location"],
                cfg["varying_parameter_values"],
                cfg["varying_cif_parameter"],
            )

            full.process(
                cfg["experiment_location"],
                full.reference_res,
                full.stats_location,
                cfg["refinements_to_check"],
                cfg["tolerance"],
                cfg["maximum_cycles"],
            )

            full.analyse(
                full.reference_res,
                cfg["experiment_location"],
                full.results_location,
                "cx-asap",
                full.chemical_formula,
                full.crystal_habit,
                full.crystal_colour,
                full.max_dimension,
                full.mid_dimension,
                full.min_dimension,
                cfg["structural_analysis_bonds"],
                cfg["structural_analysis_angles"],
                cfg["structural_analysis_torsions"],
                cfg["cif_parameters"],
                cfg["atoms_for_analysis"],
                cfg["varying_cif_parameter"],
                False,
                "instrument.cif",
                True,
                cfg["ADP_analysis"],
            )

            copy_logs(full.results_location)

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


##########-Pipeline General with Extra Setup-##########

"""This pipeline will refine a series of .hkl/.ins files based on a common   
    reference and provide finalised .cif files, checkCIF reports and structural analysis.
    This version of the script will set up folders for you based on your configuration file.
    You will be prompted to place your structure files in the folders created by the code.
    This version of the script also runs off the assumption you do not have a completed 
    reference cif with all instrument parameters. Instead, you will be able to enter   
    them into the configuration file.
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'
        This will define the value in the 'if/elif' statements
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
    """


@click.command("pipeline-general-extra", short_help="General pipeline with extra setup")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_general_extra(dependencies, files, configure, run):
    """This pipeline will refine a series of .hkl/.ins files based on a common
    reference and provide finalised .cif files, checkCIF reports and structural analysis.
    This version of the script will set up folders for you based on your configuration file.
    You will be prompted to place your structure files in the folders created by the code.
    This version of the script also runs off the assumption you do not have a completed
    reference cif with all instrument parameters. Instead, you will be able to enter
    them into the configuration file.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- SHELXL")
        click.echo("- platon")
        click.echo("- shredCIF")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a reference .ins or .res file")
        click.echo(" - a series of .hkl files")
        click.echo(" - a series .ins files corresponding to the .hkl files")
        click.echo(
            "\nYour .hkl/.ins files can be anywhere, but you will need to put them in the folders the code creates"
        )
        click.echo(
            "Your reference .ins/.res file should be located anywhere outside of the individual .hkl/.ins folders"
        )
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - atoms_for_analysis: enter the atom labels for graphical structural analysis as a list (best suited to small numbers to avoid over cluttering graphs"
        )
        click.echo(
            " - cif_parameters: enter the parameters to be extracted from your cifs - the default parameters are usually fine"
        )
        click.echo(
            " - experiment_location: enter the full path to the folder containing your data sets"
        )
        click.echo(" - maximum_cycles: enter the max number of cycles shelxl can run")
        click.echo(
            " - reference_location: enter the full path to your reference .ins/.res file"
        )
        click.echo(
            " - refinements_to_check: enter the number of refinements you want to check for convergence"
        )
        click.echo(
            " - structural_analysis_bonds: enter True for bond length analysis, otherwise enter False"
        )
        click.echo(
            " - structural_analysis_angles: enter True for angle analysis, otherwise enter False"
        )
        click.echo(
            " - structural_analysis_torsions: enter True for torsion analysis, otherwise enter False - note that this will have required the CONF command in your reference .ins/.res file"
        )
        click.echo(
            " - tolerance: enter the desired mean shift for the number of refinements in refinements_to_check"
        )
        click.echo(
            " - varying_cif_parameter: enter the parameter in your cif files that is varying in proper cif syntax (ie _diffrn_ambient_temperature)"
        )
        click.echo(
            " - varying_parameter_values: enter the values of your varying parameter as a list with errors. Ie if you did a VT experiment at 100(2)K and 200(2)K, enter 100(2) and 200(2) in the provided list\n"
        )
        click.echo("IMPORTANT NOTE: there are also a series of cif parameters listed")
        click.echo("To minimise checkCIF alerts, fill these in")
        click.echo(
            "You can remove or add cif parameters as long as correct CIF syntax is used!\n"
        )

        fields = yaml_extraction("pipeline-general-extra")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-general-extra")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            full = General_Pipeline()

            full.make_dirs(cfg["experiment_location"])

            full.file_tree_setup(
                cfg["reference_location"],
                cfg["experiment_location"],
                cfg["varying_parameter_values"],
                cfg["varying_cif_parameter"],
            )

            full.process(
                cfg["experiment_location"],
                full.reference_res,
                full.stats_location,
                cfg["refinements_to_check"],
                cfg["tolerance"],
                cfg["maximum_cycles"],
            )

            full.analyse(
                full.reference_res,
                cfg["experiment_location"],
                full.results_location,
                "cx-asap",
                full.chemical_formula,
                full.crystal_habit,
                full.crystal_colour,
                full.max_dimension,
                full.mid_dimension,
                full.min_dimension,
                cfg["structural_analysis_bonds"],
                cfg["structural_analysis_angles"],
                cfg["structural_analysis_torsions"],
                cfg["cif_parameters"],
                cfg["atoms_for_analysis"],
                cfg["varying_cif_parameter"],
                False,
                "instrument.cif",
                True,
            )

            copy_logs(full.results_location)

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


###-----Intensity Compare-----####

"""This pipeline will analyse the intensities of a single .hkl file    
    based on reflection conditions.
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag 
        This will set that parameter as 'TRUE', while the others are 'FALSE'
        This will define the value in the 'if/elif' statements
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
    """

@click.command(
    "module-intensity-compare", short_help="Compares intensities of hkl files"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def module_intensity_compare(dependencies, files, configure, run):
    """This pipeline will analyse the intensities of a single .hkl file

    based on reflection conditions.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- xprep")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a .hkl file")
        click.echo("\nYour .hkl file can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(" - a_axis: enter the a-axis of the unit cell")
        click.echo(" - alpha: enter the alpha angle of the unit cell")
        click.echo(" - b_axis: enter the b-axis of the unit cell")
        click.echo(" - beta: enter the beta angle of the unit cell")
        click.echo(" - c_axis: enter the c-axis of the unit cell")

        click.echo(" - gamma: enter the gamma angle of the unit cell")

        click.echo(
            " - h_condition_1: condition for group 1 of h-indicies, enter n for no condition"
        )
        click.echo(
            " - h_condition_2: condition for group 2 of h-indicies, enter n for no condition"
        )
        click.echo(
            " - k_condition_1: condition for group 1 of k-indicies, enter n for no condition"
        )
        click.echo(
            " - k_condition_2: condition for group 2 of k-indicies, enter n for no condition"
        )
        click.echo(
            " - l_condition_1: condition for group 1 of l-indicies, enter n for no condition"
        )
        click.echo(
            " - l_condition_2: condition for group 2 of l-indicies, enter n for no condition"
        )
        click.echo(
            " - h+k_condition_1: condition for group 1 of h+k indicies, enter n for no condition"
        )
        click.echo(
            " - h+k_condition_2: condition for group 2 of h+k indicies, enter n for no condition"
        )
        click.echo(
            " - h+l_condition_1: condition for group 1 of h+l indicies, enter n for no condition"
        )
        click.echo(
            " - h+l_condition_2: condition for group 2 of h+l indicies, enter n for no condition"
        )
        click.echo(
            " - k+l_condition_1: condition for group 1 of k+l indicies, enter n for no condition"
        )
        click.echo(
            " - k+l_condition_2: condition for group 2 of k+l indicies, enter n for no condition"
        )
        click.echo(
            " - h+k+l_condition_1: condition for group 1 of h+k+l indicies, enter n for no condition"
        )
        click.echo(
            " - h+k+l_condition_2: condition for group 2 of h+k+l indicies, enter n for no condition"
        )

        click.echo(
            " - include_h_0_condition_1: enter true if you want to include reflections where h = 0 for group 1, otherwise enter false"
        )
        click.echo(
            " - include_h_0_condition_2: enter true if you want to include reflections where h = 0 for group 2, otherwise enter false"
        )
        click.echo(
            " - include_k_0_condition_1: enter true if you want to include reflections where k = 0 for group 1, otherwise enter false"
        )
        click.echo(
            " - include_k_0_condition_2: enter true if you want to include reflections where k = 0 for group 2, otherwise enter false"
        )
        click.echo(
            " - include_l_0_condition_1: enter true if you want to include reflections where l = 0 for group 1, otherwise enter false"
        )
        click.echo(
            " - include_l_0_condition_2: enter true if you want to include reflections where l = 0 for group 2, otherwise enter false"
        )
        click.echo(" - xprep_file_name: enter the full path to your hkl file")
        click.echo(
            "\nNote that your conditions must be in the form Xn+Y, where X and Y are integers, and + can also be a -"
        )

        fields = yaml_extraction("module-intensity-compare")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("module-intensity-compare")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            xprep = Intensity_Compare()

            xprep.analyse_reflections(
                cfg["h_condition_1"],
                cfg["k_condition_1"],
                cfg["l_condition_1"],
                cfg["h_condition_2"],
                cfg["k_condition_2"],
                cfg["l_condition_2"],
                cfg["h+k_condition_1"],
                cfg["h+k_condition_2"],
                cfg["h+l_condition_1"],
                cfg["h+l_condition_2"],
                cfg["k+l_condition_1"],
                cfg["k+l_condition_2"],
                cfg["h+k+l_condition_1"],
                cfg["h+k+l_condition_2"],
                cfg["xprep_file_name"],
                cfg["a_axis"],
                cfg["b_axis"],
                cfg["c_axis"],
                cfg["alpha"],
                cfg["beta"],
                cfg["gamma"],
                cfg["include_h_0_condition_1"],
                cfg["include_k_0_condition_1"],
                cfg["include_l_0_condition_1"],
                cfg["include_h_0_condition_2"],
                cfg["include_k_0_condition_2"],
                cfg["include_l_0_condition_2"],
            )

            copy_logs(pathlib.Path(cfg["xprep_file_name"]).parent)

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


###-----Intensity Compare-----####

"""This pipeline will investigate a series of .hkl files and compare the   
    intensities based on different reflection conditions.
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'
        This will define the value in the 'if/elif' statements
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
    """


@click.command(
    "pipeline-intensity-compare",
    short_help="Compares intensities of multiple hkl files",
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_intensity_compare(dependencies, files, configure, run):
    """This pipeline will investigate a series of .hkl files and compare the

    intensities based on different reflection conditions.
    """

    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- xprep")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a series of .hkl files in a single folder")
        click.echo("\nThis folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(" - a_axis: enter the a-axis of the unit cell")
        click.echo(" - alpha: enter the alpha angle of the unit cell")
        click.echo(" - b_axis: enter the b-axis of the unit cell")
        click.echo(" - beta: enter the beta angle of the unit cell")
        click.echo(" - c_axis: enter the c-axis of the unit cell")
        click.echo(" - data_location: enter the full path to your folder of hkl files")

        click.echo(" - gamma: enter the gamma angle of the unit cell")

        click.echo(
            " - h_condition_1: condition for group 1 of h-indicies, enter n for no condition"
        )
        click.echo(
            " - h_condition_2: condition for group 2 of h-indicies, enter n for no condition"
        )
        click.echo(
            " - k_condition_1: condition for group 1 of k-indicies, enter n for no condition"
        )
        click.echo(
            " - k_condition_2: condition for group 2 of k-indicies, enter n for no condition"
        )
        click.echo(
            " - l_condition_1: condition for group 1 of l-indicies, enter n for no condition"
        )
        click.echo(
            " - l_condition_2: condition for group 2 of l-indicies, enter n for no condition"
        )
        click.echo(
            " - h+k_condition_1: condition for group 1 of h+k indicies, enter n for no condition"
        )
        click.echo(
            " - h+k_condition_2: condition for group 2 of h+k indicies, enter n for no condition"
        )
        click.echo(
            " - h+l_condition_1: condition for group 1 of h+l indicies, enter n for no condition"
        )
        click.echo(
            " - h+l_condition_2: condition for group 2 of h+l indicies, enter n for no condition"
        )
        click.echo(
            " - k+l_condition_1: condition for group 1 of k+l indicies, enter n for no condition"
        )
        click.echo(
            " - k+l_condition_2: condition for group 2 of k+l indicies, enter n for no condition"
        )
        click.echo(
            " - h+k+l_condition_1: condition for group 1 of h+k+l indicies, enter n for no condition"
        )
        click.echo(
            " - h+k+l_condition_2: condition for group 2 of h+k+l indicies, enter n for no condition"
        )

        click.echo(
            " - include_h_0_condition_1: enter true if you want to include reflections where h = 0 for group 1, otherwise enter false"
        )
        click.echo(
            " - include_h_0_condition_2: enter true if you want to include reflections where h = 0 for group 2, otherwise enter false"
        )
        click.echo(
            " - include_k_0_condition_1: enter true if you want to include reflections where k = 0 for group 1, otherwise enter false"
        )
        click.echo(
            " - include_k_0_condition_2: enter true if you want to include reflections where k = 0 for group 2, otherwise enter false"
        )
        click.echo(
            " - include_l_0_condition_1: enter true if you want to include reflections where l = 0 for group 1, otherwise enter false"
        )
        click.echo(
            " - include_l_0_condition_2: enter true if you want to include reflections where l = 0 for group 2, otherwise enter false"
        )
        click.echo(
            "\nNote that your conditions must be in the form Xn+Y, where X and Y are integers, and + can also be a -"
        )

        fields = yaml_extraction("pipeline-intensity-compare")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")
        check, cfg = configuration_check("pipeline-intensity-compare")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            xprep = Intensity_Pipeline()

            xprep.multiple_intensity(
                cfg["data_location"],
                cfg["h_condition_1"],
                cfg["k_condition_1"],
                cfg["l_condition_1"],
                cfg["h_condition_2"],
                cfg["k_condition_2"],
                cfg["l_condition_2"],
                cfg["h+k_condition_1"],
                cfg["h+k_condition_2"],
                cfg["h+l_condition_1"],
                cfg["h+l_condition_2"],
                cfg["k+l_condition_1"],
                cfg["k+l_condition_2"],
                cfg["h+k+l_condition_1"],
                cfg["h+k+l_condition_2"],
                cfg["a_axis"],
                cfg["b_axis"],
                cfg["c_axis"],
                cfg["alpha"],
                cfg["beta"],
                cfg["gamma"],
                cfg["include_h_0_condition_1"],
                cfg["include_k_0_condition_1"],
                cfg["include_l_0_condition_1"],
                cfg["include_h_0_condition_2"],
                cfg["include_k_0_condition_2"],
                cfg["include_l_0_condition_2"],
            )

            copy_logs(cfg["data_location"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


###------Cif Merge Module-------###

"""This module will merge one instrument CIF with one structure cif.
    This is useful if it was not done automatically through your structure    
    solution/refinement software.
    CheckCIF will also be run to quickly validate the output.
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 

    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command("module-cif-merge", short_help="Merge two CIFs")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def module_cif_merge(dependencies, files, configure, run):
    """This module will merge one instrument CIF with one structure cif.
    This is useful if it was not done automatically through your structure
    solution/refinement software.
    CheckCIF will also be run to quickly validate the output.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- platon")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - instrument cif")
        click.echo(" - structure cif")
        click.echo("\nThese files can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - instrument_cif: enter the full path to your CIF file containing instrument data"
        )
        click.echo(
            " - new_cif: enter the full path to your CIF file containing a structure you want to add instrument data into"
        )

        fields = yaml_extraction("module-cif-merge")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("module-cif-merge")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            merge = Cif_Merge()
            merge.import_CIFs(
                pathlib.Path(cfg["instrument_cif"]), pathlib.Path(cfg["new_cif"])
            )
            merge.merge_CIFs()
            merge.write_out(
                pathlib.Path(cfg["new_cif"]).parent,
                "combined.cif",
                "check_CIF.chk",
                pathlib.Path(cfg["new_cif"]).name,
            )
            merge.validate_CIFs("combined.cif")

            copy_logs(pathlib.Path(cfg["new_cif"]).parent)

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


#####------Make Cif Module----######

"""This module will create an instrument CIF from a fully completed structural CIF file.
    It might be useful. It might not be. The option is there for you to use.   
    Most of these click functions are specifying text output to commandline     
    The main coding functions are checking the input of the yaml file    
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'       
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command("module-make-instrument-cif", short_help="Generate an instrument cif")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def module_instrument_cif_generation(dependencies, files, configure, run):
    """This module will create an instrument CIF from a fully completed structural CIF file.
    It might be useful. It might not be. The option is there for you to use.
    """
    if dependencies:
        click.echo("\nYou do not require any additional software for this module.")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - fully completed CIF file")
        click.echo("\nThis file can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - reference_cif: full path to the completed cif file you want to extract the instrument parameters out of"
        )

        fields = yaml_extraction("module-make-instrument-cif")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("module-make-instrument-cif")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            instrument_cif = Instrument_CIF()
            instrument_cif.read_reference_cif(cfg["reference_cif"])
            instrument_cif.make_instrument_cif()

            copy_logs(pathlib.Path(cfg["reference_cif"]))

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


#####------Cif Combine Pipeline------#####

"""This pipeline will combine a series of CIFs into one file.
    No edits will be made, but a checkCIF report will be given.   
    Most of these click functions are specifying text output to commandline    
    The main coding functions are checking the input of the yaml file    
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command("pipeline-cif-combine", short_help="Combines CIFs")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_cif_combine(dependencies, files, configure, run):
    """This pipeline will combine a series of CIFs into one file.
    No edits will be made, but a checkCIF report will be given.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- platon")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a series of CIF files in a single folder")
        click.echo("\nThis folder can be located anywhere")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - location_of_cifs: the full path to the folder containing all CIF files"
        )

        fields = yaml_extraction("pipeline-cif-combine")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-cif-combine")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            combine = Cif_Combine()
            combine.combine_cifs_single_folder(cfg["location_of_cifs"])

            copy_logs(cfg["location_of_cifs"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


####------Full CIF Pipeline------#####

"""This pipeline will merge a series of CIF files with instrument CIFs.
    All CIFs will then be combined and run through checkCIF.   
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag       
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command("pipeline-cif", short_help="Edits and compiles CIFs")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_cif(dependencies, files, configure, run):
    """This pipeline will merge a series of CIF files with instrument CIFs.
    All CIFs will then be combined and run through checkCIF.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- platon")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(
            " - a series of folders, each containing a new CIF and an instrument CIF"
        )
        click.echo("\nThis series of folders should be in one parent folder")
        click.echo("This parent folder can be located anywhere ")
        click.echo(
            "Your instrument CIFs will either need the same name, or a common file ending that is not just .cif. For example, Rigaku instruments have .cif_od as their instrument CIF"
        )
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(" - chemical_formula: enter the chemical formula")
        click.echo(" - crystal_habit: enter the habit of your crystal")
        click.echo(" - crystal_colour: describe the colour of your crystal")
        click.echo(
            " - experiment_location: full path to the folder containing all folders with CIFs"
        )
        click.echo(
            " - instrument_ending: if your instrument CIFs have different names but a common ending, list it here (ie .cif_od). If not, enter False"
        )
        click.echo(
            " - instrument_file: if your instrument CIFs have identical names and endings, enter it here. If not, enter False"
        )
        click.echo(" - max_crystal_dimension: largest dimension of your crystal")
        click.echo(" - middle_crystal_dimension: middle dimension of your crystal")
        click.echo(" - min_crystal_dimension: smallest dimension of your crystal")
        click.echo(" - structure_solution: list the structure solution software used")

        click.echo(
            "\nNote that one of instrument_ending and instrument_file must be false. If neither apply to your dataset, then you cannot use this pipeline"
        )

        fields = yaml_extraction("pipeline-cif")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-cif")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            cifs = CIF_Compile_Pipeline()
            cifs.configure(
                cfg["experiment_location"],
                cfg["structure_solution"],
                cfg["chemical_formula"],
                cfg["crystal_habit"],
                cfg["crystal_colour"],
                cfg["max_crystal_dimension"],
                cfg["middle_crystal_dimension"],
                cfg["min_crystal_dimension"],
                cfg["instrument_ending"],
                cfg["instrument_file"],
            )
            cifs.compile_cifs(cfg["experiment_location"])

            copy_logs(cfg["experiment_location"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


######----- Rigaku VT Pipeline------####

"""This pipeline will perform a fully automated analysis for variable temperature 
    collections performed on a rigaku diffractometer.
    The only additional file you will need is your reference structure.
    Completed CIF files and a simple analysis will be output.   
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag 
        This will set that parameter as 'TRUE', while the others are 'FALSE'
        This will define the value in the 'if/elif' statements
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""


@click.command(
    "pipeline-rigaku-vt", short_help="Full pipeline for rigaku VT collections"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_rigaku_vt(dependencies, files, configure, run):
    """This pipeline will perform a fully automated analysis for variable temperature
    collections performed on a rigaku diffractometer.
    The only additional file you will need is your reference structure.
    Completed CIF files and a simple analysis will be output.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- SHELXL")
        click.echo("- PLATON")

    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a reference .ins or .res file")
        click.echo(" - a VT dataset from a rigaku diffractometer")
        click.echo(
            "\nYour VT dataset is likely named after the radiation source you used (ie Mo for molybdenum)"
        )
        click.echo("Your VT dataset can be located anywhere.")
        click.echo(
            "Your reference files should be located outside of the VT dataset folder."
        )
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - ADP_analysis: enter True for adp analysis, otherwise enter False"
        )
        click.echo(
            " - atoms_for_analysis: enter the atom labels for graphical structural analysis as a list (best suited to small numbers to avoid over cluttering graphs"
        )
        click.echo(" - chemical_formula: enter the chemical formula of your crystal")
        click.echo(
            " - cif_parameters: enter the parameters to be extracted from your cifs - the default parameters are usually fine"
        )
        click.echo(" - crystal_colour: enter the colour of your crystal")
        click.echo(" - crystal_habit: enter the habit of your crystal")
        click.echo(
            " - experiment_location: enter the full path to the folder containing your VT data set"
        )
        click.echo(" - maximum_cycles: enter the max number of cycles shelxl can run")
        click.echo(
            " - max_crystal_dimension: enter the largest dimension of your crystal (in mm)"
        )
        click.echo(
            " - middle_crystal_dimension: enter the middle diimension of your crystal (in mm)"
        )
        click.echo(
            " - min_crystal_dimension: enter the smallest dimension of your crystal (in mm)"
        )
        click.echo(
            " - reference_location: enter the full path to your reference .ins/.res file"
        )
        click.echo(
            " - refinements_to_check: enter the number of refinements you want to check for convergence"
        )
        click.echo(
            " - structural_analysis_bonds: enter True for bond length analysis, otherwise enter False"
        )
        click.echo(
            " - structural_analysis_angles: enter True for angle analysis, otherwise enter False"
        )
        click.echo(
            " - structural_analysis_torsions: enter True for torsion analysis, otherwise enter False - note that this will have required the CONF command in your reference .ins/.res file"
        )
        click.echo(
            " - tolerance: enter the desired mean shift for the number of refinements in refinements_to_check"
        )

        fields = yaml_extraction("pipeline-rigaku-vt")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-rigaku-vt")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            full_VT = Synergy_VT()

            full_VT.initialise(cfg["experiment_location"])

            full_VT.process(
                cfg["experiment_location"],
                cfg["reference_location"],
                full_VT.stats_location,
                cfg["refinements_to_check"],
                cfg["tolerance"],
                cfg["maximum_cycles"],
            )

            full_VT.analyse(
                cfg["reference_location"],
                cfg["experiment_location"],
                full_VT.results_location,
                "cx-asap",
                cfg["chemical_formula"],
                cfg["crystal_habit"],
                cfg["crystal_colour"],
                cfg["max_crystal_dimension"],
                cfg["middle_crystal_dimension"],
                cfg["min_crystal_dimension"],
                cfg["structural_analysis_bonds"],
                cfg["structural_analysis_angles"],
                cfg["structural_analysis_torsions"],
                cfg["cif_parameters"],
                cfg["atoms_for_analysis"],
                cfg["ADP_analysis"],
                ".cif_od",
                False,
            )

            copy_logs(full_VT.results_location)

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


#####-----Synchrotron VT Pipeline-----#####

"""This pipeline will perform a fully automated analysis for variable temperature 
    collections performed at the Australian Synchrotron MX1/MX2 beamlines.
    It does require successful autoprocessing of the majority of your datasets.
    The only additional file you will need is your reference structure.
    Completed CIF files and a simple analysis will be output.    
    Most of these click functions are specifying text output to commandline    
    The main coding functions are checking the input of the yaml file    
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'
        This will define the value in the 'if/elif' statements
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command(
    "pipeline-aus-synch-vt", short_help="full pipeline for synchrotron vt data"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_aus_synch_vt(dependencies, files, configure, run):
    """This pipeline will perform a fully automated analysis for variable temperature
    collections performed at the Australian Synchrotron MX1/MX2 beamlines.
    It does require successful autoprocessing of the majority of your datasets.
    The only additional file you will need is your reference structure.
    Completed CIF files and a simple analysis will be output.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- xprep")
        click.echo("- shelxl")
        click.echo("- platon")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a reference .ins or .res file")
        click.echo(
            " - all of the autoprocess folders in a single folder with a common root name"
        )
        click.echo("\nThis folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - ADP_analysis: enter True for ADP analysis, otherwise enter False"
        )
        click.echo(
            " - atoms_for_analysis: enter the atom labels for graphical structural analysis as a list (best suited to small numbers to avoid over cluttering graphs"
        )
        click.echo(
            " - atoms_for_rotation_analysis: enter the labels of the atoms for mean plane analysis"
        )
        click.echo(" - chemical_formula: enter the chemical formula of your crystal")
        click.echo(
            " - cif_parameters: enter the parameters to be extracted from your cifs - the default parameters are usually fine"
        )
        click.echo(" - crystal_colour: enter the colour of your crystal")
        click.echo(" - crystal_habit: enter the habit of your crystal")
        click.echo(
            " - experiment_name: root name that is common between all autoprocess folders"
        )
        click.echo(
            " - location_of_autoprocess_folders: full path to the folder containing all autoprocess folders"
        )
        click.echo(
            " - max_crystal_dimension: enter the largest dimension of your crystal (in mm)"
        )
        click.echo(" - maximum_cycles: enter the max number of cycles shelxl can run")
        click.echo(
            " - middle_crystal_dimension: enter the middle diimension of your crystal (in mm)"
        )
        click.echo(
            " - min_crystal_dimension: enter the smallest dimension of your crystal (in mm)"
        )
        click.echo(
            " - reference_location: enter the full path to your reference .ins/.res file"
        )

        click.echo(
            " - reference_plane: enter the reference plane for mean plane analysis as a list"
        )

        click.echo(
            " - refinements_to_check: enter the number of refinements you want to check for convergence"
        )
        click.echo(
            " - structural_analysis_bonds: enter True for bond length analysis, otherwise enter False"
        )
        click.echo(
            " - structural_analysis_angles: enter True for angle analysis, otherwise enter False"
        )
        click.echo(
            " - structural_analysis_torsions: enter True for torsion analysis, otherwise enter False - note that this will have required the CONF command in your reference .ins/.res file"
        )
        click.echo(
            " - tolerance: enter the desired mean shift for the number of refinements in refinements_to_check"
        )
        click.echo(
            " - transformation_matrix: if your space group is incompatable with the matrix '1 0 0 0 1 0 0 0 1' - ie you have higer symmetry to apply, enter matrix here. Otherwise, use default. You can test this by running xprep manually for your structure and seeing what the applied matrix is"
        )

        fields = yaml_extraction("pipeline-aus-synch-vt")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-aus-synch-vt")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            full_vt_analysis = VT_Pipeline()

            full_vt_analysis.setup(
                cfg["location_of_autoprocess_folders"],
                cfg["experiment_name"],
                "variable_temperature_analysis",
                cfg["reference_location"],
                cfg["atoms_for_rotation_analysis"],
            )

            full_vt_analysis.process(
                cfg["transformation_matrix"],
                full_vt_analysis.sys["analysis_path"],
                "XDS_ASCII.HKL_p1",
                full_vt_analysis.sys["space_group"],
                cfg["chemical_formula"],
                full_vt_analysis.sys["ref_path_organised"],
                full_vt_analysis.sys["current_results_path"],
                cfg["refinements_to_check"],
                cfg["tolerance"],
                cfg["maximum_cycles"],
                cfg["reference_plane"],
            )

            ref_cell = [
                full_vt_analysis.sys["ref_a"],
                full_vt_analysis.sys["ref_b"],
                full_vt_analysis.sys["ref_c"],
                full_vt_analysis.sys["ref_alpha"],
                full_vt_analysis.sys["ref_beta"],
                full_vt_analysis.sys["ref_gamma"],
                full_vt_analysis.sys["ref_volume"],
            ]

            full_vt_analysis.analyse(
                ref_cell,
                full_vt_analysis.sys["ref_path_organised"],
                full_vt_analysis.sys["analysis_path"],
                full_vt_analysis.sys["current_results_path"],
                "known structure",
                cfg["chemical_formula"],
                cfg["crystal_habit"],
                cfg["crystal_colour"],
                cfg["max_crystal_dimension"],
                cfg["middle_crystal_dimension"],
                cfg["min_crystal_dimension"],
                True,
                True,
                False,
                cfg["cif_parameters"],
                cfg["atoms_for_analysis"],
                cfg["ADP_analysis"],
                False,
                "autoprocess.cif",
            )

            copy_logs(full_vt_analysis.sys["current_results_path"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


#####-----Module - cell transformation XDS------######

"""This module will calculate the transformation matrix between two 
    XDS_ASCII reflection files.    
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command(
    "module-xds-cell-transform", short_help="calculate the transformation matrix"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def module_xds_cell_transformation(dependencies, files, configure, run):
    """This module will calculate the transformation matrix between two
    XDS_ASCII reflection files.
    """
    if dependencies:
        click.echo("\nYou do not require any additional software in your path!\n")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - two XDS_ASCII reflection files")
        click.echo("\nThese files can be located anywhere")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(" - XDS_ASCII_File_1: full path to the first XDS_ASCII file")
        click.echo(" - XDS_ASCII_File_2: full path to the second XDS_ASCII file")

        fields = yaml_extraction("module-xds-cell-transform")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("module-xds-cell-transform")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            cell_transform = XDS_Cell_Transformation()
            matrix = cell_transform.find_transformation(
                cfg["XDS_ASCII_File_1"], cfg["XDS_ASCII_File_2"]
            )
            print("The transformation matrix is: ")
            print(matrix)

            copy_logs(pathlib.Path(cfg["XDS_ASCII_File_1"]).parent)

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


######----Module XDS-----######

"""This module will reprocess a single dataset using XDS for  
    frames in a .h5 or .img format.
    If you wish to reprocess a series of datasets automatically, 
    you will need the xds reprocessing pipeline.    
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file    
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command(
    "module-xds-reprocess", short_help="reprocess a single dataset using XDS"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def module_xds_reprocess(dependencies, files, configure, run):
    """This module will reprocess a single dataset using XDS for
    frames in a .h5 or .img format.
    If you wish to reprocess a series of datasets automatically,
    you will need the xds reprocessing pipeline.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- XDS")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a folder which contains another folder in it called 'img'")
        click.echo(" - all of the required frames in the img folder")
        click.echo(" - a template XDS.INP file")
        click.echo("\nThe parent folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - experiment_location: full path to the location of your folder containing the img folder"
        )
        click.echo(" - XDS_INP_path: full path to your XDS.INP file")
        click.echo(" - xds_template_name: name of your dataset")

        fields = yaml_extraction("module-xds-reprocess")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("module-xds-reprocess")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            xds = XDS_reprocess()
            xds.run_xds(
                cfg["xds_template_name"],
                cfg["XDS_INP_path"],
                cfg["experiment_location"],
                1,
            )

            copy_logs(cfg["experiment_location"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


######------Module XPREP------#####

"""This module will use xprep to automatically process a single XDS_ASCII dataset.
    It can also perform a matrix transformation.
    To automatically process a series of datasets use the xprep pipeline
    Most of these click functions are specifying text output to commandline
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command("module-xprep", short_help="process a single dataset using xprep")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def module_xprep(dependencies, files, configure, run):
    """This module will use xprep to automatically process a single XDS_ASCII dataset.
    It can also perform a matrix transformation.
    To automatically process a series of datasets use the xprep pipeline
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- xprep")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - an XDS_ASCII reflection file")
        click.echo("\nThis file can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(" - chemical formula: enter the chemical formula of your crystal")
        click.echo(
            " - experiment_location: enter the full path to the folder containing your XDS_ASCII file"
        )
        click.echo(" - space group: enter the space group in xprep-readable form")
        click.echo(
            " - transformation matrix: enter the matrix for transformation (use 1 0 0 0 1 0 0 0 1 for no transformation)"
        )
        click.echo(" - xprep_file_name: name of your XDS_ASCII file")

        fields = yaml_extraction("module-xprep")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("module-xprep")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            xprep = XPREP()
            xprep.run(
                cfg["transformation_matrix"],
                cfg["experiment_location"],
                cfg["xprep_file_name"],
                cfg["space_group"],
                cfg["chemical_formula"],
                "shelxl",
            )

            copy_logs(cfg["experiment_location"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


#####-----Pipeline XDS------#####

"""This module will reprocess a series of datasets using XDS for   
    frames in a .h5 or .img format.    
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command(
    "pipeline-xds-reprocess", short_help="reprocess multiple datasets using XDS"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_xds_reprocess(dependencies, files, configure, run):
    """This module will reprocess a series of datasets using XDS for
    frames in a .h5 or .img format.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- XDS")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a template XDS.INP file")
        click.echo(
            " - A series of folders each with a 'img' folder linking to the frames"
        )
        click.echo(
            " - These folders should be named the same as the dataset they correspond to"
        )
        click.echo(" - These folders should all be located in the same place")
        click.echo("\nThe parent folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - experiment_location: full path to the location of your series of folders containing img folders"
        )
        click.echo(" - XDS_INP_path: full path to your XDS.INP file")

        fields = yaml_extraction("pipeline-xds-reprocess")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-xds-reprocess")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            XDS = XDS_Pipeline()
            XDS.multiple_xds(cfg["experiment_location"], cfg["XDS_INP_path"])

            copy_logs(cfg["experiment_location"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


######------Pipeline XPREP------######

"""This module will use xprep to automatically process a single XDS_ASCII dataset.
    It can also perform a matrix transformation.   
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command("pipeline-xprep", short_help="process multiple datasets using xprep")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_xprep(dependencies, files, configure, run):
    """This module will use xprep to automatically process a single XDS_ASCII dataset.
    It can also perform a matrix transformation.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- xprep")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(
            " - a series of folders each containing an XDS_ASCII file (all have the same name)"
        )
        click.echo(" - these folders should be in the same parent location")
        click.echo("\nThis parent folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(" - chemical formula: enter the chemical formula of your crystal")
        click.echo(
            " - experiment_location: enter the full path to the folder containing all dataset folders"
        )
        click.echo(" - space group: enter the space group in xprep-readable form")
        click.echo(
            " - transformation matrix: enter the matrix for transformation (use 1 0 0 0 1 0 0 0 1 for no transformation)"
        )

        click.echo(" - xprep_file_name: name of your XDS_ASCII file")

        fields = yaml_extraction("pipeline-xprep")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-xprep")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            xprep = XPREP_Pipeline()
            xprep.multiple_xprep(
                cfg["transformation_matrix"],
                cfg["experiment_location"],
                cfg["xprep_file_name"],
                cfg["space_group"],
                cfg["chemical_formula"],
            )

            copy_logs(cfg["experiment_location"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


#####----- Module Cell Analysis------#####

"""This module will analyse changes in the unit cell from a .csv file.
    If you have already run the cif reading module, then the file should be in the correct format.  
    Most of these click functions are specifying text output to commandline     
    The main coding functions are checking the input of the yaml file    
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command("module-cell-analysis", short_help="analyse changes in unit cells")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def module_cell_analysis(dependencies, files, configure, run):
    """This module will analyse changes in the unit cell from a .csv file.
    If you have already run the cif reading module, then the file should be in the correct format.
    """
    if dependencies:
        click.echo("\nYou do not require any additional software in your path!\n")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a .csv file with full unit cell information")
        click.echo(
            " - the headings of your unit cell columns must be in proper .cif format"
        )
        click.echo(" For example: a-axis is written as '_cell_length_a'")
        click.echo("\nThis file can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(" - csv_location: full path to your .csv file")
        click.echo(
            " - reference_unit_cell: enter the neutral unit cell for datasets to be compared to"
        )
        click.echo(
            " - x_axis_header: enter the column title for your independent variable (ie temperature or pressure"
        )

        fields = yaml_extraction("module-cell-analysis")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("module-cell-analysis")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            analysis = Cell_Deformation()
            analysis.import_data(cfg["csv_location"], cfg["reference_unit_cell"])
            analysis.calculate_deformations()
            analysis.graphical_analysis(cfg["x_axis_header"], cfg["x_axis_header"])

            # This is defined outside of the class, because when doing multiple analysis for a series of mapping experiments, the dictionary will be different, however, as this is a standalone module and not a pipeline for multiple structures, it is defined here

            try:
                test1 = analysis.df["_diffrn_reflns_av_R_equivalents"]
                test2 = analysis.df["_refine_ls_R_factor_gt"]
                test3 = analysis.df["_diffrn_measured_fraction_theta_full"]
            except:
                # self.logger.critical('No data quality statistics found in imported .csv file')
                logging.critical(
                    __name__
                    + " : No data quality statistics found in imported .csv file"
                )
                print("Error! Check logs")
                exit()

            y_dict = {
                "R1": analysis.df["_refine_ls_R_factor_gt"],
                "Rint": analysis.df["_diffrn_reflns_av_R_equivalents"],
                "Completeness": analysis.df["_diffrn_measured_fraction_theta_full"],
            }

            analysis.quality_analysis(
                cfg["x_axis_header"], y_dict, analysis.cfg["x_axis_header"]
            )

            copy_logs(pathlib.Path(cfg["csv_location"]).parent)

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


######-------Module Cif Read--------######

"""This module will extract data from a series of CIF files and    
    output them into a .csv file.
    Examples of parameters for output include (but are not limited to) 
    unit cell parameters and quality statistics.    
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file    
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command("module-cif-read", short_help="extract parameters from cif files")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def module_cif_read(dependencies, files, configure, run):
    """This module will extract data from a series of CIF files and
    output them into a .csv file.
    Examples of parameters for output include (but are not limited to)
    unit cell parameters and quality statistics.
    """
    if dependencies:
        click.echo("\nYou do not require any additional software in your path!\n")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a series of .cif files located in a single folder")
        click.echo("\nThis folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - ADP_analysis: enter 'true' if you want to extract ADP information, otherwise enter 'false'"
        )
        click.echo(
            " - cif_parameters: these are the parameters that will be extracted from the cif - default ones are usually enough - note that any additional ones must be written in exact cif format"
        )
        click.echo(
            " - structural_analysis_bonds: enter 'true' if you want to extract bond information, otherwise enter 'false'"
        )
        click.echo(
            " - structural_analysis_angles: enter 'true' if you want to extract angle information, otherwise enter 'false'"
        )
        click.echo(
            " - structural_analysis_torsions: enter 'true' if you want to extract torsion information, otherwise enter 'false' - note that cif files will only contain this information if you refined your structures with the 'CONF' command"
        )

        fields = yaml_extraction("module-cif-read")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("module-cif-read")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            analysis = CIF_Read()
            analysis.configure(cfg["cif_parameters"])
            analysis.get_data(
                cfg["folder_containing_cifs"],
                cfg["structural_analysis_bonds"],
                cfg["structural_analysis_angles"],
                cfg["structural_analysis_torsions"],
                cfg["ADP_analysis"],
            )
            analysis.data_output()

            copy_logs(cfg["folder_containing_cifs"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


######----- Module Rotation Planes ------#####

"""For a single dataset refined with the MPLA command, the angle to a reference plane can be calculated.
    To automatically analyse a series of files, use the rotation pipeline.  
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file    
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command(
    "module-rotation-planes", short_help="calculate MPLA angle to reference plane"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def module_rotation_planes(dependencies, files, configure, run):
    """For a single dataset refined with the MPLA command, the angle to a reference plane can be calculated.
    To automatically analyse a series of files, use the rotation pipeline.
    """
    if dependencies:
        click.echo("\nYou do not require any additional software in your path!\n")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(
            " - a .lst file output after refinement in SHELXL with the MPLA command"
        )
        click.echo("\nThis file can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - lst_file_location: enter the full path to your lst file for analysis"
        )
        click.echo(
            " - reference_plane: fill out the list with the three numbers that form your reference crystallographic plane. For example, to compare to the (100) plane, enter the three numbers '1', '0', and '0' in the three positions."
        )

        fields = yaml_extraction("module-rotation-planes")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("module-rotation-planes")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            rotation_analysis = Rotation()
            rotation_analysis.configure(cfg["reference_plane"])
            rotation_analysis.analysis(
                cfg["lst_file_location"],
                1,
                pathlib.Path(cfg["lst_file_location"]).parent,
            )

            copy_logs(pathlib.Path(cfg["lst_file_location"]).parent)

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


#######------- Module Structure -------#######

"""This module will examine structural data in .csv files and separate/graph the important parts of interest.
    For example, if you are investigating the spin-crossover properties of a metal complex, you might be most interested in the M-L bond lengths.
    By specifying the metal atom, you graphs will be provided and a separated spreadsheet of those bonds will be provided.
    This module currently requires the use of the 'cif_read' module prior.   
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file    
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command(
    "module-structural-analysis", short_help="analyse bond/angle/torsion data"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def module_structural_analysis(dependencies, files, configure, run):
    """This module will examine structural data in .csv files and separate/graph the important parts of interest.
    For example, if you are investigating the spin-crossover properties of a metal complex, you might be most interested in the M-L bond lengths.
    By specifying the metal atom, you graphs will be provided and a separated spreadsheet of those bonds will be provided.
    This module currently requires the use of the 'cif_read' module prior.
    """
    if dependencies:
        click.echo("\nYou do not require any additional software in your path!\n")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a .csv file of your bond data")
        click.echo(" - a .csv file of your angle data")
        click.echo(" - a .csv file of your torsion data")
        click.echo(" Note that do you not need all of these files present")
        click.echo("\nThese files can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - angle_data: enter the full path to your angle .csv file (or false if you do not have one)"
        )
        click.echo(
            " - atoms_for_analysis: enter the label of the atoms you are most interested in"
        )
        click.echo(
            " - bond_data: enter the full path to your bond .csv file (or false if you do not have one)"
        )

        click.echo(
            " - torsion_data: enter the full path to your torsion .csv file (or false if you do not have one)"
        )

        fields = yaml_extraction("module-structural-analysis")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("module-structural-analysis")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            geometry = Structural_Analysis()
            geometry.import_and_analyse(
                cfg["bond_data"],
                cfg["angle_data"],
                cfg["torsion_data"],
                cfg["atoms_for_analysis"],
            )

            copy_logs(pathlib.Path(cfg["bond_data"]).parent)

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")

#####------ Pipeline Rotation Planes-----#######

"""For a series of datasets refined with the MPLA command, the angle    
    to a reference plane can be calculated and compared.    
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file    
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command(
    "pipeline-rotation-planes", short_help="calculate multiple MPLA angles to reference"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_rotation_planes(dependencies, files, configure, run):
    """For a series of datasets refined with the MPLA command, the angle
    to a reference plane can be calculated and compared.
    """
    if dependencies:
        click.echo("\nYou do not require any additional software in your path!\n")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(
            " - a series of .lst files in separate folders contained in a single parent folder"
        )
        click.echo("\nThis parent folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - experiment_location: full path to the parent folder containing a series of folders with .lst files inside"
        )
        click.echo(
            " - reference_plane: fill out the list with the three numbers that form your reference crystallographic plane. For example, to compare to the (100) plane, enter the three numbers '1', '0', and '0' in the three positions."
        )

        fields = yaml_extraction("pipeline-rotation-planes")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-rotation-planes")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            multi_rotation = Rotation_Pipeline()
            multi_rotation.analysis(
                cfg["experiment_location"],
                cfg["reference_plane"],
                cfg["experiment_location"]
            )

            copy_logs(cfg["experiment_location"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


#######------Pipeline varying parameter ---------######

"""This pipeline will analyse .cif files for a dynamic experiment where one    
    parameter in the cif files is changing. For example, this might be a     
    change in temperature, pressure, etc.
    It will output graphs displaying changes in unit cell parameters and defined structural changes.
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file    
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command("pipeline-variable-analysis", short_help="general analysis of cif files")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_variable_analysis(dependencies, files, configure, run):
    """This pipeline will analyse .cif files for a dynamic experiment where one
    parameter in the cif files is changing. For example, this might be a
    change in temperature, pressure, etc.
    It will output graphs displaying changes in unit cell parameters and defined structural changes.
    """
    if dependencies:
        click.echo("\nYou do not require any additional software in your path!\n")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a series of .cif files located in a single folder")
        click.echo("\nThis folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - ADP_analysis: enter 'true' if you want to extract ADP information, otherwise enter 'false'"
        )
        click.echo(
            " - atoms_for_analysis: enter the label of the atoms you are most interested in"
        )
        click.echo(
            " - cif_parameters: these are the parameters that will be extracted from the cif - default ones are usually enough - note that any additional ones must be written in exact cif format"
        )
        click.echo(
            " - experiment_location: enter the full path to the folder which contains your .cif files"
        )
        click.echo(" - reference_unit_cell: enter the path to a reference .ins file")
        click.echo(
            " - structural_analysis_bonds: enter 'true' if you want to extract bond information, otherwise enter 'false'"
        )
        click.echo(
            " - structural_analysis_angles: enter 'true' if you want to extract angle information, otherwise enter 'false'"
        )
        click.echo(
            " - structural_analysis_torsions: enter 'true' if you want to extract torsion information, otherwise enter 'false' - note that cif files will only contain this information if you refined your structures with the 'CONF' command"
        )
        click.echo(
            " - varying_cif_parameter: enter the parameter in your .cif files that is changing. Make sure you use proper .cif syntax"
        )

        fields = yaml_extraction("pipeline-variable-analysis")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-variable-analysis")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            variable_analysis = Variable_Analysis_Pipeline()
            variable_analysis.analyse_data(
                cfg["reference_unit_cell"],
                cfg["experiment_location"],
                cfg["cif_parameters"],
                cfg["atoms_for_analysis"],
                cfg["varying_cif_parameter"],
                cfg["structural_analysis_bonds"],
                cfg["structural_analysis_angles"],
                cfg["structural_analysis_torsions"],
                cfg["ADP_analysis"],
            )

            copy_logs(cfg["experiment_location"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


####-----Pipeline Variable Position Analysis-----####

"""This pipeline will analyse .cif files for a variable    
    position experiment.
    It will output graphs displaying changes in unit cell 
    parameters and defined structural changes.    
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file    
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command(
    "pipeline-position-analysis", short_help="analysis of variable position cifs"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_position_analysis(dependencies, files, configure, run):
    """This pipeline will analyse .cif files for a variable
    position experiment.
    It will output graphs displaying changes in unit cell
    parameters and defined structural changes.
    """
    if dependencies:
        click.echo("\nYou do not require any additional software in your path!\n")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a series of .cif files located in a single folder")
        click.echo("\nThis folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - ADP_analysis: enter 'true' if you want to extract ADP information, otherwise enter 'false'"
        )
        click.echo(
            " - atoms_for_analysis: enter the label of the atoms you are most interested in"
        )
        click.echo(
            " - atoms_for_rotation_analysis: enter the atoms you used in the MPLA command"
        )
        click.echo(
            " - cif_parameters: these are the parameters that will be extracted from the cif - default ones are usually enough - note that any additional ones must be written in exact cif format"
        )
        click.echo(
            " - experiment_location: enter the full path to your folder containing your cif files"
        )
        click.echo(
            " - mapping_step_size: enter the step-size for the collection (in um)"
        )
        click.echo(
            " - min_pixels: enter the values for minimum_pixels_in_a_spot as a list for XDS proessing"
        )
        click.echo(
            " - reference_plane: fill out the list with the three numbers that form your reference crystallographic plane. For example, to compare to the (100) plane, enter the three numbers '1', '0', and '0' in the three positions."
        )
        click.echo(
            " - reference_unit_cell: enter the neutral unit cell for datasets to be compared to"
        )
        click.echo(
            " - sepmin: enter the values for sepmin as a list for XDS processing"
        )
        click.echo(
            " - spot_maximum_centroid: enter the values for spot_maximum_centroid as a list for XDS processing"
        )
        click.echo(
            " - strong_pixels: enter the values for strong_pixels as a list for XDS processing"
        )

        click.echo(
            " - structural_analysis_bonds: enter 'true' if you want to extract bond information, otherwise enter 'false'"
        )
        click.echo(
            " - structural_analysis_angles: enter 'true' if you want to extract angle information, otherwise enter 'false'"
        )
        click.echo(
            " - structural_analysis_torsions: enter 'true' if you want to extract torsion information, otherwise enter 'false' - note that cif files will only contain this information if you refined your structures with the 'CONF' command"
        )

        click.echo(
            " - wedge_angles: enter the wedge angles as a list for XDS processing"
        )

        fields = yaml_extraction("pipeline-position-analysis")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-position-analysis")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            vp_analysis = VP_Analysis_Pipeline()
            vp_analysis.analyse_data(
                cfg["reference_unit_cell"],
                cfg["experiment_location"],
                cfg["cif_parameters"],
                cfg["atoms_for_analysis"],
                cfg["atoms_for_plane"],
                cfg["mapping_step_size"],
                cfg["spot_maximum_centroid"],
                cfg["min_pixels"],
                cfg["strong_pixels"],
                cfg["sepmin"],
                cfg["wedge_angles"],
                cfg["reference_plane"],
                cfg["structural_analysis_bonds"],
                cfg["structural_analysis_angles"],
                cfg["structural_analysis_torsions"],
                cfg["ADP_analysis"],
            )

            copy_logs(cfg["experiment_location"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


#####------Pipeline Variable Temperature Analysis-----####

"""This pipeline will analyse .cif files for a variable    
    temperature experiment.
    It will output graphs displaying changes in unit cell parameters 
    and defined structural changes.    
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file    
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command(
    "pipeline-temperature-analysis", short_help="analysis of variable temperature cifs"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_temperature_analysis(dependencies, files, configure, run):
    """This pipeline will analyse .cif files for a variable
    temperature experiment.
    It will output graphs displaying changes in unit cell parameters
    and defined structural changes.
    """
    if dependencies:
        click.echo("\nYou do not require any additional software in your path!\n")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a series of .cif files located in a single folder")
        click.echo("\nThis folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - ADP_analysis: enter 'true' if you want to extract ADP information, otherwise enter 'false'"
        )
        click.echo(
            " - atoms_for_analysis: enter the label of the atoms you are most interested in"
        )
        click.echo(
            " - cif_parameters: these are the parameters that will be extracted from the cif - default ones are usually enough - note that any additional ones must be written in exact cif format"
        )
        click.echo(
            " - experiment_location: enter the full path to the folder which contains your .cif files"
        )
        click.echo(
            " - reference_unit_cell: enter the neutral unit cell for datasets to be compared to"
        )
        click.echo(
            " - structural_analysis_bonds: enter 'true' if you want to extract bond information, otherwise enter 'false'"
        )
        click.echo(
            " - structural_analysis_angles: enter 'true' if you want to extract angle information, otherwise enter 'false'"
        )
        click.echo(
            " - structural_analysis_torsions: enter 'true' if you want to extract torsion information, otherwise enter 'false' - note that cif files will only contain this information if you refined your structures with the 'CONF' command"
        )

        fields = yaml_extraction("pipeline-temperature-analysis")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-temperature-analysis")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            vt_analysis = VT_Analysis_Pipeline()
            vt_analysis.analyse_data(
                cfg["reference_unit_cell"],
                cfg["experiment_location"],
                cfg["cif_parameters"],
                cfg["atoms_for_analysis"],
                cfg["structural_analysis_bonds"],
                cfg["structural_analysis_angles"],
                cfg["structural_analysis_torsions"],
                cfg["ADP_analysis"],
            )

            copy_logs(cfg["experiment_location"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


### -------- Module Platon Squeeze --------###

"""This module will run squeeze on a single .ins file via PLATON.   
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command("module-platon-squeeze", short_help="Run Platon Squeeze")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def module_platon_squeeze(dependencies, files, configure, run):
    """This module will run squeeze on a single .ins file via PLATON."""
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- PLATON")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - .ins/.hkl files in a single folder")
        click.echo("\nThis folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(" - file_name: full path to your .ins file")

        fields = yaml_extraction("module-platon-squeeze")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("module-platon-squeeze")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            squeeze = Platon_Squeeze()
            squeeze.run_squeeze(cfg["file_name"])

            copy_logs(pathlib.Path(cfg["file_name"]).parent)

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


### ------ Pipeline Platon Squeeze ---- ###

"""This pipeline will run squeeze on a series of .ins files via PLATON.   
    It will do this is a separate directory to retain the original structures 
    for comparison.
    If you wish to use this in a larger pipeline, it is recommended 
    you use pipeline-refinement, followed by this squeeze pipeline,     
    then finally pipeline-general on the new squeezed folder.     
    Make sure you use a squeezed reference for the general pipeline!    
    Most of these click functions are specifying text output to commandline
    The main coding functions are checking the input of the yaml file    
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command(
    "pipeline-platon-squeeze", short_help="Run squeeze over multiple structures"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_platon_squeeze(dependencies, files, configure, run):
    """This pipeline will run squeeze on a series of .ins files via PLATON.
    It will do this is a separate directory to retain the original structures
    for comparison.
    If you wish to use this in a larger pipeline, it is recommended
    you use pipeline-refinement, followed by this squeeze pipeline,
    then finally pipeline-general on the new squeezed folder.
    Make sure you use a squeezed reference for the general pipeline!
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- PLATON")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a series of .hkl files")
        click.echo(" - a series .ins files corresponding to the .hkl files")
        click.echo(
            "\nYour .hkl/.ins files should be in separate folders located in a single parent folder"
        )
        click.echo("\nThis folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - experiment_location: enter the full path to your folder containing all dataset folders"
        )

        fields = yaml_extraction("pipeline-platon-squeeze")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-platon-squeeze")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            squeeze = Squeeze_Pipeline()
            squeeze.new_squeeze_directory(cfg["experiment_location"])
            squeeze.multi_squeeze(squeeze.new_location)

            copy_logs(cfg["experiment_location"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


###------Module ADP Analysis-------###

"""This module will analyse the ADPs that have been previously   
    extracted from CIF files. It is recommended that you run this 
    AFTER module-cif-read OR use a pipeline with this included    
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file    
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""


@click.command("module-ADP-analysis", short_help="analyse ADPs extracted from CIFs")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def module_adp_analysis(dependencies, files, configure, run):
    """This module will analyse the ADPs that have been previously
    extracted from CIF files. It is recommended that you run this
    AFTER module-cif-read OR use a pipeline with this included
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- None :)")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(
            " - A .csv file with your extracted ADP Data (recommended you run module-cif-read first)"
        )
        click.echo(
            " - A .csv file with your extracted unit cell Data (recommended you run module-cif-read first!)"
        )
        click.echo("\nThis file can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - cell_path: Enter the path to the .csv file containing extracted unit cell data (recommended you run module-cif-read first!)"
        )
        click.echo(
            " - csv_path: Enter the path to the .csv file containing extracted ADP data (recommended you run module-cif-read first)"
        )

        fields = yaml_extraction("module-adp-analysis")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("module-adp-analysis")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            adps = ADP_analysis()
            adps.analyse_data(cfg["csv_path"], cfg["cell_path"])

            copy_logs(pathlib.Path(cfg["csv_path"]).parent)

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


### ------ Pipeline AS Brute ---- ###

"""This pipeline will run xprep and then shelxt on a    
    series of XDS_ASCII.HKL_p1 files from the Australian 
    Synchrotron autoprocesing process. It will produce cxasap.ins 
    and cxasap.hkl if xprep is successful which will then be solved by 
    shelxt. A default cell contents is used.
    Structures that fail at any step will need closer attention with a 
    combination of XDS and manual xprep.
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'
        This will define the value in the 'if/elif' statements
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""

@click.command(
    "pipeline-AS-Brute", short_help="Run xprep/shelxt over multiple structures"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_AS_Brute(dependencies, files, configure, run):
    """This pipeline will run xprep and then shelxt on a
    series of XDS_ASCII.HKL_p1 files from the Australian
    Synchrotron autoprocesing process. It will produce cxasap.ins
    and cxasap.hkl if xprep is successful which will then be solved by
    shelxt. A default cell contents is used.
    Structures that fail at any step will need closer attention with a
    combination of XDS and manual xprep.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- XPREP")
        click.echo("- SHELXT")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a folder containing your autoprocessed data ")
        click.echo("\nThis folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - experiment_location: enter the full path to your folder containing all dataset folders"
        )

        fields = yaml_extraction("pipeline-AS-Brute")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-AS-Brute")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            asbrute = AS_Brute()
            asbrute.initialise(cfg["experiment_location"])
            asbrute.xprepreduce(cfg["experiment_location"])
            asbrute.solve(cfg["experiment_location"])
            asbrute.report(cfg["experiment_location"])

            copy_logs(cfg["experiment_location"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


# ------Module Molecule Reconstruction------#

"""This module will reconstruct the crystal structure in new    
    unit cells given a known rate of change of the unit cell axes. 
    Useful to construct theoretical structures to run calculations for 
    comparison.
    The output of this code is a series of theoretical .res files 
    which can be imported into various software and converted to the 
    required file types for computation.
    Note that this module is currently only compatable with molecules    
    where the starting atom lies on a special position.
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""


@click.command(
    "module-molecule-reconstruction", short_help="reconstruct molecules in new cells"
)
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def module_molecule_reconstruction(dependencies, files, configure, run):
    """This module will reconstruct the crystal structure in new
    unit cells given a known rate of change of the unit cell axes.
    Useful to construct theoretical structures to run calculations for
    comparison.
    The output of this code is a series of theoretical .res files
    which can be imported into various software and converted to the
    required file types for computation.
    Note that this module is currently only compatable with molecules
    where the starting atom lies on a special position.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- nothing :) ")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - reference .res file")
        click.echo("\nThis file can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(" - a_gradient: enter the rate of change of the a-axis")
        click.echo(" - alpha_gradient: enter the rate of change of the alpha angle")
        click.echo(" - atom_list: enter each atom as a list on a new line")
        click.echo(" - b_gradient: enter the rate of change of the b-axis")
        click.echo(" - beta_gradient: enter the rate of change of the beta angle")
        click.echo(
            " - bond_list: enter each bond in the molecule as a list on a new line. Example syntax: Cu1-O1"
        )
        click.echo(" - c_gradient: enter the rate of change of the c-axis")
        click.echo(" - gamma_gradient: enter the rate of change of the gamma angle")
        click.echo(
            " - max_position: enter the maximum position (this value is in the same units as x for your gradient"
        )
        click.echo(
            " - min_position: enter the minimum position (this value is in the same units as x for your gradient"
        )
        click.echo(
            " - position_step_size: enter the step size between positions for your theoretical structures"
        )
        click.echo(" - reference_path: enter the path to your reference .res file")
        click.echo(
            " - starting_atom: enter the atom that the reconstruction should start from"
        )
        click.echo(
            " - starting_coordinates: enter the coordinates of the starting atom"
        )

        fields = yaml_extraction("module-molecule-reconstruction")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("module-molecule-reconstruction")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            reconstruction = Molecule_Reconstruction(
                cfg["reference_path"],
                cfg["atom_list"],
                cfg["bond_list"],
                cfg["starting_atom"],
                cfg["starting_coordinates"],
                cfg["a_gradient"],
                cfg["b_gradient"],
                cfg["c_gradient"],
                cfg["alpha_gradient"],
                cfg["beta_gradient"],
                cfg["gamma_gradient"],
                cfg["max_position"],
                cfg["min_position"],
                cfg["position_step_size"],
            )
            reconstruction.read_in_reference()
            reconstruction.find_construction_order()
            reconstruction.find_internal_vectors()
            for item in reconstruction.positions:
                reconstruction.calculate_new_cell(item)
                reconstruction.calculate_fractional_coordinates()
                reconstruction.write_res(item)

            copy_logs(pathlib.Path(cfg["reference_path"]).parent)

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


### ------ Pipeline SHELXT AUTO ---- ###

"""This pipeline will run shelxt on a series of .ins and .hkl files.    
    It will do this is a separate directory to retain the original structures 
    for comparison.
    It will only retain the first (_a) shelxt result.
    If you wish to use this in a larger pipeline, it is recommended 
    you use pipeline-refinement, followed by this shelxt pipeline, 
    then finally pipeline-general on the new squeezed folder. 
    Make sure you use a squeezed reference for the general pipeline!
    Most of these click functions are specifying text output to commandline 
    The main coding functions are checking the input of the yaml file
    and setting up the corresponding class and calling its functions 
    Args:
        The user will enter one of the four arguments as a flag         
        This will set that parameter as 'TRUE', while the others are 'FALSE'        
        This will define the value in the 'if/elif' statements        
        dependencies (bool): will check for dependencies
        files (bool): will show the user what files are required
        configure (bool): will set up the yaml for the user to fill out
        run (bool): will execute the chosen module/pipeline 
"""


@click.command("pipeline-shelxt-auto", short_help="Run shelxt over multiple structures")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_shelxt_auto(dependencies, files, configure, run):
    """This pipeline will run shexlt on a series of .ins and .hkl files.
    It will do this is a separate directory to retain the original structures
    for comparison.
    It will retain only the first result and rename the .res file as a .ins file.
    If you wish to use this in a larger pipeline, it is recommended
    you use pipeline-refinement, followed by this shelxt pipeline,
    then finally pipeline-general on the new shelxt folder.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- SHELXT")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a series of .hkl files")
        click.echo(" - a series .ins files corresponding to the .hkl files")
        click.echo(
            "\nYour .hkl/.ins files should be in separate folders located in a single parent folder"
        )
        click.echo("\nThis folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(
            " - experiment_location: enter the full path to your folder containing all dataset folders"
        )

        fields = yaml_extraction("pipeline-shelxt-auto")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")

        check, cfg = configuration_check("pipeline-shelxt-auto")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            shelxt = SHELXT_Pipeline_auto()
            shelxt.new_shelxt_directory(cfg["experiment_location"])
            shelxt.multi_shelxt(shelxt.new_location)

            copy_logs(cfg["experiment_location"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


###Pipeline SPECIAL FOR STEPH (AS Brute)###
# """This pipeline is an individual version of the AS-Brute pipeline
#    It is designed to function alongside their autoprocessing pipeline
#    Probably no use for users at home
#    Will run xprep, and shelxt on a single dataset and then compile the results
#    Args:
#        The user will enter one of the four arguments as a flag
#        This will set that parameter as 'TRUE', while the others are 'FALSE'
#        This will define the value in the 'if/elif' statements
#        dependencies (bool): will check for dependencies
#        files (bool): will show the user what files are required
#        configure (bool): will set up the yaml for the user to fill out
#        run (bool): will execute the chosen module/pipeline
#    """

@click.command("pipeline-AS-Brute-individual", short_help="for AS use only")
@click.option("--dependencies", is_flag=True, help="view the software dependencies")
@click.option("--files", is_flag=True, help="view the required input files")
@click.option("--configure", is_flag=True, help="generate your conf.yaml file")
@click.option("--run", is_flag=True, help="run the code!")
def pipeline_AS_Brute_individual(dependencies, files, configure, run):
    """FOR AUSTRALIAN SYNCHROTRON USE ONLY
    This pipeline is an individual version of AS-Brute
    It is designed to be integrated into the AS autoprocessing pipeline
    Does not have much applicability to a single user as it only does xprep and shelxt
    on a single dataset.
    """
    if dependencies:
        click.echo("\nYou require the below software in your path:")
        click.echo("- XPREP")
        click.echo("- SHELXT")
    elif files:
        click.echo("\nYou require the below files:")
        click.echo(" - a dataset folder from the Australian Synchrotron")
        click.echo("\nThis folder can be located anywhere ")
    elif configure:
        click.echo("\nWriting a file called conf.yaml in the cx_asap folder...\n")
        click.echo("You will need to fill out the parameters.")
        click.echo("Descriptions are listed below:")
        click.echo(" - chemical_formula: enter the chemical formula of your crystal")
        click.echo(
            " - experiment_location: enter the full path to your folder containing all dataset folders"
        )

        fields = yaml_extraction("pipeline-AS-Brute-individual")
        yaml_creation(fields)

    elif run:
        click.echo("\nChecking to see if experiment configured....\n")
        check, cfg = configuration_check("pipeline-AS-Brute-individual")

        if check == False:
            click.echo("Make sure you fill in the configuration file!")
            click.echo(
                "If you last ran a different code, make sure you reconfigure for the new script!"
            )
            click.echo("Re-run configuration for description of each parameter\n")
        else:
            click.echo("READY TO RUN SCRIPT!\n")
            reset_logs()
            brute = AS_Brute_Single()
            sadabs_folders = ['/sadabs_m', '/sadabs_s', '/sadabs_w']
            print(cfg["experiment_location"])
            print(os.getcwd())
            brute.initialise(cfg["experiment_location"])
            if cfg["chemical_formula"] == 0:
                brute.xprepreduce(cfg["experiment_location"], "cxasap")
            else:
                brute.xprepreduce(cfg["experiment_location"], "cxasap", cfg["chemical_formula"])
            brute.solve(cfg["experiment_location"], "cxasap")
            brute.move_files(cfg["experiment_location"], cfg["experiment_location"], "cxasap")
            brute.report(cfg["experiment_location"], "")
            
            for item in sadabs_folders:
                if cfg["chemical_formula"] == 0:
                    brute.xprepreduce(cfg["experiment_location"] + item, "cxasap_" + item.strip('/'))
                else:
                    brute.xprepreduce(cfg["experiment_location"] + item, "cxasap_" + item.strip('/'), cfg["chemical_formula"])
                brute.solve(cfg["experiment_location"] + item, "cxasap_" + item.strip('/'))
                brute.move_files(cfg["experiment_location"], cfg["experiment_location"] + item, "cxasap_" + item.strip('/'))
                brute.report(cfg["experiment_location"], item.strip('/'))
                
            
            copy_logs(cfg["experiment_location"])

        output_message()

    else:
        click.echo("Please select an option. To view options, add --help")


os_name = platform.system()
BadOS = False
if os_name == "Windows":
    BadOS = True

# Add or remove modules from here that are windows compatable

windows_modules = [
    test,
    errors,
    module_refinement,
    pipeline_refinement,
    pipeline_general,
    pipeline_general_extra,
    module_cif_merge,
    module_instrument_cif_generation,
    pipeline_cif_combine,
    pipeline_cif,
    module_cell_analysis,
    module_cif_read,
    module_rotation_planes,
    module_structural_analysis,
    pipeline_temperature_analysis,
    pipeline_variable_analysis,
    module_adp_analysis,
    module_platon_squeeze,
    pipeline_platon_squeeze,
]

windows_modules_dev = [
    module_intensity_compare,
    pipeline_intensity_compare,
    pipeline_rigaku_vt,
    module_molecule_reconstruction,
    pipeline_shelxt_auto,
]

if BadOS == True:

    # Modules for master branch ###

    for item in windows_modules:
        cli.add_command(item)

    # Modules for dev branch ###

    for item in windows_modules_dev:
        cli.add_command(item)

else:

    ### Modules for master branch ###

    cli.add_command(test)
    cli.add_command(errors)
    cli.add_command(module_refinement)
    cli.add_command(pipeline_refinement)
    cli.add_command(pipeline_general)
    cli.add_command(pipeline_general_extra)
    cli.add_command(module_cif_merge)
    cli.add_command(module_instrument_cif_generation)
    cli.add_command(pipeline_cif_combine)
    cli.add_command(pipeline_cif)
    cli.add_command(module_cell_analysis)
    cli.add_command(module_cif_read)
    cli.add_command(module_rotation_planes)
    cli.add_command(module_structural_analysis)
    cli.add_command(pipeline_temperature_analysis)
    cli.add_command(pipeline_variable_analysis)
    cli.add_command(module_adp_analysis)

    ### Modules for dev branch ###

    cli.add_command(pipeline_vp)
    cli.add_command(module_intensity_compare)
    cli.add_command(pipeline_intensity_compare)
    cli.add_command(pipeline_rigaku_vt)
    cli.add_command(pipeline_aus_synch_vt)
    cli.add_command(module_xds_cell_transformation)
    cli.add_command(module_xds_reprocess)
    cli.add_command(module_xprep)
    cli.add_command(pipeline_xds_reprocess)
    cli.add_command(pipeline_xprep)
    cli.add_command(pipeline_rotation_planes)
    cli.add_command(pipeline_position_analysis)
    cli.add_command(pipeline_AS_Brute)
    cli.add_command(module_molecule_reconstruction)
    cli.add_command(pipeline_shelxt_auto)
    cli.add_command(module_platon_squeeze)
    cli.add_command(pipeline_platon_squeeze)
    cli.add_command(pipeline_AS_Brute_individual)


def run() -> None:
    """This function sets up the logs and runs the commandline interface
    ALLLLLLLLLL of the above click commands are merely setting up the commandline interface
    It is actually executed by this function (which is called by typing 'cxasap' into the terminal
    Ie this is the function that the setup.py files defines as an entry point
    Returns:
        None
    """
    # Fun function but if used means that 'TAB to auto-complete' doesn't work. Can't handle the print and sleep statements
    os_test()
    # Set up global logs

    log_location = (
        pathlib.Path(os.path.abspath(__file__)).parent / "error_logs/error_output.txt"
    )

    logging.basicConfig(filename=log_location, level=logging.INFO)
    logging.info("CX-ASAP Started")

    cli()
