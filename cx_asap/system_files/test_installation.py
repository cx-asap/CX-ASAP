#!/usr/bin/env python3

###################################################################################################
# -------------------------------------CX-ASAP: Test Installation----------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

import errno
import logging
import os
import pathlib
import platform
import shutil
import sys
import time

import yaml
from cif_validation.modules.cif_merge import Cif_Merge
from cif_validation.modules.instrument_cif_generation import Instrument_CIF
from cif_validation.pipelines.cif_pipeline import CIF_Compile_Pipeline
from data_refinement.modules.refinement import Structure_Refinement
from data_refinement.pipelines.refine_pipeline import Refinement_Pipeline
from overall_pipelines.cxasap_pipeline import General_Pipeline
from post_refinement_analysis.modules.cell_analysis import Cell_Deformation
from post_refinement_analysis.modules.cif_read import CIF_Read
from post_refinement_analysis.modules.structural_analysis import Structural_Analysis
from post_refinement_analysis.pipelines.variable_cif_parameter import (
    Variable_Analysis_Pipeline,
)
from system_files.utils import Config, Nice_YAML_Dumper

# ----------Class Definition----------#


class Test:
    def __init__(self) -> None:

        """Initialises the class

        Sets up the yaml parameters input by the user

        Also defines the location for the system yaml file

        which stores a yaml of code-only parameters accessible throughout

        the software package
        """

        pass

        # Setup yaml files and logger

        config = Config(True)

        self.sys = config.sys
        self.sys_path = config.sys_path

        self.outcome = False

        self.messages = []

    def test_shelxl_module(self) -> None:

        """Tests if module-refinement is working

        Uses the included test data Sets

        Uses yaml parameters in the sys.yaml file
        """

        print("Testing module-refinement..... ")

        time.sleep(3)

        self.sys["module-refinement"]["structure_location"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["module-refinement"]["structure_location"])
            ).absolute()
        )

        self.sys["module-refinement"]["reference_path"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["module-refinement"]["reference_path"])
            ).absolute()
        )

        os.chdir(pathlib.Path(self.sys["module-refinement"]["structure_location"]).parent.parent)
        for item in os.listdir():
            if item == "CIF_Analysis" or item == "Refinement_Statistics":
                os.chdir(item)
                for i in os.listdir():
                    if os.path.isdir(i) == True:
                        os.chdir(i)
                        for k in os.listdir():
                            if os.path.isdir(k) == True:
                                os.chdir(k)
                                for j in os.listdir():
                                    os.remove(j)
                                os.chdir("..")
                                os.rmdir(k)
                            else:
                                os.remove(k)

                        os.chdir("..")
                        os.rmdir(i)
                    else:
                        os.remove(i)
                os.chdir("..")
                os.rmdir(item)

        with open(self.sys_path, "w") as f:
            yaml.dump(
                self.sys,
                f,
                default_flow_style=False,
                Dumper=Nice_YAML_Dumper,
                sort_keys=False,
            )

        shelxl = Structure_Refinement(test_mode=True)

        try:
            shelxl.run_shelxl(
                self.sys["module-refinement"]["structure_location"],
                self.sys["module-refinement"]["reference_path"],
                self.sys["module-refinement"]["refinements_to_check"],
                self.sys["module-refinement"]["tolerance"],
                self.sys["module-refinement"]["maximum_cycles"],
            )
        except Exception as error:
            # self.logger.info(f'module-refinement failed with error: {error}')
            logging.info(f"module-refinement failed with error: {error}")
            message = "module-refinement could not execute: open error log for details"
            self.outcome = False

        else:
            self.outcome = True
            message = "module-refinement functioning correctly for test dataset"

        print(message)
        self.messages += [message]

    def test_shelxl_pipeline(self) -> None:

        """Tests if pipeline-refinement is working

        Uses the included test data Sets

        Uses yaml parameters in the sys.yaml file
        """

        print("Testing pipeline-refinement..... ")

        time.sleep(3)

        self.sys["pipeline-refinement"]["experiment_location"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["pipeline-refinement"]["experiment_location"])
            ).absolute()
        )

        self.sys["pipeline-refinement"]["reference_path"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["pipeline-refinement"]["reference_path"])
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

        shelxl = Refinement_Pipeline(test_mode=True)

        try:
            shelxl.multiple_refinement(
                self.sys["pipeline-refinement"]["experiment_location"],
                self.sys["pipeline-refinement"]["reference_path"],
                self.sys["pipeline-refinement"]["experiment_location"],
                self.sys["pipeline-refinement"]["refinements_to_check"],
                self.sys["pipeline-refinement"]["tolerance"],
                self.sys["pipeline-refinement"]["maximum_cycles"],
            )
        except Exception as error:
            logging.info(f"pipeline-refinement failed with error: {error}")
            message = "pipeline-refinement could not execute: open error log for details"
            self.outcome = False
        else:
            self.outcome = True
            message = "pipeline-refinement functioning correctly for test dataset"

        print(message)
        self.messages += [message]

    def test_instrument_cif_generation_module(self) -> None:

        """Tests if module-make-instrument-cif is working

        Uses the included test data Sets

        Uses yaml parameters in the sys.yaml file
        """

        print("Testing instrument_cif_generation")

        time.sleep(3)

        self.sys["module-make-instrument-cif"]["reference_cif"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["module-make-instrument-cif"]["reference_cif"])
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

        instrument_cif = Instrument_CIF(test_mode=True)

        try:
            instrument_cif.read_reference_cif(self.sys["module-make-instrument-cif"]["reference_cif"])
            instrument_cif.make_instrument_cif()
        except Exception as error:
            logging.info(f"module-make-instrument-cif failed with error: {error}")
            message = "module-make-instrument-cif could not execute: open error log for details"
            self.outcome = False
        else:
            self.outcome = True
            message = "module-make-instrument-cif functioning correctly for test dataset"

        print(message)
        self.messages += [message]

    def test_cif_module(self) -> None:

        """Tests if module-cif-merge is working

        Uses the included test data Sets

        Uses yaml parameters in the sys.yaml file
        """

        print("Testing cif_merge")

        time.sleep(3)

        self.sys["module-cif-merge"]["instrument_cif"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent / (self.sys["module-cif-merge"]["instrument_cif"])
            ).absolute()
        )

        self.sys["module-cif-merge"]["new_cif"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent / (self.sys["module-cif-merge"]["new_cif"])
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

        finalise_cifs = Cif_Merge(test_mode=True)

        try:
            finalise_cifs.import_CIFs(
                pathlib.Path(self.sys["module-cif-merge"]["instrument_cif"]),
                pathlib.Path(self.sys["module-cif-merge"]["new_cif"]),
            )
            finalise_cifs.merge_CIFs()
            finalise_cifs.write_out(
                pathlib.Path(self.sys["module-cif-merge"]["new_cif"]).parent,
                "combined.cif",
                "check_CIF.chk",
                "200.cif",
            )
            finalise_cifs.validate_CIFs("combined.cif")

        except Exception as error:
            logging.info(f"module-cif-merge failed with error: {error}")
            message = "module-cif-merge could not execute: open error log for details"
            self.outcome = False
        else:
            self.outcome = True
            message = "module-cif-merge functioning correctly for test dataset"

        current_dir = os.getcwd()
        os.chdir(pathlib.Path(self.sys["module-cif-merge"]["new_cif"]).parent)
        for item in os.listdir():
            if item.startswith("combined"):
                os.remove(item)
        os.chdir(current_dir)

        print(message)
        self.messages += [message]

    def test_cif_pipeline(self) -> None:

        """Tests if pipeline-cif is working

        Uses the included test data Sets

        Uses yaml parameters in the sys.yaml file
        """

        print("Testing cif_pipeline")

        time.sleep(3)

        self.sys["pipeline-cif"]["experiment_location"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["pipeline-cif"]["experiment_location"])
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

        current_dir = os.getcwd()
        os.chdir(self.sys["pipeline-cif"]["experiment_location"])
        for item in os.listdir():
            if item == "combined.cif":
                os.remove(item)
        experiment_location = os.getcwd()
        for item in os.listdir():
            if os.path.isdir(item) == True:
                shutil.copy(
                    pathlib.Path(experiment_location).parent / "ref" / self.sys["pipeline-cif"]["instrument_file"],
                    item,
                )
        os.chdir(current_dir)

        cifs = CIF_Compile_Pipeline(test_mode=True)

        try:
            cifs.configure(
                self.sys["pipeline-cif"]["experiment_location"],
                self.sys["pipeline-cif"]["structure_solution"],
                self.sys["pipeline-cif"]["chemical_formula"],
                self.sys["pipeline-cif"]["crystal_habit"],
                self.sys["pipeline-cif"]["crystal_colour"],
                self.sys["pipeline-cif"]["max_crystal_dimension"],
                self.sys["pipeline-cif"]["middle_crystal_dimension"],
                self.sys["pipeline-cif"]["min_crystal_dimension"],
                self.sys["pipeline-cif"]["instrument_ending"],
                self.sys["pipeline-cif"]["instrument_file"],
            )
            cifs.compile_cifs(self.sys["pipeline-cif"]["experiment_location"])
        except Exception as error:
            logging.info(f"pipeline-cif failed with error: {error}")
            message = "pipeline-cif could not execute: open error log for details"
            self.outcome = False
        else:
            self.outcome = True
            message = "pipeline-cif functioning correctly for test dataset"

        print(message)
        self.messages += [message]

    def test_cif_read_module(self) -> None:

        """Tests if module-cif-read is working

        Uses the included test data Sets

        Uses yaml parameters in the sys.yaml file
        """

        print("Testing cif_read_module")

        time.sleep(3)

        self.sys["module-cif-read"]["folder_containing_cifs"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["module-cif-read"]["folder_containing_cifs"])
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

        analysis = CIF_Read(test_mode=True)

        try:

            analysis.configure(self.sys["module-cif-read"]["cif_parameters"])
            analysis.get_data(
                self.sys["module-cif-read"]["folder_containing_cifs"],
                self.sys["module-cif-read"]["structural_analysis_bonds"],
                self.sys["module-cif-read"]["structural_analysis_angles"],
                self.sys["module-cif-read"]["structural_analysis_torsions"],
                self.sys["module-cif-read"]["ADP_analysis"],
            )
            analysis.data_output()

        except Exception as error:
            self.logger.info(f"module-cif-read failed with error: {error}")
            logging.info(f"module-cif-read failed with error: {error}")
            message = "module-cif-read could not execute: open error log for details"
            self.outcome = False
        else:
            self.outcome = True
            message = "module-cif-read functioning correctly for test dataset"

        print(message)
        self.messages += [message]

    def test_cell_analysis_module(self) -> None:

        """Tests if module-cell-analysis is working

        Uses the included test data Sets

        Uses yaml parameters in the sys.yaml file
        """

        print("Testing cell_analysis_module")

        time.sleep(3)

        self.sys["module-cell-analysis"]["csv_location"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["module-cell-analysis"]["csv_location"])
            ).absolute()
        )

        self.sys["module-cell-analysis"]["reference_unit_cell"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["module-cell-analysis"]["reference_unit_cell"])
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

        analysis = Cell_Deformation(test_mode=True)

        try:

            analysis.import_data(
                self.sys["module-cell-analysis"]["csv_location"],
                self.sys["module-cell-analysis"]["reference_unit_cell"],
            )
            analysis.calculate_deformations()
            analysis.graphical_analysis(
                self.sys["module-cell-analysis"]["x_axis_header"],
                self.sys["module-cell-analysis"]["x_axis_header"],
            )

            # This is defined outside of the class, because when doing multiple analysis for a series of mapping experiments, the dictionary will be different, however, as this is a standalone module and not a pipeline for multiple structures, it is defined here

            try:
                test1 = analysis.df["_diffrn_reflns_av_R_equivalents"]
                test2 = analysis.df["_refine_ls_R_factor_gt"]
                test3 = analysis.df["_diffrn_measured_fraction_theta_full"]
            except:
                # self.logger.critical('No data quality statistics found in imported .csv file')
                logging.critical("No data quality statistics found in imported .csv file")
                print("Error! Check logs")
                exit()

            y_dict = {
                "R1": analysis.df["_refine_ls_R_factor_gt"],
                "Rint": analysis.df["_diffrn_reflns_av_R_equivalents"],
                "Completeness": analysis.df["_diffrn_measured_fraction_theta_full"],
            }

            analysis.quality_analysis(
                self.sys["module-cell-analysis"]["x_axis_header"],
                y_dict,
                self.sys["module-cell-analysis"]["x_axis_header"],
            )

        except Exception as error:
            logging.info(f"module-cell-analysis failed with error: {error}")
            message = "module-cell-analysis could not execute: open error log for details"
            self.outcome = False
        else:
            self.outcome = True
            message = "module-cell-analysis functioning correctly for test dataset"

        print(message)
        self.messages += [message]

    def test_structural_analysis_module(self) -> None:

        """Tests if module-structural-analysis is working

        Uses the included test data Sets

        Uses yaml parameters in the sys.yaml file
        """

        print("Testing module-structural-analysis")

        time.sleep(3)

        self.sys["module-structural-analysis"]["bond_data"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["module-structural-analysis"]["bond_data"])
            ).absolute()
        )

        self.sys["module-structural-analysis"]["angle_data"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["module-structural-analysis"]["angle_data"])
            ).absolute()
        )

        self.sys["module-structural-analysis"]["torsion_data"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["module-structural-analysis"]["torsion_data"])
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

        geometry = Structural_Analysis(test_mode=True)

        try:
            geometry.import_and_analyse(
                self.sys["module-structural-analysis"]["bond_data"],
                self.sys["module-structural-analysis"]["angle_data"],
                self.sys["module-structural-analysis"]["torsion_data"],
                self.sys["module-structural-analysis"]["atoms_for_analysis"],
            )
        except Exception as error:
            logging.info(f"module-structural-analysis failed with error: {error}")
            message = "module-structural-analysis could not execute: open error log for details"
            self.outcome = False
        else:
            self.outcome = True
            message = "module-structural-analysis functioning correctly for test dataset"

        print(message)
        self.messages += [message]

    def test_variable_cif_parameter_pipeline(self) -> None:

        """Tests if pipeline-variable-analysis is working

        Uses the included test data Sets

        Uses yaml parameters in the sys.yaml file
        """

        print("Testing pipeline_variable_cif_parameter")

        time.sleep(3)

        self.sys["pipeline-variable-analysis"]["experiment_location"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["pipeline-variable-analysis"]["experiment_location"])
            ).absolute()
        )

        self.sys["pipeline-variable-analysis"]["reference_unit_cell"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["pipeline-variable-analysis"]["reference_unit_cell"])
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

        vt_analysis = Variable_Analysis_Pipeline(test_mode=True)

        try:
            vt_analysis.analyse_data(
                self.sys["pipeline-variable-analysis"]["reference_unit_cell"],
                self.sys["pipeline-variable-analysis"]["experiment_location"],
                self.sys["pipeline-variable-analysis"]["cif_parameters"],
                self.sys["pipeline-variable-analysis"]["atoms_for_analysis"],
                self.sys["pipeline-variable-analysis"]["varying_cif_parameter"],
                self.sys["pipeline-variable-analysis"]["structural_analysis_bonds"],
                self.sys["pipeline-variable-analysis"]["structural_analysis_angles"],
                self.sys["pipeline-variable-analysis"]["structural_analysis_torsions"],
                self.sys["pipeline-variable-analysis"]["ADP_analysis"],
            )
        except Exception as error:
            logging.info(f"pipeline-variable-analysis failed with error: {error}")
            message = "pipeline-variable-analysis could not execute: open error log for details"
            self.outcome = False
        else:
            self.outcome = True
            message = "pipeline-variable-analysis functioning correctly for test dataset"

        print(message)
        self.messages += [message]

    def test_general_pipeline(self) -> None:

        """Tests if pipeline-general is working

        Uses the included test data Sets

        Uses yaml parameters in the sys.yaml file
        """

        print("Testing pipeline-general..... ")

        time.sleep(3)

        current_dir = os.getcwd()
        os.chdir(self.sys["pipeline-cif"]["experiment_location"])
        for item in os.listdir():
            if os.path.isdir(item) == False:
                os.remove(item)
            elif os.path.isdir(item) == True and "_Individual" in item:
                os.chdir(item)
                for i in os.listdir():
                    os.remove(i)
                os.chdir("..")
                os.rmdir(item)
            elif os.path.isdir(item) == True and "ADP" in item:
                os.chdir(item)
                for i in os.listdir():
                    os.remove(i)
                os.chdir("..")
                os.rmdir(item)
        os.chdir(current_dir)

        self.sys["pipeline-general"]["experiment_location"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["pipeline-general"]["experiment_location"])
            ).absolute()
        )

        self.sys["pipeline-general"]["reference_cif_location"] = str(
            pathlib.Path(
                pathlib.Path(os.path.abspath(__file__)).parent.parent
                / (self.sys["pipeline-general"]["reference_cif_location"])
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

        full = General_Pipeline(test_mode=True)

        try:

            full.make_dirs(self.sys["pipeline-general"]["experiment_location"])

            full.reference_extract(
                self.sys["pipeline-general"]["reference_cif_location"],
                self.sys["pipeline-general"]["experiment_location"],
                self.sys["pipeline-general"]["varying_parameter_values"],
                self.sys["pipeline-general"]["varying_cif_parameter"],
            )

            full.process(
                self.sys["pipeline-general"]["experiment_location"],
                full.reference_res,
                full.stats_location,
                self.sys["pipeline-general"]["refinements_to_check"],
                self.sys["pipeline-general"]["tolerance"],
                self.sys["pipeline-general"]["maximum_cycles"],
            )

            full.analyse(
                full.reference_res,
                self.sys["pipeline-general"]["experiment_location"],
                full.results_location,
                "cx-asap",
                full.chemical_formula,
                full.crystal_habit,
                full.crystal_colour,
                full.max_dimension,
                full.mid_dimension,
                full.min_dimension,
                self.sys["pipeline-general"]["structural_analysis_bonds"],
                self.sys["pipeline-general"]["structural_analysis_angles"],
                self.sys["pipeline-general"]["structural_analysis_torsions"],
                self.sys["pipeline-general"]["cif_parameters"],
                self.sys["pipeline-general"]["atoms_for_analysis"],
                self.sys["pipeline-general"]["varying_cif_parameter"],
                False,
                "instrument.cif",
                True,
                self.sys["pipeline-general"]["ADP_analysis"],
            )
        except Exception as error:
            logging.info(f"pipeline-refinement failed with error: {error}")
            message = "pipeline-general could not execute: open error log for details"
            self.outcome = False
        else:
            self.outcome = True
            message = "pipeline-general functioning correctly for test dataset"

        print(message)
        self.messages += [message]

    def test_third_party(self) -> None:
        """
        Function to test that third party executables required for cxasap are
        in the executable path: platon.exe, shelxl, shredcif
        Note: xds and xprep executables not yet required (beta release)
        """

        host_os = platform.system()
        if host_os == "Windows":
            required = ["platon.exe", "shelxl.exe", "shredcif.exe"]
            separator = "\\"
        else:
            required = ["platon", "shelxl", "shredcif"]
            separator = "/"

        exec_path = os.get_exec_path()

        matched = []
        for path in exec_path:
            for exec in required:
                available = os.access(path + separator + exec, os.X_OK)
                if available:
                    matched.append(exec)

        missing_execs = set(required).difference(matched)
        if missing_execs:
            message = "You are missing the following executables from your path: {}".format(missing_execs)
            self.outcome = False
        else:
            message = "Required third party softwares are available in your executable path"
            self.outcome = True

        print(message)
        self.messages += [message]

    def test_available_code(self) -> None:

        """Runs through all the tests available

        To set up a test for an additional module, add it to a function above

        and then add in a call to the test below

        Prints out a summary of the testing at the end

        """

        print(
            "WARNING: This test will take some time - it is ensuring that each module/mini-pipeline/overall pipeline is functioning as intended."
        )
        print("RECOMMENDATION: This is the perfect time for a tea/coffee break :)")

        time.sleep(5)

        self.test_third_party()

        if self.outcome == True:
            self.test_shelxl_module()
        if self.outcome == True:
            self.test_shelxl_pipeline()
        if self.outcome == True:
            self.test_instrument_cif_generation_module()
        if self.outcome == True:
            self.test_cif_module()
        if self.outcome == True:
            self.test_cif_pipeline()
        if self.outcome == True:
            self.test_cif_read_module()
        if self.outcome == True:
            self.test_cell_analysis_module()
        if self.outcome == True:
            self.test_structural_analysis_module()
        if self.outcome == True:
            self.test_variable_cif_parameter_pipeline()
        if self.outcome == True:
            self.test_general_pipeline()

        print("------------------------------------------------------- ")
        print("Summary of testing: ")
        print("------------------------------------------------------- ")
        for item in self.messages:
            print(item)

        print("------------------------------------------------------- ")
        print("Outcome of testing: ")
        print("------------------------------------------------------- ")

        if self.outcome == True:
            print("All tests for available modules completed successfully! Hooray!")
            print("You have correctly installed the software! Good job - love your work!")
            print("Now it is time for crystallography and data processing fun! :D ")
            print(
                "NOTE THAT TESTING IS NOT CURRENTLY SUPPORTED FOR MODULES AND PIPELINES INCLUDING DATA REDUCTION INCLUDING XPREP"
            )

        else:
            print("Uh oh you still have errors. Please see above summary!")


if __name__ == "__main__":
    tester = Test()
    tester.test_available_code()
