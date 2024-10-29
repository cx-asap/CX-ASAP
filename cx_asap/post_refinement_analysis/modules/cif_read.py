#!/usr/bin/env python3

###################################################################################################
# ----------------------------------------CX-ASAP: CIF Read----------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config, Directory_Browse
from CifFile import ReadCif
import yaml
import pandas as pd
import pathlib
import logging
from typing import Tuple

# ----------Class Definition----------#


class CIF_Read:
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

        self.structures_in_cif = []
        self.successful_positions = []
        self.results = {}
        self.errors = {}

        # Sets these to 0 to reset from previous runs

        self.sys["Structures_in_each_CIF"] = self.structures_in_cif
        self.sys["Successful_Positions"] = self.successful_positions

        with open(self.sys_path, "w") as f:
            yaml.dump(
                self.sys,
                f,
                default_flow_style=False,
                Dumper=Nice_YAML_Dumper,
                sort_keys=False,
            )

    def configure(self, search_items: list) -> None:
        """Sets up dictionaries and data frames to put data into

        Headers in dictionaries are desired parameters from conf.yaml

        Args:
            search_items (list): list of cif parameters to extract
        """

        # Pulls parameters from the configuration file as necessary, and uses it to set up an empty dataframe

        self.search_items = search_items

        for item in self.search_items:
            self.results[item] = []
            self.errors[item] = []
        self.data = pd.DataFrame()
        self.bond_data = pd.DataFrame()
        self.angle_data = pd.DataFrame()
        self.torsion_data = pd.DataFrame()
        self.hbond_data = pd.DataFrame()
        self.temp_df = pd.DataFrame()
        self.adp_data = pd.DataFrame()

    def parameter_tidy(self, raw: str, item: str) -> None:
        """This tidies the output from the CIF and separates

        the errors and the values and places them in separate dictionaires

        Args:
            raw (str): string extracted from CIF for tidying
            item (str): CIF parameter (ie '_cell_length_a')
        """

        if "(" in raw:
            temp = raw.split("(")
            temp2 = temp[1].strip(")")
            self.results[item].append(float(temp[0]))
            if "." in raw:
                temp3 = temp[0].split(".")
                self.errors[item].append(int(temp2) * 10 ** -(int(len(temp3[1]))))
            else:
                self.errors[item].append(int(temp2))
        else:
            try:
                self.results[item].append(float(raw))
            except ValueError:
                self.results[item].append(raw)
            self.errors[item].append(0)

    def generate_cif_list(self, df: "pd.Dataframe", counter: list) -> Tuple[list, list]:
        """This will generate a list of cifs that matches the length of

        data extracted from each CIF

        For example, if one cif has 10 datasets in it,

        then each parameter will appear 10 times

        This function will make a new list that has the name of the cif

        10 times to match the length of the data

        Args:
            df ("pd.Dataframe"): dataframe used to store the CIF data
            counter (list): length of this list will be the number of CIFs in the
                            combined CIF file

        Returns:
            longer_cif_list (list): list of CIFs matching the length of the df
            longer_data_blocks (list): list of datablocks (ie individual structures
                                        within the CIF) matching the length of the df
        """

        longer_cif_list = []
        longer_data_blocks = []
        if len(df) != 0:
            for index, item in enumerate(counter):
                i = 0
                while i < item:
                    i += 1
                    longer_cif_list.append(self.cif_list[index])
                    longer_data_blocks.append(self.data_blocks[index])

        return longer_cif_list, longer_data_blocks

    def get_data(
        self,
        location: str,
        bonds: bool = False,
        angles: bool = False,
        torsions: bool = False,
        hbonds: bool = False,
        adp: bool = False,
        varying_parameter: str = "_diffrn_ambient_temperature",
    ) -> None:
        """Searches through all folders in current working directory for CIFs

        For each CIF identified, extracts desired parameters

        Args:
            location (str): full path to the folder containing all CIFs
            bonds (bool): whether or not bond analysis should be run
            angles (bool): whether or not angle analysis should be run
            torsions (bool): whether or not torsion analysis should be run
            hbonds (bool): whether or not Hbond analysis should be run
            adps (bool): whether or not ADP analysis should be run

        """

        # This function searches through all of the folders in the current working directory for a cif file

        self.tree_browse = Directory_Browse(pathlib.Path(location), self.test_mode)

        self.tree_browse.enter_directory_multiple(pathlib.Path(location), ".cif")

        # For all found cif_files:

        for index, item in enumerate(self.tree_browse.item_files):
            cif_file = self.tree_browse.item_files[index].absolute()

            # extracts the desired cif parameters, as well as how many structures per cif and which positions were successful

            (
                temp_data,
                structures_in_cif_tmp,
                successful_positions_tmp,
            ) = self.data_harvest(cif_file, self.search_items, varying_parameter)

            self.structural_analysis(
                cif_file, bonds, angles, torsions, hbonds, varying_parameter
            )
            self.adp_analysis(cif_file, adp)

            # self.data = self.data.append(temp_data)
            self.data = pd.concat([self.data, temp_data])
            self.structures_in_cif.append(structures_in_cif_tmp)
            for item in successful_positions_tmp:
                self.successful_positions.append(item.strip("structure_"))

        self.sys["Structures_in_each_CIF"] = self.structures_in_cif
        self.sys["Successful_Positions"] = self.successful_positions

        with open(self.sys_path, "w") as f:
            yaml.dump(
                self.sys,
                f,
                default_flow_style=False,
                Dumper=Nice_YAML_Dumper,
                sort_keys=False,
            )

    def structural_analysis(
        self,
        cif_file: str,
        bonds: bool = False,
        angles: bool = False,
        torsions: bool = False,
        hbonds: bool = False,
        varying_parameter: str = "_diffrn_ambient_temperature",
    ) -> None:
        """Extracts structural information from CIF

        Args:
            cif_file (str): full path to the CIF for analysis
            bonds (bool): whether or not bond analysis should be run
            angles (bool): whether or not angle analysis should be run
            torsions (bool): whether or not torsion analysis should be run
            hbonds (bool): whether or not Hbond analysis should be run
        """

        # separate function for structural analysis as it is not always required

        bond_paras = [
            "_geom_bond_atom_site_label_1",
            "_geom_bond_atom_site_label_2",
            "_geom_bond_distance",
            varying_parameter,
        ]
        angle_paras = [
            "_geom_angle_atom_site_label_1",
            "_geom_angle_atom_site_label_2",
            "_geom_angle_atom_site_label_3",
            "_geom_angle",
            varying_parameter,
        ]
        torsion_paras = [
            "_geom_torsion_atom_site_label_1",
            "_geom_torsion_atom_site_label_2",
            "_geom_torsion_atom_site_label_3",
            "_geom_torsion_atom_site_label_4",
            "_geom_torsion",
            varying_parameter,
        ]
        hbond_paras = [
            "_geom_hbond_atom_site_label_D",
            "_geom_hbond_atom_site_label_H",
            "_geom_hbond_atom_site_label_A",
            "_geom_hbond_distance_DH",
            "_geom_hbond_distance_HA",
            "_geom_hbond_distance_DA",
            "_geom_hbond_angle_DHA",
            varying_parameter,
        ]

        # harvests data for structural information

        if bonds == True:
            (
                temp_data_bonds,
                structures_in_cif_tmp_bonds,
                successful_positions_tmp_bonds,
            ) = self.data_harvest(cif_file, bond_paras, varying_parameter)
            # self.bond_data = self.bond_data.append(temp_data_bonds)
            self.bond_data = pd.concat([self.bond_data, temp_data_bonds])
        if angles == True:
            (
                temp_data_angles,
                structures_in_cif_tmp_angles,
                successful_positions_tmp_angles,
            ) = self.data_harvest(cif_file, angle_paras, varying_parameter)
            # self.angle_data = self.angle_data.append(temp_data_angles)
            self.angle_data = pd.concat([self.angle_data, temp_data_angles])
        if torsions == True:
            (
                temp_data_torsions,
                structures_in_cif_tmp_torsions,
                successful_positions_tmp_torsions,
            ) = self.data_harvest(cif_file, torsion_paras, varying_parameter)
            # self.torsion_data = self.torsion_data.append(temp_data_torsions)
            self.torsion_data = pd.concat([self.torsion_data, temp_data_torsions])
        if hbonds == True:
            (
                temp_data_hbonds,
                structures_in_cif_tmp_hbonds,
                successful_positions_tmp_hbonds,
            ) = self.data_harvest(cif_file, hbond_paras, varying_parameter)
            # self.hbond_data = self.hbond_data.append(temp_data_hbonds)
            self.hbond_data = pd.concat([self.hbond_data, temp_data_hbonds])

    def adp_analysis(
        self,
        cif_file: str,
        adp: bool = False,
        varying_parameter: str = "_diffrn_ambient_temperature",
    ) -> None:
        """Extracts ADP information from CIF

        Args:
            cif_file (str): full path to the CIF for analysis
            adp (bool): whether or not ADP analysis should be run
        """

        adps = [
            "_atom_site_aniso_label",
            "_atom_site_aniso_U_11",
            "_atom_site_aniso_U_22",
            "_atom_site_aniso_U_33",
            "_atom_site_aniso_U_23",
            "_atom_site_aniso_U_13",
            "_atom_site_aniso_U_12",
        ]

        if adp == True:
            (
                temp_data_adps,
                structures_in_cif_tmp_adps,
                successful_positions_tmp_adps,
            ) = self.data_harvest(cif_file, adps, varying_parameter)
            # self.adp_data = self.adp_data.append(temp_data_adps)
            self.adp_data = pd.concat([self.adp_data, temp_data_adps])

    def data_harvest(
        self,
        cif_file: str,
        search_items: list,
        varying_parameter: str = "_diffrn_ambient_temperature",
    ) -> Tuple["pd.DataFrame", int, list]:
        """Extracts all other desired parameters from CIF

        Ie cell information and quality statistics

        NOT structural or ADP information

        Args:
            cif_file (str): full path to the CIF for analysis
            search_items (list): list of items to extract from CIF
                                in proper CIF syntax (ie "_cell_length_a")

        Returns:
            temp_df ("pd.DataFrame"): data from single CIF analysed
                                     gets merged into a larger dataframe containing
                                    data from multiple CIF files
            number_of_structures (int): how many datablocks are in CIF
            data_blocks (list): headers of datablocks in the CIF
        """

        for item in search_items:
            self.results[item] = []
            self.errors[item] = []
        self.temp_df = pd.DataFrame()

        # Use of the PyCifRW library for easy parsing of CIF Files

        cif = ReadCif(cif_file.name)

        # Identifies datablocks within the CIF File

        self.data_blocks = cif.keys()
        number_of_structures = len(self.data_blocks)

        structure_analysis_counter = {}

        for item in search_items:
            logging.info(
                __name__
                + " : File "
                + cif_file.stem
                + " found. Searching for parameter "
                + item
            )

            # This parameter is blanked every loop because each data_block is searched multiple times

            self.cif_list = []
            structure_analysis_counter[item] = []

            for experiment in self.data_blocks:
                try:
                    raw = cif[experiment][item]
                except:
                    logging.critical("Failed to find " + item + " in " + cif_file.stem)
                    print("Critical Failure - see error log for details")
                    exit()

                logging.info(
                    __name__ + " : File " + cif_file.stem + " - found parameter " + item
                )
                self.cif_list.append(cif_file.stem)

                # Checks for errors and separates values / converts error to actual value

                if type(raw) == list:
                    structure_analysis_counter[item] += [len(raw)]
                else:
                    structure_analysis_counter[item] += [1]

                if type(raw) != list:
                    self.parameter_tidy(raw, item)
                else:
                    for param in raw:
                        self.parameter_tidy(param, item)

        # All of this annoying code is to take into account the fact that one datablock will give one temperature, but many many bond lengths, so this really just makes sure that the lenghts of the temperature and cif name lists are of the correct length

        equivalent = True

        test_val = list(structure_analysis_counter.values())[0]

        for i in structure_analysis_counter:
            if structure_analysis_counter[i] != test_val:
                equivalent = False

        for item in search_items:
            if item == varying_parameter and equivalent == False:
                numbers_to_multiply = self.results[item]
                errors_to_multiply = self.errors[item]
                self.results[item] = []
                self.errors[item] = []
                for index, i in enumerate(test_val):
                    self.results[item] += [numbers_to_multiply[index]] * i
                    self.errors[item] += [errors_to_multiply[index]] * i

        self.temp_df["CIF_File"] = self.cif_list

        for para in search_items:
            if len(self.results[para]) != len(self.cif_list):
                self.temp_df = pd.DataFrame()

        for para in search_items:
            self.temp_df[para] = self.results[para]
            self.temp_df[para + "_error"] = self.errors[para]

        if len(self.temp_df) != len(self.cif_list):
            (
                self.temp_df["CIF_File"],
                self.temp_df["Data_Block"],
            ) = self.generate_cif_list(self.temp_df, test_val)

        return self.temp_df, number_of_structures, self.data_blocks

    def data_output(self) -> None:
        """Outputs all data to .csv files"""

        self.data.to_csv("CIF_Parameters.csv", index=None)
        if len(self.bond_data) != 0:
            self.bond_data.to_csv("Bond_Lengths.csv", index=None)
        if len(self.angle_data) != 0:
            self.angle_data.to_csv("Bond_Angles.csv", index=None)
        if len(self.torsion_data) != 0:
            self.torsion_data.to_csv("Bond_Torsions.csv", index=None)
        if len(self.hbond_data) != 0:
            self.hbond_data.to_csv("HBond_details.csv", index=None)
        if len(self.adp_data) != 0:
            self.adp_data.to_csv("ADPs.csv", index=None)
