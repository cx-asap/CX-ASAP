#!/usr/bin/env python3

#############################################################################################################
# ----------------------------------------CX-ASAP: structural_analysis---------------------------------------#
# -----Authors: Amy J. Thompson, Tyler J. Philp, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price-----#
# -------------------------------------Python Implementation by AJT and TJP----------------------------------#
# ----------------------------------------Project Design by JRP and JKC--------------------------------------#
# -------------------------------------Valuable Coding Support by KMS & DJE----------------------------------#
#############################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config, Grapher
import os
import pathlib
import pandas as pd
import logging

# ----------Class Definition----------#


class Structural_Analysis:
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

        self.bond_paras = [
            "_geom_bond_atom_site_label_1",
            "_geom_bond_atom_site_label_2",
            "_geom_bond_distance",
        ]
        self.angle_paras = [
            "_geom_angle_atom_site_label_1",
            "_geom_angle_atom_site_label_2",
            "_geom_angle_atom_site_label_3",
            "_geom_angle",
        ]
        self.torsion_paras = [
            "_geom_torsion_atom_site_label_1",
            "_geom_torsion_atom_site_label_2",
            "_geom_torsion_atom_site_label_3",
            "_geom_torsion_atom_site_label_4",
            "_geom_torsion",
        ]
        self.hbond_paras = [
            "_geom_hbond_atom_site_label_D",
            "_geom_hbond_atom_site_label_H",
            "_geom_hbond_atom_site_label_A",
            "_geom_hbond_distance_DH",
            "_geom_hbond_distance_HA",
            "_geom_hbond_distance_DA",
            "_geom_hbond_angle_DHA",
        ]

    def import_and_analyse(
        self,
        bond_csv: str,
        angle_csv: str,
        torsion_csv: str,
        hbond_csv: str,
        atoms_for_analysis: list,
        location: str = False,
        prefix: str = "",
        flexible: bool = False,
        varying_parameter: str = "Data_Block",
    ) -> None:
        """Imports data from .csv files for analysis

        Args:
            bond_csv (str): full path to the .csv file with bond info
                            (false if analysis not wanted)
            angle_csv (str): full path to the .csv file with angle info
                            (false if analysis not wanted)
            torsion_csv (str): full path to the .csv file with torsion info
                            (false if analysis not wanted)
            hbond_csv (str): full path to the .csv file with torsion info
                (false if analysis not wanted)
            atoms_for_analysis (list): list of important atoms to separate
            location (str): full path to the location containing all .csv files
                            (if used by other cxasap pipeline)
            prefix (str): Adds a prefix to folder names
            flexible (bool): True for flexible mapping experiment, false if anything else
            varying_parameter (str): Which parameter is varying (ie _diffrn_ambient_temperature)
        """
        # Imports data from .csv files

        if location == False:
            try:
                self.location = pathlib.Path(bond_csv).parent
            except:
                try:
                    self.location = pathlib.Path(angle_csv).parent
                except:
                    try:
                        self.location = pathlib.Path(torsion_csv).parent
                    except:
                        try:
                            self.location = pathlib.Path(hbond_csv).parent
                        except:
                            pass
        else:
            self.location = location

        if bond_csv != False:
            bond_df = pd.read_csv(pathlib.Path(bond_csv))

            os.chdir(pathlib.Path(bond_csv).parent)

            self.structural_analysis(
                bond_df,
                self.bond_paras,
                "Bonds",
                "Bond_Lengths.csv",
                "Individual_Bond_Length_Data",
                "Length (Angstroms)",
                atoms_for_analysis,
                prefix,
                flexible,
                varying_parameter,
            )

        if angle_csv != False:
            angle_df = pd.read_csv(pathlib.Path(angle_csv))

            os.chdir(pathlib.Path(angle_csv).parent)

            self.structural_analysis(
                angle_df,
                self.angle_paras,
                "Angles",
                "Bond_Angles.csv",
                "Individual_Angle_Data",
                "Angle ($^\circ$)",
                atoms_for_analysis,
                prefix,
                flexible,
                varying_parameter,
            )
        if torsion_csv != False:
            try:
                torsion_df = pd.read_csv(pathlib.Path(torsion_csv))
            except FileNotFoundError:
                logging.info(
                    __name__
                    + " : No Torsion data - likely because structures not refined with CONF instruction"
                )
            else:
                os.chdir(pathlib.Path(torsion_csv).parent)

                self.structural_analysis(
                    torsion_df,
                    self.torsion_paras,
                    "Torsions",
                    "Bond_Torsions.csv",
                    "Individual_Torsion_Data",
                    "Angle ($^\circ$)",
                    atoms_for_analysis,
                    prefix,
                    flexible,
                    varying_parameter,
                )
        if hbond_csv != False:
            try:
                hbond_df = pd.read_csv(pathlib.Path(hbond_csv))
            except FileNotFoundError:
                logging.info(
                    __name__
                    + " : No Hbond data - likely because structures not refined with HTAB instruction"
                )
            else:
                os.chdir(pathlib.Path(hbond_csv).parent)

                self.structural_analysis(
                    hbond_df,
                    self.hbond_paras,
                    "Hbonds",
                    "Hbond_details.csv",
                    "Individual_Hbond_Data",
                    "Length (Angstroms)",
                    atoms_for_analysis,
                    prefix,
                    flexible,
                    varying_parameter,
                )

    def structural_analysis(
        self,
        df: "pd.DataFrame",
        components: list,
        structure_type: str,
        file_name: str,
        folder_name: str,
        y_unit: str,
        atoms_for_analysis: list,
        prefix: str,
        flexible: bool = False,
        varying_parameter: str = "Data_Block",
    ) -> None:
        """Performs structural analysis for bonds, angles, hbonds and torsions

        ONLY does the ones the user has choosen in the conf.yaml file

        Separates out the 'important' structural details based on which

        atoms the user wants to analyse

        Also will make a folder of 'individual' .csv files, which will separate

        into each bond/angle/torsion individually

        Finally graphs them

        *See comments in code for how it does this, there were some very annoying

        and fiddly things*

        Args:
            df ("pd.DataFrame"): dataframe containing the relevant data
            components (list): different CIF parameters involved in the structure parameter
                                see self.angle_paras for example
            structure_type (str): "Bonds", "Angles", "Hbonds" or "Torsions"
            file_name (str): name of the .csv file with data in it
            folder_name (str): name of the folder where individual data goes
            y_unit (str): label for y-axis
            atoms_for_analysis (list): list of important atoms to separate
            prefix (str): Adds a prefix to folder names
            flexible (bool): True for flexible mapping experiment, false if anything else
            varying_parameter (str): Which parameter is varying (ie _diffrn_ambient_temperature)
        """

        os.chdir(self.location)

        folder_name = prefix + "_" + folder_name

        graph = Grapher(self.test_mode)

        important_df = pd.DataFrame()

        if df.empty == False:
            # Separates out important atoms by looking for them in any column and merging into one dataframe

            for item in atoms_for_analysis:
                temp_df = df[df.eq(item).any(1)]
                important_df = pd.concat([important_df, temp_df], axis=0)

            # Need to make a new column of the indices - the above code will give you double ups if both atoms in a bond are "important"

            # When concatanating the dataframes, it keeps the indicies of the original dataframes

            # Meaning... doubleup indices AND out of order indices

            # Both are bad, so make a new column of the indicies, use it to drop duplicates, and then resets the index of the important dataframe

            important_df["index_2"] = list(important_df.index)

            important_df = important_df.drop_duplicates(subset=["index_2"])

            important_df = important_df.drop(["index_2"], axis=1)

            important_df = important_df.reset_index(drop=True)

            # Making symmetry equivalent bonds (ie same atom 1 and atom 2) distinguishable

            column_names = important_df.columns.values.tolist()

            # Put the atoms together into a 'joined' list to make it easier to sort through

            if structure_type == "Bonds":
                important_df["Joined"] = (
                    important_df[column_names[0]] + important_df[column_names[2]]
                )
            elif structure_type == "Angles":
                important_df["Joined"] = (
                    important_df[column_names[0]]
                    + important_df[column_names[2]]
                    + important_df[column_names[4]]
                )
            elif structure_type == "Torsions":
                important_df["Joined"] = (
                    important_df[column_names[0]]
                    + important_df[column_names[2]]
                    + important_df[column_names[4]]
                    + important_df[column_names[6]]
                )
            elif structure_type == "Hbonds":
                important_df["Joined"] = (
                    important_df[column_names[0]]
                    + important_df[column_names[2]]
                    + important_df[column_names[4]]
                )
            else:
                logging.info("Something went weird.")

            dup = important_df.duplicated(["Joined", varying_parameter], keep=False)

            # The below function counts the number of each group of duplicates

            list_dup = important_df.pivot_table(
                columns=["Joined", varying_parameter], aggfunc="size"
            ).to_dict()

            counter = 0

            new_column = []

            # This appends a suffix to each duplicated bond in the joined column based on how many there are

            for j, i in enumerate(dup):
                bond = important_df["Joined"][j]

                # NOTE HERE VARIABLE 'TEMPERATURE' CAN BE ANY PARAMETER BUT I DIDN'T WANT TO CHANGE WHOLE CODE

                temperature = important_df[varying_parameter][j]

                if i == True:
                    temp = important_df["Joined"][j]
                    new_column.append(important_df["Joined"][j] + "_" + str(counter))
                    counter += 1
                else:
                    new_column.append(important_df["Joined"][j])
                try:
                    if counter == list_dup[(bond, temperature)]:
                        counter = 0
                except UnboundLocalError:
                    pass

            important_df["Joined"] = new_column

            important_df.to_csv(prefix + "_Important_" + file_name, index=None)

            # Individual CSVs

            if structure_type == "Bonds":
                df["Joined"] = df[column_names[0]] + df[column_names[2]]
            elif structure_type == "Angles":
                df["Joined"] = (
                    df[column_names[0]] + df[column_names[2]] + df[column_names[4]]
                )
            elif structure_type == "Torsions":
                df["Joined"] = (
                    df[column_names[0]]
                    + df[column_names[2]]
                    + df[column_names[4]]
                    + df[column_names[6]]
                )
            elif structure_type == "Hbonds":
                df["Joined"] = (
                    df[column_names[0]] + df[column_names[2]] + df[column_names[4]]
                )
            else:
                logging.info("Something went weird.")

            discrete_atoms = list(dict.fromkeys(df["Joined"]))

            try:
                os.mkdir(folder_name)
            except FileExistsError:
                pass

            os.chdir(folder_name)

            for item in discrete_atoms:
                separated_df = df[df.eq(item).any(1)]
                separated_df.to_csv(
                    structure_type + "_" + str(item) + ".csv", index=None
                )

            os.chdir("..")

            # Make Graphs

            if flexible == True:
                important_df["_diffrn_ambient_temperature"] = list(
                    range(0, len(important_df["_diffrn_ambient_temperature"]))
                )
                important_df.rename(
                    columns={"_diffrn_ambient_temperature": "number"}, inplace=True
                )

            y_data = []
            y_headers = []
            x_data = []

            grouped = important_df.groupby("Joined")

            for i, g in grouped:
                y_headers.append(i)
                y_data.append(list(g[column_names[-6]]))
                if flexible == True:
                    x_data.append(list(g["number"]))
                    x_unit = "Structure number"
                x_data.append(list(g[varying_parameter]))
                if column_names[-4] == "_diffrn_ambient_temperature":
                    x_unit = "Temperature (K)"
                elif (
                    flexible == False
                    and column_names[-4] != "_diffrn_ambient_temperature"
                ):
                    x_unit = varying_parameter

            graph.single_scatter_graph(
                x_data[0],
                y_data,
                x_unit,
                y_unit,
                structure_type,
                prefix + "_" + structure_type + ".png",
                y_headers,
            )
