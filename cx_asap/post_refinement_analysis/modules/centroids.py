#!/usr/bin/env python3

###################################################################################################
# ----------------------------------------CX-ASAP: centroids--------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Config
from system_files.crystal_math import orthonorm_matrix, parse_symm_op
from post_refinement_analysis.modules.lst_read import LST_Read
import pathlib
import os
import logging
import numpy as np
import pandas as pd

# ----------Class Definition----------#


class Centroids:
    """Calculates centroid positions and inter-centroid distances from .lst files.

    Atom fractional coordinates are parsed from the embedded .res block of the
    SHELXL .lst file. Distances are computed in Cartesian coordinates using the
    unit cell orthonormalisation matrix.

    This class can also extract the inter-plane angle reported by SHELXL when
    two MPLA commands are present.
    """

    def __init__(self, test_mode: bool = False) -> None:
        """Initialises the class.

        Args:
            test_mode (bool): if True, skips conf.yaml loading
        """

        self.test_mode = test_mode

        config = Config(self.test_mode)
        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

        self.lst_reader = LST_Read(self.test_mode)

    def grab_cell(self, file_name: str) -> None:
        """Reads the unit cell parameters from the CELL line of a .lst file.

        Args:
            file_name (str): full path to the .lst file
        """

        self.bad_flag = False

        with open(file_name, "rt") as f:
            split_line = []
            for line in f:
                if line.startswith(" CELL") or line.startswith("CELL"):
                    split_line = line.split()
                    break

        if len(split_line) >= 8:
            self.cell_params = [float(split_line[i]) for i in range(2, 8)]
        else:
            logging.warning(
                __name__ + " : Could not read CELL line from " + str(file_name)
            )
            self.bad_flag = True

    def calculate_centroid(
        self, coords: dict, atom_names: list, symmetry: str = None
    ) -> "np.ndarray | None":
        """Calculates the centroid (mean fractional position) of a group of atoms.

        If a symmetry operation string is provided, each atom's fractional
        coordinates are transformed by that operation before averaging.

        Args:
            coords (dict): fractional coordinates from LST_Read.extract_atom_coordinates()
            atom_names (list): atom labels to include in the centroid calculation
            symmetry (str): optional SHELXL/CIF symmetry operation string,
                            e.g. "-x+1/2, y+1/2, -z+1/2"

        Returns:
            centroid (np.ndarray): 1D array [x, y, z] in fractional coordinates,
                                   or None if no atoms were found
        """

        if isinstance(atom_names, str):
            atom_names = atom_names.split()

        R, t = parse_symm_op(symmetry) if symmetry else (None, None)

        positions = []
        for name in atom_names:
            key = name.upper()
            if key in coords:
                pos = np.array(coords[key], dtype=float)
                if R is not None:
                    pos = np.dot(R.T, pos) + t
                positions.append(pos)
            else:
                logging.warning(
                    __name__ + f" : Atom {name} not found in coordinate list"
                )

        if not positions:
            logging.critical(
                __name__ + " : No valid atoms found for centroid calculation"
            )
            return None

        return np.mean(positions, axis=0)

    def centroid_distance(
        self,
        coords: dict,
        atom_list_1: list,
        atom_list_2: list,
        symmetry_1: str = None,
        symmetry_2: str = None,
    ) -> float:
        """Calculates the distance in Angstroms between the centroids of two atom groups.

        Centroids are calculated in fractional coordinates then converted to Cartesian
        space using the unit cell orthonormalisation matrix before computing the
        Euclidean distance.

        Args:
            coords (dict): fractional coordinates from LST_Read.extract_atom_coordinates()
            atom_list_1 (list): atom labels for the first group
            atom_list_2 (list): atom labels for the second group
            symmetry_1 (str): optional symmetry operation string for centroid 1
            symmetry_2 (str): optional symmetry operation string for centroid 2

        Returns:
            distance (float): distance in Angstroms between the two centroids,
                              or 0.0 if either centroid could not be calculated
        """

        c1 = self.calculate_centroid(coords, atom_list_1, symmetry_1)
        c2 = self.calculate_centroid(coords, atom_list_2, symmetry_2)

        if c1 is None or c2 is None:
            return 0.0

        M = orthonorm_matrix(self.cell_params)

        # fractional row vector -> Cartesian: cart = frac @ M.T
        cart1 = np.dot(c1, M.T)
        cart2 = np.dot(c2, M.T)

        return float(np.linalg.norm(cart1 - cart2))

    def find_centroid_distance(
        self,
        file_name: str,
        atom_list_1: list,
        atom_list_2: list,
        symmetry_1: str = None,
        symmetry_2: str = None,
    ) -> float:
        """Reads a .lst file and calculates the centroid-to-centroid distance.

        Args:
            file_name (str): full path to the .lst file
            atom_list_1 (list): atom labels for the first group
            atom_list_2 (list): atom labels for the second group
            symmetry_1 (str): optional symmetry operation string for centroid 1
            symmetry_2 (str): optional symmetry operation string for centroid 2

        Returns:
            distance (float): distance in Angstroms between the two centroids
        """

        data = self.lst_reader.read(file_name)
        coords = self.lst_reader.extract_atom_coordinates(data)
        return self.centroid_distance(
            coords, atom_list_1, atom_list_2, symmetry_1, symmetry_2
        )

    def analysis_centroid_distance(
        self,
        lst_name: str,
        structure_number: int,
        results_path: str,
        atom_list_1: list,
        atom_list_2: list,
        label: str = "Centroid Distance",
        symmetry_1: str = None,
        symmetry_2: str = None,
    ) -> None:
        """Calculates the centroid distance between two atom groups and appends to a .csv.

        The output csv filename is derived from the label parameter.

        Args:
            lst_name (str): full path to the .lst file
            structure_number (int): structure number as the independent variable
            results_path (str): full path to the output results directory
            atom_list_1 (list): atom labels for the first group
            atom_list_2 (list): atom labels for the second group
            label (str): column header for the distance in the output csv;
                         also used to derive the csv filename
            symmetry_1 (str): optional symmetry operation string for centroid 1
            symmetry_2 (str): optional symmetry operation string for centroid 2
        """

        if lst_name == "":
            logging.info(__name__ + " : Refinement failed, no structure to analyse")
            return

        self.grab_cell(pathlib.Path(lst_name))

        if self.bad_flag:
            return

        distance = self.find_centroid_distance(
            pathlib.Path(lst_name), atom_list_1, atom_list_2, symmetry_1, symmetry_2
        )

        df = pd.DataFrame({"Structure": [structure_number], label: [distance]})
        csv_name = label.lower().replace(" ", "_") + ".csv"

        os.chdir(results_path)
        try:
            old_data = pd.read_csv(csv_name)
        except FileNotFoundError:
            df.to_csv(csv_name, index=None)
        else:
            new_df = pd.concat([old_data, df])
            new_df.to_csv(csv_name, index=None)
