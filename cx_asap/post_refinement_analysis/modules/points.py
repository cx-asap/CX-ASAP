#!/usr/bin/env python3

###################################################################################################
# -----------------------------------------CX-ASAP: points----------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Config
from system_files.crystal_math import (
    parse_symm_op,
    fractional_to_cartesian,
    distance_between_points,
    angle_between_points,
    torsion_between_points,
    point_to_plane_distance,
)
from post_refinement_analysis.modules.lst_read import LST_Read
import pathlib
import os
import logging
import numpy as np
import pandas as pd

# ----------Class Definition----------#


class PointGeometryEngine:
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

    def _normalise_atom_names(
        self, atom_names: "str | list[str] | tuple | set | np.ndarray"
    ) -> list:
        """Normalises atom label input into a list of strings."""

        if isinstance(atom_names, str):
            atom_names = atom_names.split()
        elif isinstance(atom_names, (list, tuple, set, np.ndarray)):
            atom_names = list(atom_names)
        else:
            message = (
                "atom_names must be a list/tuple/set/ndarray of labels or "
                "a whitespace-delimited string"
            )
            logging.error(__name__ + f" : {message}. Got {type(atom_names).__name__}")
            raise TypeError(message)

        return [str(item) for item in atom_names if str(item).strip()]

    def _extract_positions(
        self,
        coords: dict,
        atom_names: "str | list[str] | tuple | set | np.ndarray",
        symmetry: str = None,
    ) -> list:
        """Gets transformed fractional coordinates for each requested atom."""

        atom_names = self._normalise_atom_names(atom_names)
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

        return positions

    def calculate_point_group_center(
        self,
        coords: dict,
        atom_names: "str | list[str] | tuple | set | np.ndarray",
        symmetry: str = None,
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

        positions = self._extract_positions(coords, atom_names, symmetry)

        if not positions:
            logging.critical(
                __name__ + " : No valid atoms found for centroid calculation"
            )
            return None

        return np.mean(positions, axis=0)

    def calculate_point(
        self,
        coords: dict,
        atom_names: "str | list[str] | tuple | set | np.ndarray",
        symmetry: str = None,
        group_center_round_dp: int = None,
    ) -> "np.ndarray | None":
        """Calculates a single point from atom input.

        One atom label resolves to that atom position.
        Multiple atom labels resolve to the centroid of those atoms.
        """

        atom_names = self._normalise_atom_names(atom_names)
        if not atom_names:
            logging.critical(__name__ + " : No atoms supplied for point calculation")
            return None

        if len(atom_names) == 1:
            positions = self._extract_positions(coords, atom_names, symmetry)
            if not positions:
                logging.critical(
                    __name__ + " : Could not resolve atom for point calculation"
                )
                return None
            return positions[0]

        centroid = self.calculate_point_group_center(coords, atom_names, symmetry)
        if centroid is None:
            return None

        if group_center_round_dp is not None:
            return np.round(centroid, group_center_round_dp)

        return centroid

    def _point_uses_group_center(
        self, atom_names: "str | list[str] | tuple | set | np.ndarray"
    ) -> bool:
        """Returns True when a point definition resolves via centroid averaging."""

        return len(self._normalise_atom_names(atom_names)) > 1

    def _fractional_to_cartesian(self, point_frac: "np.ndarray | list") -> np.ndarray:
        """Converts one fractional point into Cartesian coordinates."""

        return fractional_to_cartesian(point_frac, self.cell_params)

    def point_distance(
        self,
        coords: dict,
        point_1_atoms: "str | list[str] | tuple | set | np.ndarray",
        point_2_atoms: "str | list[str] | tuple | set | np.ndarray",
        symmetry_1: str = None,
        symmetry_2: str = None,
        group_center_round_dp: int = None,
    ) -> float:
        """Calculates Cartesian distance between two points."""

        p1 = self.calculate_point(coords, point_1_atoms, symmetry_1, group_center_round_dp)
        p2 = self.calculate_point(coords, point_2_atoms, symmetry_2, group_center_round_dp)

        if p1 is None or p2 is None:
            return 0.0

        cart1 = self._fractional_to_cartesian(p1)
        cart2 = self._fractional_to_cartesian(p2)
        return distance_between_points(cart1, cart2)

    def point_angle(
        self,
        coords: dict,
        point_1_atoms: "str | list[str] | tuple | set | np.ndarray",
        point_2_atoms: "str | list[str] | tuple | set | np.ndarray",
        point_3_atoms: "str | list[str] | tuple | set | np.ndarray",
        symmetry_1: str = None,
        symmetry_2: str = None,
        symmetry_3: str = None,
        group_center_round_dp: int = None,
    ) -> float:
        """Calculates angle (degrees) formed by points 1-2-3 at point 2."""

        p1 = self.calculate_point(coords, point_1_atoms, symmetry_1, group_center_round_dp)
        p2 = self.calculate_point(coords, point_2_atoms, symmetry_2, group_center_round_dp)
        p3 = self.calculate_point(coords, point_3_atoms, symmetry_3, group_center_round_dp)

        if p1 is None or p2 is None or p3 is None:
            return 0.0

        c1 = self._fractional_to_cartesian(p1)
        c2 = self._fractional_to_cartesian(p2)
        c3 = self._fractional_to_cartesian(p3)

        try:
            return angle_between_points(c1, c2, c3)
        except ValueError as error:
            logging.warning(__name__ + f" : {error}")
            return 0.0

    def point_torsion(
        self,
        coords: dict,
        point_1_atoms: "str | list[str] | tuple | set | np.ndarray",
        point_2_atoms: "str | list[str] | tuple | set | np.ndarray",
        point_3_atoms: "str | list[str] | tuple | set | np.ndarray",
        point_4_atoms: "str | list[str] | tuple | set | np.ndarray",
        symmetry_1: str = None,
        symmetry_2: str = None,
        symmetry_3: str = None,
        symmetry_4: str = None,
        group_center_round_dp: int = None,
    ) -> float:
        """Calculates torsion angle (degrees) for points 1-2-3-4."""

        point_1_list = self._normalise_atom_names(point_1_atoms)
        point_2_list = self._normalise_atom_names(point_2_atoms)
        point_3_list = self._normalise_atom_names(point_3_atoms)
        point_4_list = self._normalise_atom_names(point_4_atoms)

        p1 = self.calculate_point(
            coords, point_1_list, symmetry_1, group_center_round_dp
        )
        p2 = self.calculate_point(
            coords, point_2_list, symmetry_2, group_center_round_dp
        )
        p3 = self.calculate_point(
            coords, point_3_list, symmetry_3, group_center_round_dp
        )
        p4 = self.calculate_point(
            coords, point_4_list, symmetry_4, group_center_round_dp
        )

        if p1 is None or p2 is None or p3 is None or p4 is None:
            return 0.0

        c1 = self._fractional_to_cartesian(p1)
        c2 = self._fractional_to_cartesian(p2)
        c3 = self._fractional_to_cartesian(p3)
        c4 = self._fractional_to_cartesian(p4)

        try:
            torsion_value = torsion_between_points(c1, c2, c3, c4)
        except ValueError as error:
            logging.warning(__name__ + f" : {error}")
            return 0.0

        has_symmetry = any(
            item is not None for item in [symmetry_1, symmetry_2, symmetry_3, symmetry_4]
        )
        has_centroid_point = any(
            len(item) > 1
            for item in [point_1_list, point_2_list, point_3_list, point_4_list]
        )

        if has_symmetry and has_centroid_point and abs(torsion_value) < 90.0:
            return torsion_value - 180.0 if torsion_value > 0.0 else torsion_value + 180.0

        return torsion_value

    def point_plane_distance(
        self,
        coords: dict,
        point_atoms: "str | list[str] | tuple | set | np.ndarray",
        plane_atoms: "str | list[str] | tuple | set | np.ndarray",
        point_symmetry: str = None,
        plane_symmetry: str = None,
        group_center_round_dp: int = None,
    ) -> float:
        """Calculates absolute distance from a point to a best-fit plane."""

        point_frac = self.calculate_point(
            coords, point_atoms, point_symmetry, group_center_round_dp
        )
        if point_frac is None:
            return 0.0

        plane_positions = self._extract_positions(coords, plane_atoms, plane_symmetry)
        if len(plane_positions) < 3:
            logging.warning(__name__ + " : Need at least 3 atoms to define a plane")
            return 0.0

        point_cart = self._fractional_to_cartesian(point_frac)
        plane_cart = np.array(
            [self._fractional_to_cartesian(item) for item in plane_positions],
            dtype=float,
        )

        try:
            return point_to_plane_distance(point_cart, plane_cart)
        except ValueError as error:
            logging.warning(__name__ + f" : {error}")
            return 0.0

    def _append_measurements_csv(
        self,
        csv_name: str,
        structure_number: int,
        values: "dict[str, float]",
        results_path: str,
    ) -> None:
        """Appends one structure worth of labelled measurements to a CSV."""

        if not values:
            return

        row = {"Structure": [structure_number]}
        for label, value in values.items():
            row[label] = [value]

        df = pd.DataFrame(row)

        os.chdir(results_path)
        try:
            old_data = pd.read_csv(csv_name)
        except FileNotFoundError:
            df.to_csv(csv_name, index=None)
        else:
            new_df = pd.concat([old_data, df])
            new_df.to_csv(csv_name, index=None)

    def analyse_point_geometry(
        self,
        lst_name: str,
        structure_number: int,
        results_path: str,
        distance_definitions: "list[dict]" = None,
        angle_definitions: "list[dict]" = None,
        torsion_definitions: "list[dict]" = None,
        plane_distance_definitions: "list[dict]" = None,
        mercury_output: bool = False,
    ) -> None:
        """Calculates configured point-geometry values and appends CSV outputs."""

        if lst_name == "":
            logging.info(__name__ + " : Refinement failed, no structure to analyse")
            return

        self.grab_cell(pathlib.Path(lst_name))
        if self.bad_flag:
            return

        data = self.lst_reader.read(lst_name)
        coords = self.lst_reader.extract_atom_coordinates(data)

        distance_values = {}
        mercury_distance_values = {}
        angle_values = {}
        mercury_angle_values = {}
        torsion_values = {}
        mercury_torsion_values = {}
        plane_distance_values = {}
        mercury_plane_distance_values = {}

        for index, definition in enumerate(distance_definitions or []):
            if not isinstance(definition, dict):
                continue
            label = definition.get("label", f"Distance_{index + 1}")
            distance_values[label] = self.point_distance(
                coords,
                definition.get("point_1_atoms", []),
                definition.get("point_2_atoms", []),
                definition.get("point_1_symmetry") or None,
                definition.get("point_2_symmetry") or None,
            )
            if mercury_output and any(
                self._point_uses_group_center(definition.get(key, []))
                for key in ["point_1_atoms", "point_2_atoms"]
            ):
                mercury_distance_values[label] = self.point_distance(
                    coords,
                    definition.get("point_1_atoms", []),
                    definition.get("point_2_atoms", []),
                    definition.get("point_1_symmetry") or None,
                    definition.get("point_2_symmetry") or None,
                    group_center_round_dp=3,
                )

        for index, definition in enumerate(angle_definitions or []):
            if not isinstance(definition, dict):
                continue
            label = definition.get("label", f"Angle_{index + 1}")
            angle_values[label] = self.point_angle(
                coords,
                definition.get("point_1_atoms", []),
                definition.get("point_2_atoms", []),
                definition.get("point_3_atoms", []),
                definition.get("point_1_symmetry") or None,
                definition.get("point_2_symmetry") or None,
                definition.get("point_3_symmetry") or None,
            )
            if mercury_output and any(
                self._point_uses_group_center(definition.get(key, []))
                for key in ["point_1_atoms", "point_2_atoms", "point_3_atoms"]
            ):
                mercury_angle_values[label] = self.point_angle(
                    coords,
                    definition.get("point_1_atoms", []),
                    definition.get("point_2_atoms", []),
                    definition.get("point_3_atoms", []),
                    definition.get("point_1_symmetry") or None,
                    definition.get("point_2_symmetry") or None,
                    definition.get("point_3_symmetry") or None,
                    group_center_round_dp=3,
                )

        for index, definition in enumerate(torsion_definitions or []):
            if not isinstance(definition, dict):
                continue
            label = definition.get("label", f"Torsion_{index + 1}")
            torsion_values[label] = self.point_torsion(
                coords,
                definition.get("point_1_atoms", []),
                definition.get("point_2_atoms", []),
                definition.get("point_3_atoms", []),
                definition.get("point_4_atoms", []),
                definition.get("point_1_symmetry") or None,
                definition.get("point_2_symmetry") or None,
                definition.get("point_3_symmetry") or None,
                definition.get("point_4_symmetry") or None,
            )
            if mercury_output and any(
                self._point_uses_group_center(definition.get(key, []))
                for key in [
                    "point_1_atoms",
                    "point_2_atoms",
                    "point_3_atoms",
                    "point_4_atoms",
                ]
            ):
                mercury_torsion_values[label] = self.point_torsion(
                    coords,
                    definition.get("point_1_atoms", []),
                    definition.get("point_2_atoms", []),
                    definition.get("point_3_atoms", []),
                    definition.get("point_4_atoms", []),
                    definition.get("point_1_symmetry") or None,
                    definition.get("point_2_symmetry") or None,
                    definition.get("point_3_symmetry") or None,
                    definition.get("point_4_symmetry") or None,
                    group_center_round_dp=3,
                )

        for index, definition in enumerate(plane_distance_definitions or []):
            if not isinstance(definition, dict):
                continue
            label = definition.get("label", f"Point_Plane_Distance_{index + 1}")
            plane_distance_values[label] = self.point_plane_distance(
                coords,
                definition.get("point_atoms", []),
                definition.get("plane_atoms", []),
                definition.get("point_symmetry") or None,
                definition.get("plane_symmetry") or None,
            )
            if mercury_output and self._point_uses_group_center(
                definition.get("point_atoms", [])
            ):
                mercury_plane_distance_values[label] = self.point_plane_distance(
                    coords,
                    definition.get("point_atoms", []),
                    definition.get("plane_atoms", []),
                    definition.get("point_symmetry") or None,
                    definition.get("plane_symmetry") or None,
                    group_center_round_dp=3,
                )

        self._append_measurements_csv(
            "point_geometry_distances.csv",
            structure_number,
            distance_values,
            results_path,
        )
        if mercury_output:
            self._append_measurements_csv(
                "point_geometry_distances_mercury.csv",
                structure_number,
                mercury_distance_values,
                results_path,
            )
        self._append_measurements_csv(
            "point_geometry_angles.csv", structure_number, angle_values, results_path
        )
        if mercury_output:
            self._append_measurements_csv(
                "point_geometry_angles_mercury.csv",
                structure_number,
                mercury_angle_values,
                results_path,
            )
        self._append_measurements_csv(
            "point_geometry_torsions.csv",
            structure_number,
            torsion_values,
            results_path,
        )
        if mercury_output:
            self._append_measurements_csv(
                "point_geometry_torsions_mercury.csv",
                structure_number,
                mercury_torsion_values,
                results_path,
            )
        self._append_measurements_csv(
            "point_geometry_plane_distances.csv",
            structure_number,
            plane_distance_values,
            results_path,
        )
        if mercury_output:
            self._append_measurements_csv(
                "point_geometry_plane_distances_mercury.csv",
                structure_number,
                mercury_plane_distance_values,
                results_path,
            )

    def point_group_distance(
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

        return self.point_distance(
            coords, atom_list_1, atom_list_2, symmetry_1, symmetry_2
        )

    def find_point_group_distance(
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
        return self.point_group_distance(
            coords, atom_list_1, atom_list_2, symmetry_1, symmetry_2
        )

    def analyse_point_group_distance(
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

        distance = self.find_point_group_distance(
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


