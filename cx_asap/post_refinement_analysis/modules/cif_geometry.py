#!/usr/bin/env python3

###################################################################################################
# --------------------------------------CX-ASAP: cif_geometry--------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

from CifFile import ReadCif
from post_refinement_analysis.modules.cif_read import CIF_Read
from post_refinement_analysis.modules.points import PointGeometryEngine
from system_files.crystal_math import (
    fractional_to_cartesian,
    best_fit_plane_normal,
    cartesian_plane_normal_to_fractional,
    angle_between_vectors,
)
import pathlib
import pandas as pd
import logging
import numpy as np


class CIF_Geometry:
    """Performs symmetry-aware point and plane analysis directly from CIF files."""

    def __init__(self, test_mode: bool = False) -> None:
        """Initialise CIF geometry engine and point-geometry helper."""

        self.test_mode = test_mode
        self.point_engine = PointGeometryEngine(self.test_mode)

    @staticmethod
    def _parse_float(raw) -> float:
        """Parses a CIF numeric value, handling values with uncertainties."""

        if isinstance(raw, (int, float)):
            return float(raw)

        raw_str = str(raw).strip()
        if raw_str in ["", ".", "?"]:
            raise ValueError("Missing numeric CIF value")

        if "(" in raw_str:
            raw_str = raw_str.split("(")[0]

        return float(raw_str)

    def _extract_structure_data(self, block) -> dict:
        """Extract unit-cell and fractional atom coordinates from one CIF block."""

        cell_keys = [
            "_cell_length_a",
            "_cell_length_b",
            "_cell_length_c",
            "_cell_angle_alpha",
            "_cell_angle_beta",
            "_cell_angle_gamma",
        ]
        cell = [self._parse_float(block[item]) for item in cell_keys]

        labels = block["_atom_site_label"]
        x_vals = block["_atom_site_fract_x"]
        y_vals = block["_atom_site_fract_y"]
        z_vals = block["_atom_site_fract_z"]

        if not isinstance(labels, list):
            labels = [labels]
            x_vals = [x_vals]
            y_vals = [y_vals]
            z_vals = [z_vals]

        coords = {}
        for label, x_raw, y_raw, z_raw in zip(labels, x_vals, y_vals, z_vals):
            try:
                coords[str(label).upper()] = [
                    self._parse_float(x_raw),
                    self._parse_float(y_raw),
                    self._parse_float(z_raw),
                ]
            except ValueError:
                continue

        return {
            "block": block,
            "coords": coords,
            "cell": cell,
        }

    @staticmethod
    def _append_rows(csv_name: str, rows: list, results_directory: str) -> None:
        """Write collected row dictionaries to a CSV when rows are present."""

        if len(rows) == 0:
            return

        df = pd.DataFrame(rows)
        output_path = pathlib.Path(results_directory) / csv_name
        df.to_csv(output_path, index=None)

    @staticmethod
    def _definition_uses_mercury_rounding(
        definition: dict, keys: list, engine: PointGeometryEngine
    ) -> bool:
        """Check whether any definition key should also be written in Mercury-rounded form."""

        for key in keys:
            if engine._point_uses_centroid(definition.get(key, [])):
                return True
        return False

    @staticmethod
    def _valid_reference_plane(reference_plane: list) -> bool:
        """Validate a reference plane is a numeric [h, k, l] list."""

        if not isinstance(reference_plane, list) or len(reference_plane) != 3:
            return False
        try:
            [float(item) for item in reference_plane]
        except (TypeError, ValueError):
            return False
        return True

    def _rotation_plane_angle(
        self,
        coords: dict,
        cell_params: list,
        definition: dict,
        reference_plane: list,
    ) -> tuple:
        """Compute acute angle between a fitted plane normal and reference plane."""

        normal_frac = self._plane_normal_fractional(coords, cell_params, definition)
        ref_frac = np.array(reference_plane, dtype=float)
        angle = angle_between_vectors(normal_frac, ref_frac, fold_to_acute=True)

        return angle, normal_frac

    def _plane_normal_fractional(
        self,
        coords: dict,
        cell_params: list,
        definition: dict,
    ) -> np.ndarray:
        """Fit a plane from atoms and return its normal in fractional coordinates."""

        atoms = definition.get("plane_atoms", [])
        symmetry = definition.get("plane_symmetry") or None
        frac_points = self.point_engine._extract_positions(coords, atoms, symmetry)

        if len(frac_points) < 3:
            raise ValueError("Need at least 3 atoms to define rotation plane")

        cart_points = np.array(
            [fractional_to_cartesian(item, cell_params) for item in frac_points],
            dtype=float,
        )
        normal_cart = best_fit_plane_normal(cart_points)
        normal_frac = cartesian_plane_normal_to_fractional(normal_cart, cell_params)
        return normal_frac

    def run(
        self,
        cif_location: str,
        results_directory: str,
        varying_parameter: str,
        point_geometry_distances: list = None,
        point_geometry_angles: list = None,
        point_geometry_torsions: list = None,
        point_geometry_plane_distances: list = None,
        mercury_output: bool = False,
        reference_plane: list = None,
        rotation_plane_definitions: list = None,
        mean_plane_definitions: list = None,
        calculate_interplane_angle: bool = False,
    ) -> None:
        """Run CIF-native point geometry and rotation-plane analyses.

        For each CIF block, this writes CSV outputs for any configured distance,
        angle, torsion, point-plane, rotation-plane, and optional interplane
        calculations. Mercury-style rounded outputs are generated when
        ``mercury_output`` is enabled and definitions require companion rounded values.
        """

        point_geometry_distances = point_geometry_distances or []
        point_geometry_angles = point_geometry_angles or []
        point_geometry_torsions = point_geometry_torsions or []
        point_geometry_plane_distances = point_geometry_plane_distances or []
        rotation_plane_definitions = rotation_plane_definitions or []
        mean_plane_definitions = mean_plane_definitions or []

        files = CIF_Read._read_cif_files(cif_location)
        if len(files) == 0:
            logging.warning(__name__ + " : No .cif files found for CIF geometry analysis")
            return

        rows_distance = []
        rows_angle = []
        rows_torsion = []
        rows_plane = []
        rows_distance_mercury = []
        rows_angle_mercury = []
        rows_torsion_mercury = []
        rows_plane_mercury = []
        rows_rotation = []
        rows_interplane = []

        structure_number = 0
        for cif_file in files:
            cif_obj = ReadCif(str(cif_file))
            for block_name in cif_obj.keys():
                structure_number += 1

                try:
                    structure = self._extract_structure_data(cif_obj[block_name])
                except Exception as error:
                    logging.warning(
                        __name__
                        + " : Failed to extract structure from "
                        + str(cif_file)
                        + " block "
                        + str(block_name)
                        + " due to "
                        + str(error)
                    )
                    continue

                block = structure["block"]
                coords = structure["coords"]
                cell = structure["cell"]
                self.point_engine.cell_params = cell

                base_row = {
                    "Structure": structure_number,
                    "CIF_File": cif_file.stem,
                    "Data_Block": block_name,
                }
                try:
                    base_row[varying_parameter] = self._parse_float(block[varying_parameter])
                except Exception:
                    pass

                distance_row = dict(base_row)
                angle_row = dict(base_row)
                torsion_row = dict(base_row)
                plane_row = dict(base_row)
                distance_row_mercury = dict(base_row)
                angle_row_mercury = dict(base_row)
                torsion_row_mercury = dict(base_row)
                plane_row_mercury = dict(base_row)

                for index, definition in enumerate(point_geometry_distances):
                    if not isinstance(definition, dict):
                        continue
                    label = definition.get("label", f"Distance_{index + 1}")
                    distance_row[label] = self.point_engine.point_distance(
                        coords,
                        definition.get("point_1_atoms", []),
                        definition.get("point_2_atoms", []),
                        definition.get("point_1_symmetry") or None,
                        definition.get("point_2_symmetry") or None,
                    )
                    if mercury_output and self._definition_uses_mercury_rounding(
                        definition,
                        ["point_1_atoms", "point_2_atoms"],
                        self.point_engine,
                    ):
                        distance_row_mercury[label] = self.point_engine.point_distance(
                            coords,
                            definition.get("point_1_atoms", []),
                            definition.get("point_2_atoms", []),
                            definition.get("point_1_symmetry") or None,
                            definition.get("point_2_symmetry") or None,
                            centroid_round_dp=3,
                        )

                for index, definition in enumerate(point_geometry_angles):
                    if not isinstance(definition, dict):
                        continue
                    label = definition.get("label", f"Angle_{index + 1}")
                    angle_row[label] = self.point_engine.point_angle(
                        coords,
                        definition.get("point_1_atoms", []),
                        definition.get("point_2_atoms", []),
                        definition.get("point_3_atoms", []),
                        definition.get("point_1_symmetry") or None,
                        definition.get("point_2_symmetry") or None,
                        definition.get("point_3_symmetry") or None,
                    )
                    if mercury_output and self._definition_uses_mercury_rounding(
                        definition,
                        ["point_1_atoms", "point_2_atoms", "point_3_atoms"],
                        self.point_engine,
                    ):
                        angle_row_mercury[label] = self.point_engine.point_angle(
                            coords,
                            definition.get("point_1_atoms", []),
                            definition.get("point_2_atoms", []),
                            definition.get("point_3_atoms", []),
                            definition.get("point_1_symmetry") or None,
                            definition.get("point_2_symmetry") or None,
                            definition.get("point_3_symmetry") or None,
                            centroid_round_dp=3,
                        )

                for index, definition in enumerate(point_geometry_torsions):
                    if not isinstance(definition, dict):
                        continue
                    label = definition.get("label", f"Torsion_{index + 1}")
                    torsion_row[label] = self.point_engine.point_torsion(
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
                    if mercury_output and self._definition_uses_mercury_rounding(
                        definition,
                        [
                            "point_1_atoms",
                            "point_2_atoms",
                            "point_3_atoms",
                            "point_4_atoms",
                        ],
                        self.point_engine,
                    ):
                        torsion_row_mercury[label] = self.point_engine.point_torsion(
                            coords,
                            definition.get("point_1_atoms", []),
                            definition.get("point_2_atoms", []),
                            definition.get("point_3_atoms", []),
                            definition.get("point_4_atoms", []),
                            definition.get("point_1_symmetry") or None,
                            definition.get("point_2_symmetry") or None,
                            definition.get("point_3_symmetry") or None,
                            definition.get("point_4_symmetry") or None,
                            centroid_round_dp=3,
                        )

                for index, definition in enumerate(point_geometry_plane_distances):
                    if not isinstance(definition, dict):
                        continue
                    label = definition.get("label", f"Point_Plane_Distance_{index + 1}")
                    plane_row[label] = self.point_engine.point_plane_distance(
                        coords,
                        definition.get("point_atoms", []),
                        definition.get("plane_atoms", []),
                        definition.get("point_symmetry") or None,
                        definition.get("plane_symmetry") or None,
                    )
                    if mercury_output and self._definition_uses_mercury_rounding(
                        definition,
                        ["point_atoms"],
                        self.point_engine,
                    ):
                        plane_row_mercury[label] = self.point_engine.point_plane_distance(
                            coords,
                            definition.get("point_atoms", []),
                            definition.get("plane_atoms", []),
                            definition.get("point_symmetry") or None,
                            definition.get("plane_symmetry") or None,
                            centroid_round_dp=3,
                        )

                if len(distance_row.keys()) > len(base_row.keys()):
                    rows_distance.append(distance_row)
                if len(angle_row.keys()) > len(base_row.keys()):
                    rows_angle.append(angle_row)
                if len(torsion_row.keys()) > len(base_row.keys()):
                    rows_torsion.append(torsion_row)
                if len(plane_row.keys()) > len(base_row.keys()):
                    rows_plane.append(plane_row)

                if mercury_output:
                    if len(distance_row_mercury.keys()) > len(base_row.keys()):
                        rows_distance_mercury.append(distance_row_mercury)
                    if len(angle_row_mercury.keys()) > len(base_row.keys()):
                        rows_angle_mercury.append(angle_row_mercury)
                    if len(torsion_row_mercury.keys()) > len(base_row.keys()):
                        rows_torsion_mercury.append(torsion_row_mercury)
                    if len(plane_row_mercury.keys()) > len(base_row.keys()):
                        rows_plane_mercury.append(plane_row_mercury)

                if self._valid_reference_plane(reference_plane):
                    rotation_row = dict(base_row)
                    for index, definition in enumerate(rotation_plane_definitions):
                        if not isinstance(definition, dict):
                            continue

                        label = definition.get("label", f"MPLA_{index + 1}_Rotation_Angle")
                        try:
                            angle, normal_frac = self._rotation_plane_angle(
                                coords,
                                cell,
                                definition,
                                reference_plane,
                            )
                        except ValueError as error:
                            logging.warning(__name__ + " : " + str(error))
                            continue

                        rotation_row[label] = angle

                    if len(rotation_row.keys()) > len(base_row.keys()):
                        rows_rotation.append(rotation_row)

                if calculate_interplane_angle:
                    mean_normals = []
                    mean_labels = []
                    for index, definition in enumerate(mean_plane_definitions):
                        if not isinstance(definition, dict):
                            continue

                        try:
                            mean_normals.append(
                                self._plane_normal_fractional(coords, cell, definition)
                            )
                            mean_labels.append(
                                definition.get("label", f"Mean_Plane_{index + 1}")
                            )
                        except ValueError as error:
                            logging.warning(__name__ + " : " + str(error))
                            continue

                    if len(mean_normals) >= 2:
                        interplane = angle_between_vectors(
                            mean_normals[0], mean_normals[1], fold_to_acute=True
                        )
                        rows_interplane.append(
                            {
                                "Structure": structure_number,
                                "CIF_File": cif_file.stem,
                                "Data_Block": block_name,
                                "Mean Plane 1": mean_labels[0],
                                "Mean Plane 2": mean_labels[1],
                                "Interplane Angle": interplane,
                            }
                        )
                    elif len(mean_plane_definitions) > 0:
                        logging.warning(
                            __name__
                            + " : calculate_interplane_angle requires at least two valid mean planes"
                        )

        self._append_rows("point_geometry_distances.csv", rows_distance, results_directory)
        self._append_rows("point_geometry_angles.csv", rows_angle, results_directory)
        self._append_rows("point_geometry_torsions.csv", rows_torsion, results_directory)
        self._append_rows(
            "point_geometry_plane_distances.csv", rows_plane, results_directory
        )

        if mercury_output:
            self._append_rows(
                "point_geometry_distances_mercury.csv",
                rows_distance_mercury,
                results_directory,
            )
            self._append_rows(
                "point_geometry_angles_mercury.csv",
                rows_angle_mercury,
                results_directory,
            )
            self._append_rows(
                "point_geometry_torsions_mercury.csv",
                rows_torsion_mercury,
                results_directory,
            )
            self._append_rows(
                "point_geometry_plane_distances_mercury.csv",
                rows_plane_mercury,
                results_directory,
            )

        self._append_rows("rotation_angles.csv", rows_rotation, results_directory)
        self._append_rows("interplane_angles.csv", rows_interplane, results_directory)
