#!/usr/bin/env python3

###################################################################################################
# --------------------------------------CX-ASAP: points_pipeline-----------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Config, Directory_Browse, Grapher
from post_refinement_analysis.modules.points import PointGeometryEngine
import os
import pathlib
import pandas as pd
import logging

# ----------Class Definition----------#


class PointsPipeline:
    def __init__(self) -> None:
        """Initialises the class

        Sets up the yaml parameters input by the user

        Also defines the location for the system yaml file

        which stores a yaml of code-only parameters accessible throughout

        the software package
        """

        config = Config()

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    def create_numbered_results_directory(
        self, base_directory: str, folder_name: str
    ) -> pathlib.Path:
        """Creates a numbered results directory inside base_directory/folder_name."""

        base_path = pathlib.Path(base_directory)
        results_root = base_path / folder_name
        results_root.mkdir(exist_ok=True)

        existing_numbers = []
        for item in results_root.iterdir():
            if item.is_dir():
                try:
                    existing_numbers.append(int(item.name))
                except ValueError:
                    continue

        next_number = max(existing_numbers) + 1 if existing_numbers else 1
        results_path = results_root / str(next_number)
        results_path.mkdir()
        return results_path

    def point_group_distance_analysis(
        self,
        working_directory: str,
        atom_list_1: list,
        atom_list_2: list,
        results_directory: str,
        label: str = "Centroid Distance",
        symmetry_1: str = None,
        symmetry_2: str = None,
    ) -> None:
        """Calculates centroid-to-centroid distances for a series of .lst files
        in separate folders within a common parent folder.

        Produces a scatter graph of centroid distance vs structure number.

        Args:
            working_directory (str): full path to the parent folder containing
                                     folders with .lst files
            atom_list_1 (list): atom labels for the first point group
            atom_list_2 (list): atom labels for the second point group
            results_directory (str): full path to the output directory
            label (str): column header for the output csv and graph y-axis
            symmetry_1 (str): optional symmetry operation string for centroid 1
            symmetry_2 (str): optional symmetry operation string for centroid 2
        """

        logging.info(
            "Running centroid distance analysis with "
            f"working_directory={working_directory}, "
            f"results_directory={results_directory}, "
            f"atom_list_1={atom_list_1}, "
            f"atom_list_2={atom_list_2}, "
            f"symmetry_1={symmetry_1}, "
            f"symmetry_2={symmetry_2}"
        )

        centroid = PointGeometryEngine()
        tree = Directory_Browse(working_directory)

        for index, item in enumerate(tree.directories):
            tree.enter_directory(item, ".lst")
            centroid.analyse_point_group_distance(
                tree.item_file,
                index + 1,
                results_directory,
                atom_list_1,
                atom_list_2,
                label,
                symmetry_1,
                symmetry_2,
            )
            tree.exit_directory()

        os.chdir(results_directory)
        csv_name = label.lower().replace(" ", "_") + ".csv"

        try:
            full_data = pd.read_csv(csv_name)
            x = full_data["Structure"]
            y = full_data[label]
            graph = Grapher()
            graph.single_scatter_graph(
                x,
                y,
                "Structure Number",
                r"Distance ($\AA$)",
                label,
                csv_name.replace(".csv", ".png"),
            )
        except FileNotFoundError:
            logging.error(f"No {csv_name} file found...")

    def point_geometry_analysis(
        self,
        working_directory: str,
        results_directory: str,
        distance_definitions: "list[dict]" = None,
        angle_definitions: "list[dict]" = None,
        torsion_definitions: "list[dict]" = None,
        plane_distance_definitions: "list[dict]" = None,
        mercury_output: bool = False,
    ) -> None:
        """Calculates point-based geometry values across multiple .lst files.

        Supported metrics are:
                    - 2-point distances
          - 3-point angles
          - 4-point torsions
          - point-to-plane distances
        """

        logging.info(
            "Running point geometry analysis with "
            f"working_directory={working_directory}, "
            f"results_directory={results_directory}, "
            f"distance_definitions={len(distance_definitions or [])}, "
            f"angle_definitions={len(angle_definitions or [])}, "
            f"torsion_definitions={len(torsion_definitions or [])}, "
            f"plane_distance_definitions={len(plane_distance_definitions or [])}"
        )

        centroid = PointGeometryEngine()
        tree = Directory_Browse(working_directory)
        processed_structures = 0

        for index, item in enumerate(tree.directories):
            tree.enter_directory(item, ".lst")
            centroid.analyse_point_geometry(
                tree.item_file,
                index + 1,
                results_directory,
                distance_definitions,
                angle_definitions,
                torsion_definitions,
                plane_distance_definitions,
                mercury_output=mercury_output,
            )
            processed_structures += 1
            tree.exit_directory()

        os.chdir(results_directory)

        self._graph_point_geometry_csv(
            "point_geometry_distances.csv",
            "Point Geometry Distances",
            r"Distance ($\AA$)",
            "point_geometry_distances.png",
        )
        if mercury_output:
            self._graph_point_geometry_csv(
                "point_geometry_distances_mercury.csv",
                "Point Geometry Distances (Mercury)",
                r"Distance ($\AA$)",
                "point_geometry_distances_mercury.png",
            )
        self._graph_point_geometry_csv(
            "point_geometry_angles.csv",
            "Point Geometry Angles",
            "Angle($^\\circ$)",
            "point_geometry_angles.png",
        )
        if mercury_output:
            self._graph_point_geometry_csv(
                "point_geometry_angles_mercury.csv",
                "Point Geometry Angles (Mercury)",
                "Angle($^\\circ$)",
                "point_geometry_angles_mercury.png",
            )
        self._graph_point_geometry_csv(
            "point_geometry_torsions.csv",
            "Point Geometry Torsions",
            "Angle($^\\circ$)",
            "point_geometry_torsions.png",
        )
        if mercury_output:
            self._graph_point_geometry_csv(
                "point_geometry_torsions_mercury.csv",
                "Point Geometry Torsions (Mercury)",
                "Angle($^\\circ$)",
                "point_geometry_torsions_mercury.png",
            )
        self._graph_point_geometry_csv(
            "point_geometry_plane_distances.csv",
            "Point to Plane Distances",
            r"Distance ($\AA$)",
            "point_geometry_plane_distances.png",
        )
        if mercury_output:
            self._graph_point_geometry_csv(
                "point_geometry_plane_distances_mercury.csv",
                "Point to Plane Distances (Mercury)",
                r"Distance ($\AA$)",
                "point_geometry_plane_distances_mercury.png",
            )

        logging.info(
            "Point geometry analysis finished with "
            f"processed_structures={processed_structures}, "
            f"distances_csv={os.path.exists('point_geometry_distances.csv')}, "
            f"angles_csv={os.path.exists('point_geometry_angles.csv')}, "
            f"torsions_csv={os.path.exists('point_geometry_torsions.csv')}, "
            "plane_distances_csv="
            f"{os.path.exists('point_geometry_plane_distances.csv')}, "
            "mercury_output="
            f"{mercury_output}"
        )

    def _graph_point_geometry_csv(
        self, csv_name: str, graph_title: str, y_axis_title: str, figure_name: str
    ) -> None:
        """Creates a summary scatter graph for one point-geometry CSV file."""

        try:
            full_data = pd.read_csv(csv_name)
        except FileNotFoundError:
            return

        if "Structure" not in full_data.columns:
            return

        y_cols = [col for col in full_data.columns if col != "Structure"]
        if not y_cols:
            return

        x = full_data["Structure"]
        y_data = [list(full_data[col]) for col in y_cols]
        graph = Grapher()
        graph.single_scatter_graph(
            x,
            y_data,
            "Structure Number",
            y_axis_title,
            graph_title,
            figure_name,
            y_series_title=y_cols if len(y_cols) > 1 else None,
        )
