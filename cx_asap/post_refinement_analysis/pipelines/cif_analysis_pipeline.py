#!/usr/bin/env python3

###################################################################################################
# ----------------------------------CX-ASAP: cif_analysis_pipeline---------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

from system_files.utils import Config
from system_files.utils import Grapher
from post_refinement_analysis.modules.cif_analysis import CIF_Analysis
from CifFile import ReadCif
import pathlib
import logging
import pandas as pd


class CIF_Analysis_Pipeline:
    """Pipeline wrapper for running CIF analysis on directories or files."""

    def __init__(self, test_mode: bool = False) -> None:
        """Initialise pipeline configuration and runtime options."""

        self.test_mode = test_mode

        config = Config(self.test_mode)
        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    def create_numbered_results_directory(
        self, base_directory: str, folder_name: str
    ) -> pathlib.Path:
        """Create the next numbered results folder under ``base_directory``."""

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

    def _dataset_directories(self, experiment_location: str, cif_input_mode: str) -> list:
        """Return dataset folders for flat or nested experiment layouts."""

        mode = str(cif_input_mode).strip().lower()
        exp = pathlib.Path(experiment_location)

        if mode == "flat":
            return [exp]

        nested_dirs = [item for item in exp.iterdir() if item.is_dir()]
        if len(nested_dirs) == 0:
            return [exp]

        return nested_dirs

    @staticmethod
    def _is_results_folder(folder_name: str) -> bool:
        """Identify folders that should be excluded from CIF discovery."""

        name = folder_name.strip().lower()
        blocked_exact = {
            "cif_analysis",
            "geometry_analysis",
            "refinement_statistics",
            "results",
            "analysis",
            "ref",
            "failed_autoprocessing",
        }
        if name in blocked_exact:
            return True
        if name.startswith("_"):
            return True
        return False

    def _discover_cif_files(self, location: str) -> list:
        """Root-first CIF discovery with one-level fallback."""

        base = pathlib.Path(location)

        if base.is_file():
            if base.suffix.lower() == ".cif":
                return [base.resolve()]
            return []

        if not base.exists() or not base.is_dir():
            return []

        root_files = [item for item in sorted(base.glob("*.cif")) if item.is_file()]

        nested_files = []
        for child in sorted(base.iterdir()):
            if not child.is_dir():
                continue
            if self._is_results_folder(child.name):
                continue
            nested_files.extend(
                [item for item in sorted(child.glob("*.cif")) if item.is_file()]
            )

        if len(root_files) > 0:
            if len(nested_files) > 0:
                logging.warning(
                    __name__
                    + " : Mixed CIF layout detected (root-level and nested CIFs). "
                    + "Using root-level CIFs only; nested CIFs will be ignored."
                )
            return [item.resolve() for item in root_files]

        seen = set()
        unique_files = []
        for item in nested_files:
            key = str(item.resolve())
            if key in seen:
                continue
            seen.add(key)
            unique_files.append(item)

        return unique_files

    @staticmethod
    def _write_combined_cif(cif_files: list, combined_path: pathlib.Path) -> None:
        """Merge multiple CIF files into one combined CIF file."""

        with open(combined_path, "w") as out_file:
            for cif_file in cif_files:
                cif_obj = ReadCif(str(cif_file))
                out_file.write(cif_obj.WriteOut())

    @staticmethod
    def _write_combined_manifest(cif_files: list, manifest_path: pathlib.Path) -> None:
        """Write source-file provenance list for a combined CIF input."""

        with open(manifest_path, "w") as manifest:
            for cif_file in cif_files:
                manifest.write(str(cif_file) + "\n")

    @staticmethod
    def _combine_dataset_csvs(results_directory: str, csv_name: str) -> "pd.DataFrame | None":
        """Combines one CSV from all dataset_* folders into a single pipeline-level CSV."""

        root = pathlib.Path(results_directory)
        root_csv = root / csv_name

        # If analysis already wrote directly to the root folder, reuse that file.
        if root_csv.exists():
            try:
                return pd.read_csv(root_csv)
            except Exception as error:
                logging.warning(
                    __name__
                    + " : Could not read "
                    + str(root_csv)
                    + " due to "
                    + str(error)
                )

        frames = []

        for dataset_dir in sorted(root.glob("dataset_*")):
            csv_path = dataset_dir / csv_name
            if not csv_path.exists():
                continue

            try:
                df = pd.read_csv(csv_path)
            except Exception as error:
                logging.warning(
                    __name__
                    + " : Could not read "
                    + str(csv_path)
                    + " due to "
                    + str(error)
                )
                continue

            if len(df) == 0:
                continue

            df["Dataset"] = dataset_dir.name
            frames.append(df)

        if len(frames) == 0:
            return None

        combined = pd.concat(frames, ignore_index=True)
        combined.to_csv(root / csv_name, index=None)
        return combined

    @staticmethod
    def _plot_combined_csv(
        df: "pd.DataFrame | None",
        varying_cif_parameter: str,
        y_axis_title: str,
        graph_title: str,
        figure_name: str,
        results_directory: str,
    ) -> None:
        """Creates a points-style scatter plot for a combined CSV output."""

        if df is None or len(df) == 0:
            return

        if varying_cif_parameter in df.columns:
            x_col = varying_cif_parameter
            x_title = varying_cif_parameter
        elif "Structure" in df.columns:
            x_col = "Structure"
            x_title = "Structure Number"
        else:
            return

        metadata_cols = {"Structure", "CIF_File", "Data_Block", "Dataset", x_col}
        y_cols = [col for col in df.columns if col not in metadata_cols]
        y_cols = [col for col in y_cols if pd.api.types.is_numeric_dtype(df[col])]

        if not y_cols:
            return

        x_data = list(df[x_col])
        y_data = [list(df[col]) for col in y_cols]

        graph = Grapher()
        graph.single_scatter_graph(
            x_data,
            y_data,
            x_title,
            y_axis_title,
            graph_title,
            str(pathlib.Path(results_directory) / figure_name),
            y_series_title=y_cols if len(y_cols) > 1 else None,
        )

    def run(
        self,
        experiment_location: str,
        results_directory: str,
        cif_input_mode: str,
        cif_parameters: list,
        atoms_for_analysis: list,
        varying_cif_parameter: str,
        reference_unit_cell: str = "",
        structural_analysis_bonds: bool = False,
        structural_analysis_angles: bool = False,
        structural_analysis_torsions: bool = False,
        structural_analysis_hbonds: bool = False,
        ADP_analysis_enabled: bool = False,
        point_geometry_distances: list = None,
        point_geometry_angles: list = None,
        point_geometry_torsions: list = None,
        point_geometry_plane_distances: list = None,
        mercury_output: bool = False,
        reference_plane: list = None,
        rotation_plane_definitions: list = None,
        mean_plane_definitions: list = None,
        calculate_interplane_angle: bool = False,
        precombine_cifs: bool = False,
    ) -> None:
        """Execute a CIF analysis pipeline run for a chosen experiment location.

        When ``precombine_cifs`` is true and multiple inputs are discovered, a
        combined CIF and source manifest are created in the run output folder,
        then analysed as a single input.
        """

        analyser = CIF_Analysis(self.test_mode)

        logging.info(
            __name__
            + " : Running CIF analysis pipeline with mode="
            + str(cif_input_mode)
            + " at root location="
            + str(experiment_location)
        )

        analysis_location = str(experiment_location)

        if precombine_cifs:
            cif_files = self._discover_cif_files(str(experiment_location))
            if len(cif_files) == 0:
                logging.warning(
                    __name__
                    + " : precombine_cifs is true but no CIF files were discovered at "
                    + str(experiment_location)
                )
            elif len(cif_files) == 1:
                analysis_location = str(cif_files[0])
                logging.info(
                    __name__
                    + " : precombine_cifs requested with one CIF; analysing file directly: "
                    + analysis_location
                )
            else:
                combined_path = pathlib.Path(results_directory) / "combined_input.cif"
                manifest_path = pathlib.Path(results_directory) / "combined_input_sources.txt"
                self._write_combined_cif(cif_files, combined_path)
                self._write_combined_manifest(cif_files, manifest_path)
                analysis_location = str(combined_path)
                logging.info(
                    __name__
                    + " : precombine_cifs created combined input: "
                    + str(combined_path)
                )

        # Run once at the selected analysis location to keep output layout
        # aligned with regular CIF analysis (single results folder).
        analyser.run(
            analysis_location,
            str(results_directory),
            cif_parameters,
            atoms_for_analysis,
            varying_cif_parameter,
            reference_unit_cell=reference_unit_cell,
            structural_analysis_bonds=structural_analysis_bonds,
            structural_analysis_angles=structural_analysis_angles,
            structural_analysis_torsions=structural_analysis_torsions,
            structural_analysis_hbonds=structural_analysis_hbonds,
            ADP_analysis_enabled=ADP_analysis_enabled,
            point_geometry_distances=point_geometry_distances,
            point_geometry_angles=point_geometry_angles,
            point_geometry_torsions=point_geometry_torsions,
            point_geometry_plane_distances=point_geometry_plane_distances,
            mercury_output=mercury_output,
            reference_plane=reference_plane,
            rotation_plane_definitions=rotation_plane_definitions,
            mean_plane_definitions=mean_plane_definitions,
            calculate_interplane_angle=calculate_interplane_angle,
        )

        # Consolidate angle outputs at pipeline root to match points-pipeline UX.
        point_angles_df = self._combine_dataset_csvs(
            results_directory, "point_geometry_angles.csv"
        )
        self._plot_combined_csv(
            point_angles_df,
            varying_cif_parameter,
            "Angle($^\\circ$)",
            "Point Geometry Angles",
            "point_geometry_angles.png",
            results_directory,
        )

        rotation_angles_df = self._combine_dataset_csvs(
            results_directory, "rotation_angles.csv"
        )
        self._plot_combined_csv(
            rotation_angles_df,
            varying_cif_parameter,
            "Angle($^\\circ$)",
            "Rotation Angles",
            "rotation_angles.png",
            results_directory,
        )

        interplane_angles_df = self._combine_dataset_csvs(
            results_directory, "interplane_angles.csv"
        )
        self._plot_combined_csv(
            interplane_angles_df,
            varying_cif_parameter,
            "Angle($^\\circ$)",
            "Interplane Angles",
            "interplane_angles.png",
            results_directory,
        )
