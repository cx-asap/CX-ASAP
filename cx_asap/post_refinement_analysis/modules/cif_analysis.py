#!/usr/bin/env python3

###################################################################################################
# --------------------------------------CX-ASAP: cif_analysis--------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

from system_files.utils import Config
from post_refinement_analysis.modules.cif_read import CIF_Read
from post_refinement_analysis.modules.structural_analysis import Structural_Analysis
from post_refinement_analysis.modules.ADP_analysis import ADP_analysis
from post_refinement_analysis.modules.cell_analysis import Cell_Deformation
from post_refinement_analysis.modules.cif_geometry import CIF_Geometry
from post_refinement_analysis.modules.rotation_planes import Rotation
import pathlib
import os
import logging


class CIF_Analysis:
    """Coordinates CIF parameter extraction and optional analysis modules."""

    def __init__(self, test_mode: bool = False) -> None:
        """Initialise shared configuration for one analysis run."""

        self.test_mode = test_mode

        config = Config(self.test_mode)
        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    @staticmethod
    def _find_first_file(location: str, suffix: str) -> str:
        """Return the first file matching suffix from a file or folder path."""

        base = pathlib.Path(location)

        if base.is_file():
            if base.suffix.lower() == suffix.lower():
                return str(base)
            return ""

        if not base.exists() or not base.is_dir():
            return ""

        for item in base.iterdir():
            if item.is_file() and item.suffix.lower() == suffix.lower():
                return str(item)
        return ""

    def run(
        self,
        cif_location: str,
        results_directory: str,
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
        lst_file_location: str = "",
    ) -> None:
        """Run CIF extraction plus structural/ADP/cell/geometry analysis outputs.

        The method always writes results into ``results_directory``. Geometry and
        rotation analyses are optional and are enabled by non-empty definitions.
        Legacy LST-based rotation fallback is used only when a reference plane is
        provided and no CIF-native rotation plane definitions are configured.
        """

        os.chdir(results_directory)

        point_geometry_distances = point_geometry_distances or []
        point_geometry_angles = point_geometry_angles or []
        point_geometry_torsions = point_geometry_torsions or []
        point_geometry_plane_distances = point_geometry_plane_distances or []
        rotation_plane_definitions = rotation_plane_definitions or []
        mean_plane_definitions = mean_plane_definitions or []

        cif_data = CIF_Read(self.test_mode)
        cif_data.configure(cif_parameters)
        cif_data.get_data(
            cif_location,
            structural_analysis_bonds,
            structural_analysis_angles,
            structural_analysis_torsions,
            structural_analysis_hbonds,
            ADP_analysis_enabled,
            varying_cif_parameter,
        )
        cif_data.data_output()

        geometry = Structural_Analysis(self.test_mode)
        bond_file = (
            "Bond_Lengths.csv"
            if structural_analysis_bonds and pathlib.Path("Bond_Lengths.csv").exists()
            else False
        )
        angle_file = (
            "Bond_Angles.csv"
            if structural_analysis_angles and pathlib.Path("Bond_Angles.csv").exists()
            else False
        )
        torsion_file = (
            "Bond_Torsions.csv"
            if structural_analysis_torsions
            and pathlib.Path("Bond_Torsions.csv").exists()
            else False
        )
        hbond_file = (
            "HBond_details.csv"
            if structural_analysis_hbonds and pathlib.Path("HBond_details.csv").exists()
            else False
        )

        geometry.import_and_analyse(
            bond_file,
            angle_file,
            torsion_file,
            hbond_file,
            atoms_for_analysis,
            results_directory,
            varying_parameter=varying_cif_parameter,
        )

        if reference_unit_cell not in ["", 0, None]:
            cell = Cell_Deformation(self.test_mode)
            cell.import_data("CIF_Parameters.csv", reference_unit_cell)
            cell.calculate_deformations()
            cell.quality_analysis(
                varying_cif_parameter,
                {
                    "R1": cell.df["_refine_ls_R_factor_gt"],
                    "Rint": cell.df["_diffrn_reflns_av_R_equivalents"],
                    "Completeness": cell.df["_diffrn_measured_fraction_theta_full"],
                },
                varying_cif_parameter,
            )
            cell.graphical_analysis(varying_cif_parameter, varying_cif_parameter)

        if ADP_analysis_enabled and pathlib.Path("ADPs.csv").exists():
            adp_obj = ADP_analysis(self.test_mode)
            adp_obj.analyse_data("ADPs.csv", "CIF_Parameters.csv")

        run_cif_geometry = (
            len(point_geometry_distances) > 0
            or len(point_geometry_angles) > 0
            or len(point_geometry_torsions) > 0
            or len(point_geometry_plane_distances) > 0
            or (
                reference_plane not in [None, 0, "", [0], [0, 0, 0]]
                and len(rotation_plane_definitions) > 0
            )
            or (calculate_interplane_angle and len(mean_plane_definitions) > 0)
        )

        if run_cif_geometry:
            cif_geometry = CIF_Geometry(self.test_mode)
            cif_geometry.run(
                cif_location,
                results_directory,
                varying_cif_parameter,
                point_geometry_distances,
                point_geometry_angles,
                point_geometry_torsions,
                point_geometry_plane_distances,
                mercury_output=mercury_output,
                reference_plane=reference_plane,
                rotation_plane_definitions=rotation_plane_definitions,
                mean_plane_definitions=mean_plane_definitions,
                calculate_interplane_angle=calculate_interplane_angle,
            )

        if lst_file_location in ["", 0, None]:
            lst_file_location = self._find_first_file(cif_location, ".lst")

        # Legacy fallback: retain old MPLA-on-LST behaviour when no CIF plane definitions are supplied.
        if (
            reference_plane not in [None, 0, "", [0], [0, 0, 0]]
            and len(rotation_plane_definitions) == 0
        ):
            if lst_file_location not in ["", 0, None]:
                rot = Rotation(self.test_mode)
                rot.configure(reference_plane)
                rot.analysis(lst_file_location, 1, results_directory)
                if calculate_interplane_angle:
                    rot.analyse_interplane_angle(lst_file_location, 1, results_directory)
            else:
                logging.warning(
                    __name__
                    + " : Rotation-plane analysis requested but no CIF rotation_plane_definitions or .lst file was found"
                )
