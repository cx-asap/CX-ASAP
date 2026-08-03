#!/usr/bin/env python3

import pathlib
import unittest

from CifFile import ReadCif
from post_refinement_analysis.modules.cif_geometry import CIF_Geometry

CIF = (
    pathlib.Path(__file__).parent.parent
    / "cx_asap"
    / "test_data"
    / "data"
    / "210K"
    / "210.cif"
)

SYM = "1-x, -y, 1-z"
# CIF_Geometry._extract_structure_data uppercases all atom labels
RING_ATOMS = ["CU1", "O1", "O2", "C2", "C3", "C4"]


class TestCIFGeometryEngine(unittest.TestCase):
    def setUp(self):
        cif = ReadCif(str(CIF))
        block = cif[cif.keys()[0]]

        self.geo = CIF_Geometry(test_mode=True)
        data = self.geo._extract_structure_data(block)
        self.coords = data["coords"]
        self.geo.point_engine.cell_params = data["cell"]
        self.geo.point_engine.bad_flag = False

    def test_distance_to_ring_center(self):
        """C1 to group center of chelate ring atoms (CIF coordinates)."""
        result = self.geo.point_engine.point_distance(self.coords, ["C1"], RING_ATOMS)
        self.assertAlmostEqual(result, 2.939, places=3)

    def test_angle_through_ring_center(self):
        """C1 - ring center - C5 angle (CIF coordinates)."""
        result = self.geo.point_engine.point_angle(
            self.coords, ["C1"], RING_ATOMS, ["C5"]
        )
        self.assertAlmostEqual(result, 116.001, places=3)

    def test_torsion_with_symmetry(self):
        """O1 - Cu1 - O1(sym) - O2(sym) torsion with 1-x,-y,1-z on points 3 and 4."""
        result = self.geo.point_engine.point_torsion(
            self.coords,
            ["O1"],
            ["CU1"],
            ["O1"],
            ["O2"],
            symmetry_3=SYM,
            symmetry_4=SYM,
        )
        self.assertAlmostEqual(result, -49.732, places=3)

    def test_plane_distance_equatorial(self):
        """H1C to equatorial plane defined by Cu1, O1, O2 (CIF coordinates)."""
        result = self.geo.point_engine.point_plane_distance(
            self.coords, ["H1C"], ["CU1", "O1", "O2"]
        )
        self.assertAlmostEqual(result, 1.176, places=3)

    def test_plane_distance_equatorial_symmetry_consistent(self):
        """H1C to equatorial plane with symmetry applied is identical (Cu1 maps to itself)."""
        no_sym = self.geo.point_engine.point_plane_distance(
            self.coords, ["H1C"], ["CU1", "O1", "O2"]
        )
        with_sym = self.geo.point_engine.point_plane_distance(
            self.coords, ["H1C"], ["CU1", "O1", "O2"], plane_symmetry=SYM
        )
        self.assertAlmostEqual(no_sym, with_sym, places=3)

    def test_plane_distance_chelate_ring(self):
        """H1C to best-fit plane through the full chelate ring (CIF coordinates)."""
        result = self.geo.point_engine.point_plane_distance(
            self.coords, ["H1C"], RING_ATOMS
        )
        self.assertAlmostEqual(result, 0.953, places=3)
