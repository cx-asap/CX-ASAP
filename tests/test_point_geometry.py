#!/usr/bin/env python3

import pathlib
import unittest

from post_refinement_analysis.modules.points import PointGeometryEngine

LST = (
    pathlib.Path(__file__).parent.parent
    / "cx_asap"
    / "test_data"
    / "data"
    / "210K"
    / "210.lst"
)

SYM = "1-x, -y, 1-z"
RING_ATOMS = ["Cu1", "O1", "O2", "C2", "C3", "C4"]


@unittest.skipIf(not LST.exists(), "test data not available")
class TestPointGeometryEngine(unittest.TestCase):
    def setUp(self):
        self.engine = PointGeometryEngine(test_mode=True)
        self.engine.grab_cell(str(LST))
        data = self.engine.lst_reader.read(str(LST))
        self.coords = self.engine.lst_reader.extract_atom_coordinates(data)

    def test_distance_to_ring_center(self):
        """C1 to group center of chelate ring atoms."""
        result = self.engine.point_distance(self.coords, ["C1"], RING_ATOMS)
        self.assertAlmostEqual(result, 2.939, places=3)

    def test_angle_through_ring_center(self):
        """C1 - ring center - C5 angle."""
        result = self.engine.point_angle(self.coords, ["C1"], RING_ATOMS, ["C5"])
        self.assertAlmostEqual(result, 115.995, places=3)

    def test_torsion_with_symmetry(self):
        """O1 - Cu1 - O1(sym) - O2(sym) torsion with 1-x,-y,1-z on points 3 and 4."""
        result = self.engine.point_torsion(
            self.coords,
            ["O1"],
            ["Cu1"],
            ["O1"],
            ["O2"],
            symmetry_3=SYM,
            symmetry_4=SYM,
        )
        self.assertAlmostEqual(result, -48.913, places=3)

    def test_plane_distance_equatorial(self):
        """H1C to equatorial plane defined by Cu1, O1, O2."""
        result = self.engine.point_plane_distance(
            self.coords, ["H1C"], ["Cu1", "O1", "O2"]
        )
        self.assertAlmostEqual(result, 1.176, places=3)

    def test_plane_distance_equatorial_symmetry_consistent(self):
        """H1C to equatorial plane with symmetry applied is identical (Cu1 maps to itself)."""
        no_sym = self.engine.point_plane_distance(
            self.coords, ["H1C"], ["Cu1", "O1", "O2"]
        )
        with_sym = self.engine.point_plane_distance(
            self.coords, ["H1C"], ["Cu1", "O1", "O2"], plane_symmetry=SYM
        )
        self.assertAlmostEqual(no_sym, with_sym, places=3)

    def test_plane_distance_chelate_ring(self):
        """H1C to best-fit plane through the full chelate ring."""
        result = self.engine.point_plane_distance(
            self.coords, ["H1C"], RING_ATOMS
        )
        self.assertAlmostEqual(result, 0.953, places=3)
