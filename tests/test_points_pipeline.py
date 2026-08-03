#!/usr/bin/env python3

import os
import pathlib
import shutil
import tempfile
import unittest

import pandas as pd

from post_refinement_analysis.pipelines.points_pipeline import PointsPipeline

DATA_DIR = (
    pathlib.Path(__file__).parent.parent / "cx_asap" / "test_data" / "data"
)

_LST_AVAILABLE = (DATA_DIR / "210K" / "210.lst").exists()

DIST_DEFS = [
    {
        "label": "Ring_Center_To_C1",
        "point_1_atoms": ["C1"],
        "point_2_atoms": ["Cu1", "O1", "O2", "C2", "C3", "C4"],
    }
]
ANGLE_DEFS = [
    {
        "label": "C1_Center_C5",
        "point_1_atoms": ["C1"],
        "point_2_atoms": ["Cu1", "O1", "O2", "C2", "C3", "C4"],
        "point_3_atoms": ["C5"],
    }
]


@unittest.skipIf(not _LST_AVAILABLE, "test data not available")
class TestPointsPipeline(unittest.TestCase):
    def setUp(self):
        self._orig_cwd = os.getcwd()
        self._tmp = tempfile.mkdtemp()

    def tearDown(self):
        os.chdir(self._orig_cwd)
        shutil.rmtree(self._tmp, ignore_errors=True)

    def _run(self):
        pipe = PointsPipeline()
        pipe.point_geometry_analysis(
            str(DATA_DIR),
            self._tmp,
            distance_definitions=DIST_DEFS,
            angle_definitions=ANGLE_DEFS,
        )
        csv_d = pathlib.Path(self._tmp) / "point_geometry_distances.csv"
        csv_a = pathlib.Path(self._tmp) / "point_geometry_angles.csv"
        return pd.read_csv(csv_d), pd.read_csv(csv_a)

    def test_output_files_exist(self):
        """Both geometry CSVs are created."""
        self._run()
        self.assertTrue((pathlib.Path(self._tmp) / "point_geometry_distances.csv").exists())
        self.assertTrue((pathlib.Path(self._tmp) / "point_geometry_angles.csv").exists())

    def test_distance_csv_columns(self):
        """Distance CSV has Structure and label columns."""
        df_d, _ = self._run()
        self.assertIn("Structure", df_d.columns)
        self.assertIn("Ring_Center_To_C1", df_d.columns)

    def test_angle_csv_columns(self):
        """Angle CSV has Structure and label columns."""
        _, df_a = self._run()
        self.assertIn("Structure", df_a.columns)
        self.assertIn("C1_Center_C5", df_a.columns)

    def test_row_count(self):
        """One row per dataset folder containing a .lst file (5 datasets)."""
        df_d, df_a = self._run()
        self.assertEqual(len(df_d), 5)
        self.assertEqual(len(df_a), 5)

    def test_distance_values_in_range(self):
        """All ring-center distances are in physically reasonable range."""
        df_d, _ = self._run()
        for val in df_d["Ring_Center_To_C1"]:
            self.assertGreater(val, 2.9)
            self.assertLess(val, 3.0)

    def test_angle_values_in_range(self):
        """All C1-center-C5 angles are in physically reasonable range."""
        _, df_a = self._run()
        for val in df_a["C1_Center_C5"]:
            self.assertGreater(val, 115.0)
            self.assertLess(val, 117.0)

    def test_210K_distance_golden(self):
        """210K dataset (structure 2 in sorted order) matches golden distance value."""
        df_d, _ = self._run()
        # sorted order: 200K=1, 210K=2, 220K=3, 230K=4, 240K=5
        val_210k = df_d[df_d["Structure"] == 2]["Ring_Center_To_C1"].iloc[0]
        self.assertAlmostEqual(val_210k, 2.939, places=3)

    def test_210K_angle_golden(self):
        """210K dataset (structure 2 in sorted order) matches golden angle value."""
        _, df_a = self._run()
        val_210k = df_a[df_a["Structure"] == 2]["C1_Center_C5"].iloc[0]
        self.assertAlmostEqual(val_210k, 115.995, places=3)
