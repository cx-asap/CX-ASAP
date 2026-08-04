#!/usr/bin/env python3

import os
import pathlib
import shutil
import tempfile
import unittest

import pandas as pd

from post_refinement_analysis.pipelines.cif_analysis_pipeline import CIF_Analysis_Pipeline

DATA_DIR = pathlib.Path(__file__).parent.parent / "cx_asap" / "test_data" / "data"
GOLDEN = DATA_DIR / "CIF_Analysis" / "1"

_DATA_AVAILABLE = (DATA_DIR / "210K" / "210.cif").exists() and GOLDEN.exists()

CIF_PARAMS = [
    "_cell_length_a",
    "_cell_length_b",
    "_cell_length_c",
    "_cell_angle_alpha",
    "_cell_angle_beta",
    "_cell_angle_gamma",
    "_cell_volume",
    "_diffrn_reflns_av_R_equivalents",
    "_diffrn_measured_fraction_theta_full",
    "_diffrn_ambient_temperature",
    "_refine_ls_R_factor_gt",
]


@unittest.skipIf(not _DATA_AVAILABLE, "test data not available")
class TestCIFAnalysisPipeline(unittest.TestCase):
    def setUp(self):
        self._orig_cwd = os.getcwd()
        self._tmp = tempfile.mkdtemp()
        pipe = CIF_Analysis_Pipeline()
        pipe.run(
            str(DATA_DIR),
            self._tmp,
            "nested",
            CIF_PARAMS,
            ["Cu1"],
            "_diffrn_ambient_temperature",
            ADP_analysis_enabled=True,
        )

    def tearDown(self):
        os.chdir(self._orig_cwd)
        shutil.rmtree(self._tmp, ignore_errors=True)

    def test_cif_parameters_exists(self):
        """CIF_Parameters.csv is produced."""
        self.assertTrue((pathlib.Path(self._tmp) / "CIF_Parameters.csv").exists())

    def test_adp_csv_exists(self):
        """ADPs.csv is produced when ADP_analysis_enabled=True."""
        self.assertTrue((pathlib.Path(self._tmp) / "ADPs.csv").exists())

    def test_cif_parameters_row_count(self):
        """One row per temperature dataset (5 total)."""
        df = pd.read_csv(pathlib.Path(self._tmp) / "CIF_Parameters.csv")
        self.assertEqual(len(df), 5)

    def test_cif_parameters_temperatures(self):
        """Temperatures match the five VT datasets."""
        df = pd.read_csv(pathlib.Path(self._tmp) / "CIF_Parameters.csv")
        temps = sorted(df["_diffrn_ambient_temperature"].tolist())
        self.assertEqual(temps, [200.0, 210.0, 220.0, 230.0, 240.0])

    def test_cif_parameters_a_axis_matches_golden(self):
        """a-axis values match the stored CIF_Analysis/1 golden output."""
        produced = pd.read_csv(pathlib.Path(self._tmp) / "CIF_Parameters.csv")
        golden = pd.read_csv(GOLDEN / "CIF_Parameters.csv")
        produced_sorted = produced.sort_values("_diffrn_ambient_temperature").reset_index(drop=True)
        golden_sorted = golden.sort_values("_diffrn_ambient_temperature").reset_index(drop=True)
        for a_prod, a_gold in zip(
            produced_sorted["_cell_length_a"], golden_sorted["_cell_length_a"]
        ):
            self.assertAlmostEqual(a_prod, a_gold, places=4)

    def test_adp_row_count_matches_golden(self):
        """ADP row count matches the stored golden output."""
        produced = pd.read_csv(pathlib.Path(self._tmp) / "ADPs.csv")
        golden = pd.read_csv(GOLDEN / "ADPs.csv")
        self.assertEqual(len(produced), len(golden))
