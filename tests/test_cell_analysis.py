#!/usr/bin/env python

import unittest

import pandas as pd
from post_refinement_analysis.modules.cell_analysis import Cell_Deformation


class testCellDeformation(unittest.TestCase):
    def setUp(self):
        """
        Defines a reference cell and unit cell data to calculate the deformations
        """

        self.test = Cell_Deformation(test_mode=True)

        self.ref_cell = [
            10.272757,
            4.67377,
            11.30969,
            90.0,
            92.1849,
            90.0,
            542.6117640255889,
        ]

        data = {
            "_cell_length_a": [10.2728, 10.2751, 10.2746, 10.2812, 10.2773],
            "_cell_length_b": [4.6738, 4.6771, 4.6799, 4.6838, 4.6861],
            "_cell_length_c": [11.3097, 11.3137, 11.3211, 11.3208, 11.3335],
            "_cell_angle_alpha": [90, 90, 90, 90, 90],
            "_cell_angle_beta": [92.185, 92.164, 92.114, 92.118, 92.018],
            "_cell_angle_gamma": [90, 90, 90, 90, 90],
            "_cell_volume": [542.61, 543.33, 544, 544.79, 545.48],
        }

        self.df = pd.DataFrame(data)

    def test_deformation_calculation(self):
        """
        Checks that the output matches the values of the calcualted deformations

        Rounds to 6 decimal places
        """

        expected_a = [0.000419, 0.022808, 0.017941, 0.082188, 0.044224]
        expected_b = [0.000642, 0.071249, 0.131158, 0.214602, 0.263813]
        expected_c = [0.000088, 0.035456, 0.100887, 0.098234, 0.210527]
        expected_alpha = [0, 0, 0, 0, 0]
        expected_beta = [0.000108, -0.022672, -0.076911, -0.072572, -0.181049]
        expected_gamma = [0, 0, 0, 0, 0]
        expected_volume = [-0.000325, 0.132366, 0.255843, 0.401435, 0.528598]

        output_deformation = self.test.deformation(self.df, self.ref_cell)

        output_a = [round(item, 6) for item in output_deformation["a_axis_deformation"]]
        output_b = [round(item, 6) for item in output_deformation["b_axis_deformation"]]
        output_c = [round(item, 6) for item in output_deformation["c_axis_deformation"]]
        output_alpha = [round(item, 6) for item in output_deformation["alpha_deformation"]]
        output_beta = [round(item, 6) for item in output_deformation["beta_deformation"]]
        output_gamma = [round(item, 6) for item in output_deformation["gamma_deformation"]]
        output_vol = [round(item, 6) for item in output_deformation["volume_deformation"]]

        self.assertEqual(output_a, expected_a)
        self.assertEqual(output_b, expected_b)
        self.assertEqual(output_c, expected_c)
        self.assertEqual(output_alpha, expected_alpha)
        self.assertEqual(output_beta, expected_beta)
        self.assertEqual(output_gamma, expected_gamma)
        self.assertEqual(output_vol, expected_volume)
