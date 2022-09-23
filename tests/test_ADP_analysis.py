#!/usr/bin/env python

import unittest

import numpy as np
from post_refinement_analysis.modules.ADP_analysis import ADP_analysis


class testADPanalysis(unittest.TestCase):
    def setUp(self):
        """
        Defines test cases for testing the ADP module
        """

        self.test = ADP_analysis(test_mode=True)

        self.vector = np.array([[-0.419205, 0.72285586, -0.54931464]])
        self.G = np.array(
            [
                [1.05530420e02, 2.93994911e-15, -4.42958924e00],
                [2.93994911e-15, 2.18444064e01, 3.23669715e-15],
                [-4.42958924e00, 3.23669715e-15, 1.27909314e02],
            ]
        )
        self.vector_angle = np.array([[-0.05140024, 0.08863197, -0.06735346]])
        self.cell_axis = np.array([[1], [0], [0]])
        self.cell_length = 10.2728

    def test_vector_scaling(self):
        """
        Defines the expected output for scaling the test vector
        """

        expected_scaled_vector = np.array([[-0.05140024, 0.08863197, -0.06735346]])

        output_scaled_vector = self.test.scale_vector(self.vector, self.G)

        self.assertEqual(expected_scaled_vector.all(), output_scaled_vector.all())

    def test_angle_calculation(self):
        """
        Defines the expected output for calculating the angle between two vectors
        """

        expected_angle = 119.9
        expected_supplementary_angle = 60.1

        output_angle, output_supplementary = self.test.calculate_angle(
            self.vector_angle, self.cell_axis, self.cell_length, self.G
        )

        self.assertEqual(round(output_angle, 1), expected_angle)
        self.assertEqual(round(output_supplementary, 1), expected_supplementary_angle)
