#!/usr/bin/env python

import unittest
from post_refinement_analysis.pipelines.variable_cif_parameter import (
    Variable_Analysis_Pipeline,
)
import pandas as pd


class testBehaviour(unittest.TestCase):
    def setUp(self):
        """
        Sets up a sample data frame and the expected behaviour of the temperatures in that dataframe
        """

        self.test = Variable_Analysis_Pipeline(test_mode=True)

        data = {
            "Temperature": [
                100,
                110,
                120,
                120,
                120,
                130,
                140,
                130,
                120,
                110,
                100,
                90,
                100,
            ]
        }

        self.sample_dataframe = pd.DataFrame(data)

        self.expected_behaviour = [
            "Increasing",
            "Increasing",
            "Increasing",
            "Did Not Change",
            "Did Not Change",
            "Increasing",
            "Maxima",
            "Decreasing",
            "Decreasing",
            "Decreasing",
            "Decreasing",
            "Minima",
            "Increasing",
        ]

    def test_behaviour(self):
        """
        Checks that the function 'determine_behaviour' gives the same output as expected

        This checks for all the types of behaviour present - increasing, decreasing, minima, maxima, and did not change
        """

        output_behaviour = self.test.determine_behaviour(
            self.sample_dataframe, "Temperature"
        )

        self.assertEqual(output_behaviour, self.expected_behaviour)
