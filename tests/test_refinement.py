#!/usr/bin/env python

import logging
import statistics
import unittest

from data_refinement.modules.refinement import Structure_Refinement


class testRefinement(unittest.TestCase):
    def setUp(self):
        """
        Defines example cell + reference .ins files
        """
        self.example_new_ins = """TITL jkc18jc03bVT_07_240_00 in P2(1)/n
REM P2(1)/n (#14 in standard setting)
CELL 0.71073  10.277258   4.686109  11.333463  90.0000  92.0179  90.0000
ZERR    2.00   0.000809   0.000459   0.001020   0.0000   0.0084   0.0000
LATT  1
SYMM -x+1/2, y+1/2,-z+1/2
SFAC C H O Cu
UNIT 20.00 28.00 8.00 2.00
TREF
HKLF 4
END"""
        self.example_reference_ins = """TITL jkc18jc03bVT_11_200_00 in P2(1)/n
    200.res
    created by SHELXL-2018/3 at 10:50:40 on 29-Jun-2022
REM P2(1)/n (#14 in standard setting)
CELL 0.71073  10.272757   4.673770  11.309690  90.0000  92.1849  90.0000
ZERR    2.00   0.000643   0.000360   0.000816   0.0000   0.0068   0.0000
LATT  1
SYMM 1/2-X, 1/2+Y, 1/2-Z
SFAC C H O CU
UNIT 20 28 8 2
L.S. 10
BOND
LIST 4
CONF
ACTA
SHEL 10 0.73
FMAP 2
PLAN 20
WGHT    0.035500    0.272400
FVAR       8.16324
CU1   4    0.500000    0.000000    0.500000    10.50000    0.02009    0.01932 =
         0.02403   -0.00402    0.00608   -0.00255
C1    1    0.817490    0.558904    0.552847    11.00000    0.02154    0.03047 =
         0.03396   -0.00356    0.00322   -0.00483
AFIX 137
H1A   2    0.876973    0.422502    0.521548    11.00000   -1.50000
H1B   2    0.852765    0.632373    0.626581    11.00000   -1.50000
H1C   2    0.804671    0.713394    0.497688    11.00000   -1.50000
AFIX   0
O1    3    0.662054    0.203645    0.506334    11.00000    0.02197    0.02175 =
         0.02953   -0.00402    0.00694   -0.00237
O2    3    0.434749    0.197701    0.634362    11.00000    0.02457    0.02205 =
         0.02654   -0.00314    0.00825   -0.00334
C2    1    0.689320    0.415690    0.572747    11.00000    0.01997    0.01902 =
         0.02322    0.00298   -0.00100   -0.00015
C4    1    0.488366    0.410679    0.685935    11.00000    0.02563    0.01988 =
         0.01907    0.00076    0.00098    0.00229
C3    1    0.609127    0.525644    0.659431    11.00000    0.02586    0.02437 =
         0.02275   -0.00455    0.00112   -0.00468
AFIX  43
H3    2    0.638249    0.685133    0.701901    11.00000   -1.20000
AFIX   0
C5    1    0.413385    0.544845    0.783423    11.00000    0.03168    0.03378 =
         0.02897   -0.00916    0.00980   -0.00275
AFIX 137
H5A   2    0.374673    0.720507    0.755610    11.00000   -1.50000
H5B   2    0.471428    0.583505    0.850042    11.00000   -1.50000
H5C   2    0.346245    0.416068    0.806776    11.00000   -1.50000
AFIX   0
HKLF 4




REM  jkc18jc03bVT_11_200_00 in P2(1)/n
REM wR2 = 0.0858, GooF = S = 1.106, Restrained GooF = 1.106 for all data
REM R1 = 0.0324 for 1197 Fo > 4sig(Fo) and 0.0401 for all 1401 data
REM 72 parameters refined using 0 restraints

END"""

        self.example_combined_ins = """TITL jkc18jc03bVT_07_240_00 in P2(1)/n
REM P2(1)/n (#14 in standard setting)
CELL 0.71073  10.277258   4.686109  11.333463  90.0000  92.0179  90.0000
ZERR    2.00   0.000809   0.000459   0.001020   0.0000   0.0084   0.0000
LATT  1
SYMM 1/2-X, 1/2+Y, 1/2-Z
SFAC C H O CU
UNIT 20 28 8 2
L.S. 10
BOND
LIST 4
CONF
ACTA
SHEL 10 0.73
FMAP 2
PLAN 20
WGHT    0.035500    0.272400
FVAR       8.16324
CU1   4    0.500000    0.000000    0.500000    10.50000    0.02009    0.01932 =
         0.02403   -0.00402    0.00608   -0.00255
C1    1    0.817490    0.558904    0.552847    11.00000    0.02154    0.03047 =
         0.03396   -0.00356    0.00322   -0.00483
AFIX 137
H1A   2    0.876973    0.422502    0.521548    11.00000   -1.50000
H1B   2    0.852765    0.632373    0.626581    11.00000   -1.50000
H1C   2    0.804671    0.713394    0.497688    11.00000   -1.50000
AFIX   0
O1    3    0.662054    0.203645    0.506334    11.00000    0.02197    0.02175 =
         0.02953   -0.00402    0.00694   -0.00237
O2    3    0.434749    0.197701    0.634362    11.00000    0.02457    0.02205 =
         0.02654   -0.00314    0.00825   -0.00334
C2    1    0.689320    0.415690    0.572747    11.00000    0.01997    0.01902 =
         0.02322    0.00298   -0.00100   -0.00015
C4    1    0.488366    0.410679    0.685935    11.00000    0.02563    0.01988 =
         0.01907    0.00076    0.00098    0.00229
C3    1    0.609127    0.525644    0.659431    11.00000    0.02586    0.02437 =
         0.02275   -0.00455    0.00112   -0.00468
AFIX  43
H3    2    0.638249    0.685133    0.701901    11.00000   -1.20000
AFIX   0
C5    1    0.413385    0.544845    0.783423    11.00000    0.03168    0.03378 =
         0.02897   -0.00916    0.00980   -0.00275
AFIX 137
H5A   2    0.374673    0.720507    0.755610    11.00000   -1.50000
H5B   2    0.471428    0.583505    0.850042    11.00000   -1.50000
H5C   2    0.346245    0.416068    0.806776    11.00000   -1.50000
AFIX   0
HKLF 4




REM  jkc18jc03bVT_11_200_00 in P2(1)/n
REM wR2 = 0.0858, GooF = S = 1.106, Restrained GooF = 1.106 for all data
REM R1 = 0.0324 for 1197 Fo > 4sig(Fo) and 0.0401 for all 1401 data
REM 72 parameters refined using 0 restraints

END"""
        self.test = Structure_Refinement(test_mode=True)

    def test_merge_data(self):
        """
        Tests that the combination of cell + reference structure works
        """
        combined = self.test.merge_data(self.example_reference_ins, self.example_new_ins)
        self.assertEqual(combined, self.example_combined_ins)

    def test_convergence_check(self):
        """
        Test 1 - shifts refined, weights different - not converged
        """

        convergence = False
        refinement_failed = False
        shift = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        refinements = 8
        tolerance = 0.002
        weights_old = [0.2, 0]
        weights_new = [0.1, 0]

        expected_convergence = False
        expected_refinement_failed = False

        convergence_out, refinement_failed_out = self.test.converge(
            convergence,
            refinement_failed,
            shift,
            refinements,
            tolerance,
            weights_old,
            weights_new,
        )

        self.assertEqual(convergence_out, expected_convergence)
        self.assertEqual(refinement_failed_out, expected_refinement_failed)

        """
        Test 2 - shifts refined, weights equal - converged
        """

        convergence = False
        refinement_failed = False
        shift = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        refinements = 8
        tolerance = 0.002
        weights_old = [0.2, 0]
        weights_new = [0.2, 0]

        expected_convergence = True
        expected_refinement_failed = False

        convergence_out, refinement_failed_out = self.test.converge(
            convergence,
            refinement_failed,
            shift,
            refinements,
            tolerance,
            weights_old,
            weights_new,
        )

        self.assertEqual(convergence_out, expected_convergence)
        self.assertEqual(refinement_failed_out, expected_refinement_failed)

        """
        Test 3 - shifts not refined, weights equal - not converged
        """

        convergence = False
        refinement_failed = False
        shift = [1, 1, 1, 1, 1, 1, 1, 1]
        refinements = 8
        tolerance = 0.002
        weights_old = [0.2, 0]
        weights_new = [0.2, 0]

        expected_convergence = False
        expected_refinement_failed = False

        convergence_out, refinement_failed_out = self.test.converge(
            convergence,
            refinement_failed,
            shift,
            refinements,
            tolerance,
            weights_old,
            weights_new,
        )

        self.assertEqual(convergence_out, expected_convergence)
        self.assertEqual(refinement_failed_out, expected_refinement_failed)

        """
        Test 4 - shifts not refined, weights not equal - not converged
        """

        convergence = False
        refinement_failed = False
        shift = [1, 1, 1, 1, 1, 1, 1, 1]
        refinements = 8
        tolerance = 0.002
        weights_old = [0.2, 0]
        weights_new = [0.1, 0]

        expected_convergence = False
        expected_refinement_failed = False

        convergence_out, refinement_failed_out = self.test.converge(
            convergence,
            refinement_failed,
            shift,
            refinements,
            tolerance,
            weights_old,
            weights_new,
        )

        self.assertEqual(convergence_out, expected_convergence)
        self.assertEqual(refinement_failed_out, expected_refinement_failed)
