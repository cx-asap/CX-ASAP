#!/usr/bin/env python3

###################################################################################################
# -------------------------------------CX-ASAP: crystal_math--------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

import numpy as np

# ----------Functions----------#


def metric_matrix(cell_params: list) -> np.ndarray:
    """Calculates the metric matrix G for a unit cell.

    G encodes the inner products of the real-space basis vectors.

    Args:
        cell_params (list): [a, b, c, alpha, beta, gamma] where a/b/c are in
                            Angstroms and angles are in degrees

    Returns:
        G (np.ndarray): 3x3 real symmetric metric matrix
    """

    a, b, c = cell_params[0], cell_params[1], cell_params[2]
    alpha = np.radians(cell_params[3])
    beta = np.radians(cell_params[4])
    gamma = np.radians(cell_params[5])

    G = np.array(
        [
            [a**2, a * b * np.cos(gamma), a * c * np.cos(beta)],
            [a * b * np.cos(gamma), b**2, b * c * np.cos(alpha)],
            [a * c * np.cos(beta), b * c * np.cos(alpha), c**2],
        ]
    )

    return G


def orthonorm_matrix(cell_params: list) -> np.ndarray:
    """Builds the orthonormalisation matrix M for a unit cell.

    Converts fractional row vectors to Cartesian coordinates:
        cart_row = frac_row @ M

    The convention used is the standard crystallographic one where:
        - a is along x
        - b is in the xy plane
        - c is oriented to complete a right-handed system

    Args:
        cell_params (list): [a, b, c, alpha, beta, gamma] where a/b/c are in
                            Angstroms and angles are in degrees

    Returns:
        M (np.ndarray): 3x3 orthonormalisation matrix
    """

    a, b, c = cell_params[0], cell_params[1], cell_params[2]
    alpha = np.radians(cell_params[3])
    beta = np.radians(cell_params[4])
    gamma = np.radians(cell_params[5])

    alpha_star = np.arccos(
        (np.cos(beta) * np.cos(gamma) - np.cos(alpha)) / (np.sin(beta) * np.sin(gamma))
    )

    M = np.array(
        [
            [a, b * np.cos(gamma), c * np.cos(beta)],
            [0, b * np.sin(gamma), -c * np.sin(beta) * np.cos(alpha_star)],
            [0, 0, c * np.sin(beta) * np.sin(alpha_star)],
        ]
    )

    return M
