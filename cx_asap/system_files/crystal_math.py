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


def parse_symm_op(symm_str: str):
    """Parses a SHELXL/CIF symmetry operation string into a rotation matrix and
    translation vector (both in fractional coordinates).

    Converts strings like "-x+1/2, y+1/2, -z+1/2" into a 3x3 matrix R and a
    length-3 vector t such that:

        frac_new = frac @ R.T + t   (row-vector convention)

    Args:
        symm_str (str): symmetry operation string, e.g. "-x+1/2, y, -z+1/2"

    Returns:
        R (np.ndarray): 3x3 rotation/inversion matrix (entries 0, +1, or -1)
        t (np.ndarray): length-3 translation vector in fractional coordinates
    """

    import re

    R = np.zeros((3, 3))
    t = np.zeros(3)

    axes = {"x": 0, "y": 1, "z": 2}
    frac_re = re.compile(r"([+-]?\s*\d+)\s*/\s*(\d+)")
    coeff_re = re.compile(r"([+-]?\s*\d*\.?\d*)\s*([xyz])")

    for row, expr in enumerate(symm_str.split(",")):
        expr = expr.strip()

        # Extract and remove all fraction terms (e.g. +1/2, -1/3)
        for m in frac_re.finditer(expr):
            t[row] += float(m.group(1).replace(" ", "")) / float(m.group(2))
        expr_no_frac = frac_re.sub("", expr)

        # Extract remaining integer/float terms not attached to x/y/z
        # (bare numbers like +1 or -2 without axis letter)
        remainder = coeff_re.sub("", expr_no_frac).strip()
        remainder = remainder.replace(" ", "")
        if remainder and remainder not in ("+", "-", ""):
            try:
                t[row] += float(remainder)
            except ValueError:
                pass

        # Extract axis coefficients
        for m in coeff_re.finditer(expr_no_frac):
            coeff_str = m.group(1).replace(" ", "")
            axis = m.group(2)
            if coeff_str in ("", "+"):
                coeff = 1.0
            elif coeff_str == "-":
                coeff = -1.0
            else:
                coeff = float(coeff_str)
            R[row, axes[axis]] = coeff

    return R, t
