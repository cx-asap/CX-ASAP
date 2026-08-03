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

    # Accept a list (e.g. from YAML block syntax) or a plain string
    if isinstance(symm_str, list):
        symm_str = symm_str[0]

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


def fractional_to_cartesian(
    point_frac: "np.ndarray | list", cell_params: list
) -> np.ndarray:
    """Converts one fractional point into Cartesian coordinates."""

    M = orthonorm_matrix(cell_params)
    return np.dot(np.array(point_frac, dtype=float), M.T)


def distance_between_points(point_1: "np.ndarray | list", point_2: "np.ndarray | list") -> float:
    """Calculates Euclidean distance between two Cartesian points."""

    p1 = np.array(point_1, dtype=float)
    p2 = np.array(point_2, dtype=float)
    return float(np.linalg.norm(p1 - p2))


def angle_between_points(
    point_1: "np.ndarray | list",
    point_2: "np.ndarray | list",
    point_3: "np.ndarray | list",
) -> float:
    """Calculates the angle in degrees for points 1-2-3 at point 2."""

    p1 = np.array(point_1, dtype=float)
    p2 = np.array(point_2, dtype=float)
    p3 = np.array(point_3, dtype=float)

    v1 = p1 - p2
    v2 = p3 - p2

    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)

    if norm1 == 0.0 or norm2 == 0.0:
        raise ValueError("Zero-length vector in angle calculation")

    cos_theta = np.dot(v1, v2) / (norm1 * norm2)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_theta)))


def torsion_between_points(
    point_1: "np.ndarray | list",
    point_2: "np.ndarray | list",
    point_3: "np.ndarray | list",
    point_4: "np.ndarray | list",
) -> float:
    """Calculates torsion angle in degrees for points 1-2-3-4."""

    c1 = np.array(point_1, dtype=float)
    c2 = np.array(point_2, dtype=float)
    c3 = np.array(point_3, dtype=float)
    c4 = np.array(point_4, dtype=float)

    b0 = c2 - c1
    b1 = c3 - c2
    b2 = c4 - c3

    b1_norm = np.linalg.norm(b1)
    if b1_norm == 0.0:
        raise ValueError("Zero-length central bond in torsion calculation")

    b1_unit = b1 / b1_norm
    v = b0 - np.dot(b0, b1_unit) * b1_unit
    w = b2 - np.dot(b2, b1_unit) * b1_unit

    v_norm = np.linalg.norm(v)
    w_norm = np.linalg.norm(w)
    if v_norm == 0.0 or w_norm == 0.0:
        raise ValueError("Degenerate geometry in torsion calculation")

    x = np.dot(v, w)
    y = np.dot(np.cross(b1_unit, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def point_to_plane_distance(
    point: "np.ndarray | list", plane_points: "list[np.ndarray] | np.ndarray | list"
) -> float:
    """Calculates absolute Cartesian distance from a point to a best-fit plane."""

    p = np.array(point, dtype=float)
    plane = np.array(plane_points, dtype=float)

    if len(plane) < 3:
        raise ValueError("Need at least 3 points to define a plane")

    centroid = np.mean(plane, axis=0)
    centered = plane - centroid
    _, _, vh = np.linalg.svd(centered)
    normal = vh[-1]

    normal_norm = np.linalg.norm(normal)
    if normal_norm == 0.0:
        raise ValueError("Could not resolve plane normal")

    normal_unit = normal / normal_norm
    return float(abs(np.dot(p - centroid, normal_unit)))


def reciprocal_orthonorm_matrix(cell_params: list) -> np.ndarray:
    """Returns the reciprocal-space transform matrix M* = inv(M).

    This is useful for converting vectors between Cartesian and fractional
    representations in reciprocal space.
    """

    M = orthonorm_matrix(cell_params)
    return np.linalg.inv(M)


def cartesian_plane_normal_to_fractional(
    normal_cart: "np.ndarray | list", cell_params: list
) -> np.ndarray:
    """Converts a Cartesian plane normal into fractional representation.

    Args:
        normal_cart: plane normal components in Cartesian basis
        cell_params: [a, b, c, alpha, beta, gamma]

    Returns:
        np.ndarray: length-3 fractional vector
    """

    n_cart = np.array(normal_cart, dtype=float)
    M_star = reciprocal_orthonorm_matrix(cell_params)
    return np.dot(n_cart, M_star)


def angle_between_vectors(
    vector_1: "np.ndarray | list",
    vector_2: "np.ndarray | list",
    fold_to_acute: bool = False,
) -> float:
    """Calculates the angle in degrees between two vectors.

    Args:
        vector_1: first vector
        vector_2: second vector
        fold_to_acute: if True, folds obtuse results into [0, 90]

    Returns:
        float: angle in degrees
    """

    v1 = np.array(vector_1, dtype=float)
    v2 = np.array(vector_2, dtype=float)

    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 == 0.0 or n2 == 0.0:
        raise ValueError("Zero-length vector in angle calculation")

    cos_theta = np.dot(v1, v2) / (n1 * n2)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    angle = float(np.degrees(np.arccos(cos_theta)))

    if fold_to_acute and angle > 90.0:
        angle = 180.0 - angle

    return angle


def best_fit_plane_normal(
    plane_points: "list[np.ndarray] | np.ndarray | list",
) -> np.ndarray:
    """Returns a unit normal vector for a best-fit Cartesian plane.

    Args:
        plane_points: array-like collection of 3D Cartesian points

    Returns:
        np.ndarray: length-3 unit normal vector
    """

    plane = np.array(plane_points, dtype=float)
    if len(plane) < 3:
        raise ValueError("Need at least 3 points to define a plane")

    centroid = np.mean(plane, axis=0)
    centered = plane - centroid
    _, _, vh = np.linalg.svd(centered)
    normal = vh[-1]

    normal_norm = np.linalg.norm(normal)
    if normal_norm == 0.0:
        raise ValueError("Could not resolve plane normal")

    return normal / normal_norm
