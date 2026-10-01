#!/usr/bin/env python3

###################################################################################################
# ----------------------------------------CX-ASAP: lst_read---------------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Config
import logging
import re

# ----------Class Definition----------#


class LST_Read:
    """Parses data from SHELXL .lst output files.

    Provides methods to extract atom fractional coordinates from the embedded
    .res block and to extract plane equations and inter-plane angles from the
    Least-squares planes (MPLA) section.
    """

    # SHELXL instruction keywords — lines starting with these are not atom records
    _SHELXL_KEYWORDS = {
        "TITL",
        "CELL",
        "ZERR",
        "LATT",
        "SYMM",
        "SFAC",
        "UNIT",
        "L.S.",
        "MPLA",
        "PLAN",
        "TEMP",
        "CONF",
        "BOND",
        "LIST",
        "FMAP",
        "MORE",
        "SHEL",
        "WGHT",
        "FVAR",
        "AFIX",
        "PART",
        "RESI",
        "HKLF",
        "END",
        "SIMU",
        "DELU",
        "RIGU",
        "ISOR",
        "DFIX",
        "DANG",
        "SAME",
        "SADI",
        "CHIV",
        "FLAT",
        "DEFS",
        "BLOC",
        "OMIT",
        "TWIN",
        "BASF",
        "MERG",
        "ACTA",
        "SIZE",
        "HTAB",
        "RTAB",
        "EQIV",
        "REM",
        "EXTI",
        "SWAT",
    }

    # Matches an atom record: LABEL  SFAC_NUM  x  y  z  ...
    # SFAC_NUM is a single digit 1-9; x/y/z are decimal fractions
    _ATOM_RE = re.compile(
        r"^\s*([A-Za-z][A-Za-z0-9]{0,3})\s+([1-9])\s+"
        r"([-+]?\d*\.?\d+)\s+([-+]?\d*\.?\d+)\s+([-+]?\d*\.?\d+)"
    )

    def __init__(self, test_mode: bool = False) -> None:
        """Initialises the class.

        Args:
            test_mode (bool): if True, skips conf.yaml loading
        """

        self.test_mode = test_mode
        config = Config(self.test_mode)
        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

    def read(self, file_name: str) -> list:
        """Reads a .lst file and returns its lines.

        Args:
            file_name (str): full path to the .lst file

        Returns:
            lines (list): list of strings, one per line
        """

        with open(file_name, "rt") as f:
            return f.readlines()

    def extract_atom_coordinates(self, data: list) -> dict:
        """Extracts fractional atom coordinates from the embedded .res block of a .lst file.

        Parses atom records between the FVAR line and the HKLF line. Each atom
        record has the form:
            LABEL  SFAC_NUM  x  y  z  sof  [adp params ...]

        Continuation lines (starting with whitespace after '=') and AFIX/hydrogen
        lines are skipped automatically because they either start with whitespace
        or use SFAC_NUM 2 for H atoms (which are still captured if wanted).

        Args:
            data (list): lines of the .lst file as returned by read()

        Returns:
            coords (dict): atom label (uppercase) -> [x, y, z] as floats
        """

        coords = {}
        in_res = False

        for line in data:
            tokens = line.split()
            if not tokens:
                continue

            keyword = tokens[0].upper()

            if keyword.startswith("FVAR"):
                in_res = True
                continue

            if not in_res:
                continue

            if keyword == "HKLF":
                break

            if keyword in self._SHELXL_KEYWORDS:
                continue

            m = self._ATOM_RE.match(line)
            if m:
                label = m.group(1).upper()
                x = float(m.group(3))
                y = float(m.group(4))
                z = float(m.group(5))
                coords[label] = [x, y, z]

        return coords

    def extract_plane_normals(self, data: list) -> list:
        """Extracts the plane normal coefficients from the Least-squares planes section.

        Each MPLA plane is reported as:
            h x + k y + l z = d
        where h, k, l are the normal coefficients in crystal (fractional) coordinates.

        Args:
            data (list): lines of the .lst file as returned by read()

        Returns:
            planes (list): list of [h, k, l] float arrays, one per plane, in order
        """

        planes = []
        in_planes = False

        # Matches the plane equation line, e.g.:
        #   3.2122 (0.0266) x +  2.0735 (0.0414) y +  7.0389 (0.0168) z =  4.1556 (0.0111)
        plane_re = re.compile(r"([-+]?\s*\d*\.?\d+)\s*\([^)]*\)\s*x\s*[+-]")
        coeff_re = re.compile(r"([-+]?\s*\d*\.?\d+)\s*\([^)]*\)\s*[xyz]")

        for line in data:
            if "Least-squares planes" in line:
                in_planes = True
                continue

            if not in_planes:
                continue

            # Stop at the next major section
            if in_planes and line.strip() and line.strip()[0] == "R" and "R1" in line:
                break

            if plane_re.search(line):
                coeffs = coeff_re.findall(line)
                if len(coeffs) >= 3:
                    try:
                        h = float(coeffs[0].replace(" ", ""))
                        k = float(coeffs[1].replace(" ", ""))
                        l = float(coeffs[2].replace(" ", ""))
                        planes.append([h, k, l])
                    except ValueError:
                        logging.warning(
                            __name__
                            + f" : Could not parse plane equation: {line.strip()}"
                        )

        return planes

    def extract_interplane_angle(self, data: list) -> float:
        """Extracts the angle between consecutive MPLA planes from a .lst file.

        SHELXL outputs 'Angle to previous plane (with approximate esd) = VALUE ( ESD )'
        between each pair of consecutive MPLA plane blocks.

        Args:
            data (list): lines of the .lst file as returned by read()

        Returns:
            angle (float): the angle in degrees between the first two planes,
                           or 0.0 if not found
        """

        for line in data:
            if "Angle to previous plane" in line:
                m = re.search(r"=\s*([\d.]+)\s*\(", line)
                if m:
                    return float(m.group(1))

        logging.warning(__name__ + " : No inter-plane angle found in .lst file")
        return 0.0
