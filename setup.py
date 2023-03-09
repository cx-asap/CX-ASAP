#!/usr/bin/env python3

#####################################################################################################
# -----------------------------------------CX-ASAP: setup------------------------------------------ #
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jason R. Price & Jack K. Clegg--- #
# -----------------------------------Python Implementation by AJT---------------------------------- #
# -----------------------------------Project Design by JRP and JKC--------------------------------- #
# --------------------------------Valuable Coding Support by KMS & DJE----------------------------- #
#####################################################################################################

from setuptools import setup, find_packages
import shutil
import pathlib
import os

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="CX-ASAP",
    version="1.0.2",
    description="A Library of Tools for Automated Crystallographic Analysis with the Australian Synchrotron",
    author="Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg, Jason R. Price",
    author_email="amy.thompson2@uqconnect.edu.au",
    package_dir={"": "cx_asap"},
    packages=find_packages(where="cx_asap"),
    install_requires=requirements,
    package_data={
        "": ["*.yaml", "*.ipynb", "*.md", "*.jpg", "*.ins", "*.cif", "*.hkl"],
    },
    entry_points={"console_scripts": ["cxasap=cxasap:run"]},
)

