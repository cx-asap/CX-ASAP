#!/usr/bin/env python3

###################################################################################################
# ---------------------------------CX-ASAP: instrument cif generation------------------------------#
# ---Authors: Amy J. Thompson, Kate M. Smith, Daniel J. Eriksson, Jack K. Clegg & Jason R. Price---#
# -----------------------------------Python Implementation by AJT----------------------------------#
# -----------------------------------Project Design by JRP and JKC---------------------------------#
# --------------------------------Valuable Coding Support by KMS & DJE-----------------------------#
###################################################################################################

# ----------Required Modules----------#

from system_files.utils import Nice_YAML_Dumper, Config
from CifFile import ReadCif, CifBlock, CifFile
import pathlib
import os
import logging

# ----------Class Definition----------#


class Instrument_CIF:
    def __init__(self, test_mode: bool = False) -> None:

        """Initialises the class

        Sets up the yaml parameters input by the user

        Also defines the location for the system yaml file

        which stores a yaml of code-only parameters accessible throughout

        the software package

        Args:
            test_mode (bool): Automatically false, if true it will

                            make the functions compatible with the testing script
        """

        # Setup yaml files and logger

        config = Config(test_mode)

        self.cfg = config.cfg
        self.sys = config.sys
        self.conf_path = config.conf_path
        self.sys_path = config.sys_path

        self.instrument_parameters = [
            "_chemical_formula_moiety",
            "_chemical_absolute_configuration",
            "_exptl_absorpt_correction_T_max",
            "_exptl_absorpt_correction_T_min",
            "_exptl_absorpt_correction_type",
            "_exptl_crystal_colour",
            "_exptl_crystal_description",
            "_exptl_crystal_size_max",
            "_exptl_crystal_size_mid",
            "_exptl_crystal_size_min",
            "_diffrn_detector",
            "_diffrn_detector_area_resol_mean",
            "_diffrn_detector_type",
            "_diffrn_measurement_device",
            "_diffrn_measurement_device_type",
            "_diffrn_measurement_method",
            "_diffrn_radiation_monochromator",
            "_diffrn_radiation_probe",
            "_diffrn_radiation_type",
            "_diffrn_radiation_wavelength",
            "_diffrn_source",
            "_diffrn_source_type",
            "_computing_cell_refinement",
            "_computing_data_collection",
            "_computing_data_reduction",
            "_computing_molecular_graphics",
            "_computing_publication_material",
            "_computing_structure_refinement",
            "_computing_structure_solution",
        ]

    def read_reference_cif(self, reference: str) -> None:

        """Imports a fully completed CIF file which is used as a reference

        This is done using the PyCIFRW library

        Args:
            reference (str): Full path to the instrument cif file
        """

        reference = pathlib.Path(reference)

        os.chdir(reference.parent)

        try:
            self.reference_cif = ReadCif(reference.name)
        except:
            logging.critical(__name__ + " : Could not read reference Cif!")
            print("Error! Check logs")
            exit()

    def make_instrument_cif(self) -> None:

        """Extracts the instrument parameters defined in the _init_ function

        from the reference CIF.

        It saves these into a separated file called "instrument.cif"
        """

        data_block = self.reference_cif.first_block()
        data_contained = data_block.GetItemOrder()

        instrument_cif = CifFile()
        instrument_data_block = CifBlock()
        instrument_cif["instrument_information"] = instrument_data_block

        for item in self.instrument_parameters:
            try:
                instrument_cif["instrument_information"][item] = data_block[item]
            except KeyError:
                logging.info(
                    __name__ + " : parameter " + item + " not in reference cif"
                )

        with open("instrument.cif", "w") as f:
            f.write(instrument_cif.WriteOut())
