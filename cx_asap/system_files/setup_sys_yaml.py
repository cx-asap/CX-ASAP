# Generate sys.yaml - import statement down here just in case the user does not already have it and it gets installed for the first time in the above setup code

import os
import yaml
import pathlib


class Nice_YAML_Dumper(yaml.SafeDumper):
    def write_line_break(self, data: dict = None) -> None:
        """Makes the yaml have better formatting when edited

        Unsure of the specifics, got this code from StackOverflow
        """

        super().write_line_break(data)

        if len(self.indents) == 1:
            super().write_line_break()

sys_path = pathlib.Path.cwd() / "cx_asap" / "system_files" / "sys.yaml"

sys_params = {
    "module-refinement": {
        "structure_location": "test_data/data/200K/200.ins",
        "reference_path": "test_data/ref/ref.res",
        "refinements_to_check": 4,
        "tolerance": 0.005,
        "maximum_cycles": 10,
    },
    "pipeline-refinement": {
        "experiment_location": "test_data/data",
        "reference_path": "test_data/ref/ref.res",
        "refinements_to_check": 4,
        "tolerance": 0.005,
        "maximum_cycles": 10,
    },
    "module-cif-merge": {
        "instrument_cif": "test_data/ref/instrument.cif",
        "new_cif": "test_data/data/200K/200.cif",
    },
    "module-make-instrument-cif": {"reference_cif": "test_data/ref/ref.cif"},
    "module-cif-read": {
        "cif_parameters": [
            "_cell_length_a",
            "_cell_length_b",
            "_cell_length_c",
            "_cell_angle_alpha",
            "_cell_angle_beta",
            "_cell_angle_gamma",
            "_cell_volume",
            "_diffrn_reflns_av_R_equivalents",
            "_diffrn_measured_fraction_theta_full",
            "_diffrn_ambient_temperature",
            "_refine_ls_R_factor_gt",
        ],
        "structural_analysis_angles": True,
        "structural_analysis_bonds": True,
        "structural_analysis_torsions": True,
        "structural_analysis_hbonds": False,
        "folder_containing_cifs": "test_data/data",
        "ADP_analysis": True,
    },
    "pipeline-cif": {
        "experiment_location": "test_data/data",
        "instrument_ending": False,
        "instrument_file": "instrument.cif",
        "structure_solution": "shelxt",
        "chemical_formula": "Cu C10 H14 O4",
        "crystal_habit": "needle",
        "crystal_colour": "blue",
        "max_crystal_dimension": 0.2,
        "middle_crystal_dimension": 0.08,
        "min_crystal_dimension": 0.05,
    },
    "module-cell-analysis": {
        "csv_location": "test_data/data/CIF_Parameters.csv",
        "reference_unit_cell": "test_data/ref/ref.res",
        "x_axis_header": "_diffrn_ambient_temperature",
    },
    "module-structural-analysis": {
        "atoms_for_analysis": ["Cu1"],
        "bond_data": "test_data/data/Bond_Lengths.csv",
        "angle_data": "test_data/data/Bond_Angles.csv",
        "torsion_data": "test_data/data/Bond_Torsions.csv",
        "hbond_data": False,
    },
    "pipeline-variable-analysis": {
        "cif_parameters": [
            "_cell_length_a",
            "_cell_length_b",
            "_cell_length_c",
            "_cell_angle_alpha",
            "_cell_angle_beta",
            "_cell_angle_gamma",
            "_cell_volume",
            "_diffrn_reflns_av_R_equivalents",
            "_diffrn_measured_fraction_theta_full",
            "_diffrn_ambient_temperature",
            "_refine_ls_R_factor_gt",
        ],
        "atoms_for_analysis": ["Cu1"],
        "structural_analysis_angles": True,
        "structural_analysis_bonds": True,
        "structural_analysis_torsions": True,
        "structural_analysis_hbonds": False,
        "ADP_analysis": True,
        "varying_cif_parameter": "_diffrn_ambient_temperature",
        "experiment_location": "test_data/data",
        "reference_unit_cell": "test_data/ref/ref.res",
    },
    "pipeline-general": {
        "experiment_location": "test_data/data",
        "reference_cif_location": "test_data/ref/ref.cif",
        "refinements_to_check": 4,
        "tolerance": 0.005,
        "maximum_cycles": 10,
        "cif_parameters": [
            "_cell_length_a",
            "_cell_length_b",
            "_cell_length_c",
            "_cell_angle_alpha",
            "_cell_angle_beta",
            "_cell_angle_gamma",
            "_cell_volume",
            "_diffrn_reflns_av_R_equivalents",
            "_diffrn_measured_fraction_theta_full",
            "_diffrn_ambient_temperature",
            "_refine_ls_R_factor_gt",
        ],
        "structural_analysis_bonds": True,
        "structural_analysis_angles": True,
        "structural_analysis_torsions": True,
        "structural_analysis_hbonds": False,
        "ADP_analysis": True,
        "atoms_for_analysis": ["Cu1"],
        "varying_cif_parameter": "_diffrn_ambient_temperature",
        "varying_parameter_values": ["200(2)", "210(2)", "220(2)", "230(2)", "240(2)"],
    },
    "home_path": None,
    "analysis_path": None,
    "ref_path": None,
    "results_path": None,
    "failed_path": None,
    "ref_path_organised": None,
    "current_results_path": None,
    "process_counter": None,
    "hello there": "General Kenobi",
    "Structures_in_each_CIF": [None],
    "Successful_Positions": [None],
    "ref_a": None,
    "ref_b": None,
    "ref_c": None,
    "ref_alpha": None,
    "ref_beta": None,
    "ref_gamma": None,
    "ref_volume": None,
    "space_group": None,
    "frames_path": None,
    "XDS_inp_organised": None,
    "start_angle": None,
    "total_angle": None,
    "total_frames": None,
}

with open(sys_path, "w") as f:
    yaml.dump(
        sys_params,
        f,
        default_flow_style=False,
        Dumper=Nice_YAML_Dumper,
        sort_keys=False,
    )


# Resets error logs

error_logs = pathlib.Path.cwd() / "cx_asap" / "error_logs" / "error_output.txt"

if os.path.exists(error_logs):
    os.remove(error_logs)
