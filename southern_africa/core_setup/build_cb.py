import os
from helpers.windows_filesystem import get_dropbox_location
from gridded_sims.run.core_cb_setup import basic_gridded_config_builder, set_executable

from gridded_sims.run.site import *

# Run parameters
simulation_duration_days = 365*10

# Hardcoded:
def build_project_cb(num_cores=2):
    cb = basic_gridded_config_builder()

    project_folder = get_project_folder()
    set_executable(cb, os.path.join(project_folder, "bin/malaria_ongoing_build_185/"))
    # set_executable(cb, os.path.join(project_folder, "bin/original/"))

    cb.update_params({
        "Num_Cores": num_cores,
        "Simulation_Duration": simulation_duration_days,
        "Demographics_Filenames": ["Demographics/demo.json"],

        "Climate_Model": "CLIMATE_BY_DATA",
        "Air_Temperature_Filename": "Climate/Mozambique_30arcsec_air_temperature_daily.bin",
        "Land_Temperature_Filename": "Climate/Mozambique_30arcsec_air_temperature_daily.bin",
        "Rainfall_Filename": "Climate/Mozambique_30arcsec_rainfall_daily.bin",
        "Relative_Humidity_Filename": "Climate/Mozambique_30arcsec_relative_humidity_daily.bin",

        "Migration_Model": "FIXED_RATE_MIGRATION",
        "Migration_Pattern": "SINGLE_ROUND_TRIPS",

        "Enable_Local_Migration": 1,
        "Local_Migration_Roundtrip_Duration": 2,  # mean of exponential days-at-destination distribution
        "Local_Migration_Roundtrip_Probability": 1,  # fraction that return
        "Local_Migration_Filename": "Migration/local_migration.bin",

        'Enable_Regional_Migration': 1,
        'x_Regional_Migration': 0.03,
        'Regional_Migration_Roundtrip_Duration': 3,
        'Regional_Migration_Roundtrip_Probability': 1,
        'Regional_Migration_Filename': 'Migration/_Regional_Migration.bin',

        'Vector_Species_Names': ['arabiensis', 'funestus']
    })


    return cb


def get_project_folder():
    dropbox_folder = get_dropbox_location()
    project_folder = os.path.join(dropbox_folder, "projects/mz_magude/")
    return project_folder


def catchment_cb_params(cb, catch):
    # Set input folder to catchment-specific location
    project_folder = get_project_folder()
    cb.set_input_files_root(os.path.join(project_folder, "dtk_simulation_input/COMPS_experiments/{}/".format(catch)))

    # Set catchment-specific entomology
    catchment_ento(cb, catch)

    # Set x_local for proper amplitude of local migration
    x_local_df = pd.read_csv(os.path.join(project_folder, "migration/x_local_catch.csv"), index_col="catch")
    x_local = np.float(x_local_df.loc[catch])
    cb.update_params({"x_Local_Migration": x_local})
