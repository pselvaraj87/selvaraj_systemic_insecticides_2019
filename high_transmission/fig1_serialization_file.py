from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.utils.builders.sweep import GenericSweepBuilder
from dtk.vector.species import set_larval_habitat, set_species_param
from dtk.generic.serialization import add_SerializationTimesteps
from dtk.generic.climate import set_climate_constant
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser


if __name__ == "__main__":

    expname = 'Serialized_file_local_IVM'

    sim_duration = 365 * 5
    num_seeds = 2

    SetupParser('LOCAL')

    cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')

    builder = GenericSweepBuilder.from_dict({'Run_Number': range(num_seeds)})

    set_climate_constant(cb)
    set_species_param(cb, 'gambiae', 'Larval_Habitat_Types',
                      {"LINEAR_SPLINE": {
                          "Capacity_Distribution_Number_Of_Years": 1,
                          "Capacity_Distribution_Over_Time": {
                              "Times": [0.0, 30.417, 60.833, 91.25, 121.667, 152.083,
                                        182.5, 212.917, 243.333, 273.75, 304.167, 334.583],
                              "Values": [3, 0.8, 1.25, 0.1, 2.7, 10, 6, 35, 2.8, 1.5, 1.6, 2.1] # with 'x_Temporary_Larval_Habitat': 0.2 for EIR~10
                          },
                          "Max_Larval_Capacity": pow(10, 9)
                      }})
    set_species_param(cb, "gambiae", "Adult_Life_Expectancy", 20)
    set_species_param(cb, "gambiae", "Indoor_Feeding_Fraction", 0.9)

    cb.update_params({

                      'Age_Initialization_Distribution_Type': 'DISTRIBUTION_SIMPLE',
                      'Antigen_Switch_Rate': pow(10, -9.116590124),
                      'Base_Gametocyte_Production_Rate': 0.06150582,
                      'Base_Gametocyte_Mosquito_Survival_Rate': 0.002011099,
                      'Base_Population_Scale_Factor': 1,            # Change x_Temporary_Larval_Habitat by same factor if changing Base_Population_Scale_Factor
                      'x_Temporary_Larval_Habitat': 0.2,

                      'Birth_Rate_Dependence': 'FIXED_BIRTH_RATE',
                      "Death_Rate_Dependence": "NONDISEASE_MORTALITY_BY_AGE_AND_GENDER",

                      "Demographics_Filenames": ['Interventions/single_node_demographics.json',
                                   'Interventions/single_node_drugs_overlay.json'],

                      "Disable_IP_Whitelist": 1,

                      'Enable_Birth': 1,
                      'Enable_Default_Reporting': 1,
                      'logLevel_default': 'ERROR',
                      'Enable_Vital_Dynamics': 1,
                      'Enable_Property_Output': 1,

                      'Falciparum_MSP_Variants': 32,
                      'Falciparum_Nonspecific_Types': 76,
                      'Falciparum_PfEMP1_Variants': 1070,
                      'Gametocyte_Stage_Survival_Rate': 0.588569307,

                      'MSP1_Merozoite_Kill_Fraction': 0.511735322,
                      'Max_Individual_Infections': 3,
                      'Nonspecific_Antigenicity_Factor': 0.415111634,

                      "Simulation_Duration": sim_duration + 1,
                      "Run_Number": 0,
                      "Vector_Species_Names": ['gambiae']
                      })

    add_SerializationTimesteps(cb, [sim_duration], end_at_final=False)

    run_sim_args = {'config_builder': cb,
                    'exp_name': expname,
                    'exp_builder': builder}

    exp_manager = ExperimentManagerFactory.from_setup()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)