from malaria.study_sites.site_setup_functions import summary_report_fn
from dtk.interventions.property_change import change_individual_property
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import set_species_param
from dtk.generic.climate import set_climate_constant
from Endectocides_paper.createSimDirectoryMap import createSimDirectoryMap
from malaria.interventions.malaria_drug_campaigns import add_drug_campaign
from dtk.interventions.ivermectin import add_ivermectin
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.SetupParser import SetupParser

import numpy as np


def add_ivermectin_group(cb, coverage=1.0, agemax=10, start_days=[60, 60+365], drug_code=30, smc_code='DP5', target_group='Everyone'):

    if target_group == 'Everyone':
        add_ivermectin(cb, drug_code=drug_code, coverage=1.0, start_days=start_days,
                       trigger_condition_list=['Received_Campaign_Drugs'],
                       target_group={'agemin': 0, 'agemax': 5}
                       )
        add_ivermectin(cb, drug_code=drug_code, coverage=coverage, start_days=start_days,
                       target_group={'agemin': 5, 'agemax': 200}
                       )

    elif target_group == 'Excluding<5':
        # SMC + Ivermectin for kids aged 6 to 15 and adults > 15 w/o women of child bearing age
        add_ivermectin(cb, drug_code=drug_code, coverage=coverage, start_days=start_days,
                       target_group={'agemin': 5, 'agemax': 200}
                       )

    elif target_group == 'Excludingwomenand<5':
        add_ivermectin(cb, drug_code=drug_code, coverage=coverage, start_days=start_days,
                       target_group={'agemin': 5, 'agemax': 200, 'gender': 'Male'})
        add_ivermectin(cb, drug_code=drug_code, coverage=coverage, start_days=start_days,
                       target_group={'agemin': 51, 'agemax': 200, 'gender': 'Female'})
        add_ivermectin(cb, drug_code=drug_code, coverage=coverage, start_days=start_days,
                       target_group={'agemin': 5, 'agemax': 12, 'gender': 'Female'}
                       )

    elif target_group == 'Excludingwomen':
        add_ivermectin(cb, drug_code=drug_code, coverage=1.0, start_days=start_days,
                       trigger_condition_list=['Received_Campaign_Drugs'],
                       target_group={'agemin': 0, 'agemax': 5}
                       )
        add_ivermectin(cb, drug_code=drug_code, coverage=coverage, start_days=start_days,
                       target_group={'agemin': 5, 'agemax': 200, 'gender': 'Male'})
        add_ivermectin(cb, drug_code=drug_code, coverage=coverage, start_days=start_days,
                       target_group={'agemin': 51, 'agemax': 200, 'gender': 'Female'})
        add_ivermectin(cb, drug_code=drug_code, coverage=coverage, start_days=start_days,
                       target_group={'agemin': 5, 'agemax': 12, 'gender': 'Female'}
                       )

    elif target_group == 'ExpandedSMC':
        # SMC + Ivermectin for adults > 15
        add_ivermectin(cb, drug_code=drug_code, coverage=1.0, start_days=start_days,
                       trigger_condition_list=['Received_Campaign_Drugs'],
                       target_group={'agemin': 5, 'agemax': agemax}
                       )
        add_ivermectin(cb, drug_code=drug_code, coverage=coverage, start_days=start_days,
                       target_group={'agemin': agemax, 'agemax': 200}
                       )

    elif target_group == 'ExpandedSMCexcludingwomen':
        if agemax < 12:
            add_ivermectin(cb, drug_code=drug_code, coverage=1.0, start_days=start_days,
                           trigger_condition_list=['Received_Campaign_Drugs'],
                           target_group={'agemin': 5, 'agemax': agemax, 'gender': 'Female'},
                           )
            add_ivermectin(cb, drug_code=drug_code, coverage=1.0, start_days=start_days,
                           target_group={'agemin': agemax, 'agemax': 12, 'gender': 'Female'},
                           )

        else:
            add_ivermectin(cb, drug_code=drug_code, coverage=1.0, start_days=start_days,
                           trigger_condition_list=['Received_Campaign_Drugs'],
                           target_group={'agemin': 5, 'agemax': agemax, 'gender': 'Female'},
                           )
        add_ivermectin(cb, drug_code=drug_code, coverage=1.0, start_days=start_days,
                       trigger_condition_list=['Received_Campaign_Drugs'],
                       target_group={'agemin': 5, 'agemax': agemax, 'gender': 'Male'},
                       )
        add_ivermectin(cb, drug_code=drug_code, coverage=coverage, start_days=start_days,
                       target_group={'agemin': agemax, 'agemax': 200, 'gender': 'Male'})
        add_ivermectin(cb, drug_code=drug_code, coverage=coverage, start_days=start_days,
                       target_group={'agemin': 51, 'agemax': 200, 'gender': 'Female'})

    elif target_group == 'ExpandedSMCexcludingwomenonly':
        if agemax < 12:
            add_ivermectin(cb, drug_code=drug_code, coverage=1.0, start_days=start_days,
                           trigger_condition_list=['Received_Campaign_Drugs'],
                           target_group={'agemin': 0, 'agemax': agemax, 'gender': 'Female'},
                           )
            add_ivermectin(cb, drug_code=drug_code, coverage=1.0, start_days=start_days,
                           target_group={'agemin': agemax, 'agemax': 12, 'gender': 'Female'},
                           )

        else:
            add_ivermectin(cb, drug_code=drug_code, coverage=1.0, start_days=start_days,
                           trigger_condition_list=['Received_Campaign_Drugs'],
                           target_group={'agemin': 0, 'agemax': agemax, 'gender': 'Female'},
                           )
        add_ivermectin(cb, drug_code=drug_code, coverage=1.0, start_days=start_days,
                       trigger_condition_list=['Received_Campaign_Drugs'],
                       target_group={'agemin': 0, 'agemax': agemax, 'gender': 'Male'},
                       )
        add_ivermectin(cb, drug_code=drug_code, coverage=coverage, start_days=start_days,
                       target_group={'agemin': agemax, 'agemax': 200, 'gender': 'Male'})
        add_ivermectin(cb, drug_code=drug_code, coverage=coverage, start_days=start_days,
                       target_group={'agemin': 51, 'agemax': 200, 'gender': 'Female'})

    return {'Intervention_type': '%s+%s+IV_%i' % (smc_code, target_group, drug_code)}


def add_smc_group(cb, coverage=1.0, start_days=[60], agemax=10, drug='DP',
                  target_property_name='ReceivedIntervention', target_property_value='NoDrugs'):

    add_drug_campaign(cb, 'SMC', drug, start_days=start_days, repetitions=4, tsteps_btwn_repetitions=30,
                      coverage=coverage,
                      target_group={'agemin': 0, 'agemax': agemax})
                      # ind_property_restrictions=[{target_property_name: target_property_value}])

    return {'Coverage': coverage, 'Start': start_days[0], 'Intervention_type': 'SMC%i' %agemax}


def add_pop_intervention(cb, target_property_name='ReceivedIntervention', target_property_value='NoDrugs', coverage=1.0):
    change_individual_property(cb, target_property_name=target_property_name, target_property_value=target_property_value,
                               coverage=coverage
                               )

    return {'Coverage': coverage}


def add_summary_report(cb, start_day=0):

    summary_report_fn(start=start_day+1, interval=1.0, description='Daily_Report',
                      age_bins=[5.0, 15.0, 100.0])(cb)

    return {'Summary_report_start_day': start_day}


def make_simmap(expname, filetype='Inset'):

    simmap = createSimDirectoryMap(expname)
    # simmap = get_filepath(simmap, filetype)

    return simmap


# Get filelist of files in simmap
def get_filepath(simmap, filetype='state-'):

    filelist = []
    for path in simmap['outpath']:
        outputpath = os.path.join(path, 'output')
        for file in os.listdir(outputpath):
            if filetype in file:
                filelist.append(os.path.join(outputpath, file))

    simmap['filelist'] = filelist

    return simmap


def get_outpath_for_serialized_file(simmap, seed):

    temp = simmap[simmap['Run_Number'] == seed]

    return temp[temp['Run_Number'] == seed]['outpath'].values[0]


if __name__ == "__main__":

    expname = 'exe_test'
    serialized_file_output = ['', '']

    sim_duration = 365 * 1
    num_seeds = 2

    SetupParser.init('LOCAL')

    cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')

    intervention_days = [180]
    coverages = [x for x in np.arange(0.0, 1.01, 0.1)]
    target_groups_IVM = ['Excludingwomen']
    expanded_SMC_groups_IVM = ['ExpandedSMCexcludingwomenonly']

    # DP + Ivermectin (expanded and regular)
    SMC_IVM = [
        [
            ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
            ModFn(DTKConfigBuilder.set_param, 'Simulation_Duration', smc_start_day + 365),
            ModFn(add_smc_group,
                  start_days=[smc_start_day],
                  coverage=smc_coverage, drug=drug, agemax=agemax,),
            ModFn(add_ivermectin_group, coverage=smc_coverage, start_days=[smc_start_day], drug_code=drug_code,
                  smc_code=drug + str(agemax), target_group=target_group),
            ModFn(add_summary_report, start_day=smc_start_day),
            ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path',
                  '{path}/output'.format(path=
                                         path))
        ]
        for smc_start_day in intervention_days
        for smc_coverage in coverages
        for seed in range(num_seeds)
        for path in serialized_file_output
        for agemax in [5]
        for target_group in target_groups_IVM
        for drug in ['DP']
        for drug_code in [14, 30, 60, 90]
    ]

    expanded_SMC_IVM = [
        [
            ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
            ModFn(DTKConfigBuilder.set_param, 'Simulation_Duration', smc_start_day + 365),
            ModFn(add_smc_group,
                  start_days=[smc_start_day],
                  coverage=smc_coverage, drug=drug, agemax=agemax,),
            ModFn(add_ivermectin_group, agemax=agemax, coverage=smc_coverage, start_days=[smc_start_day], drug_code=drug_code,
                  smc_code=drug + str(agemax), target_group=target_group),
            ModFn(add_summary_report, start_day=smc_start_day),
            ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path',
                  '{path}/output'.format(path=
                                         path))
        ]
        for smc_start_day in intervention_days
        for smc_coverage in coverages
        for seed in range(num_seeds)
        for path in serialized_file_output
        for agemax in [10, 15]
        for target_group in expanded_SMC_groups_IVM
        for drug in ['DP']
        for drug_code in [14, 30, 60, 90]
    ]

    SMC = [
        [
            ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
            ModFn(DTKConfigBuilder.set_param, 'Simulation_Duration', smc_start_day + 365),
            ModFn(add_smc_group,
                  start_days=[smc_start_day],
                  coverage=smc_coverage, drug=drug, agemax=agemax,),
            ModFn(add_summary_report, start_day=smc_start_day),
            ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path',
                  '{path}/output'.format(path=
                                         path))
        ]
        for smc_start_day in intervention_days
        for smc_coverage in coverages
        for seed in range(num_seeds)
        for path in serialized_file_output
        for agemax in [5, 10]
        for drug in ['DP']
    ]

    builder = ModBuilder.from_list(expanded_SMC_IVM)

    set_climate_constant(cb)
    set_species_param(cb, 'gambiae', 'Larval_Habitat_Types',
                      {"LINEAR_SPLINE": {
                          "Capacity_Distribution_Per_Year": {
                              "Times": [0.0, 30.417, 60.833, 91.25, 121.667, 152.083,
                                        182.5, 212.917, 243.333, 273.75, 304.167, 334.583],
                              "Values": [3, 0.8,  1.25, 0.1,  2.7, 8, 4, 35, 6.8,  6.5, 2.6, 2.1]
                          },
                          "Max_Larval_Capacity": 1e9
                      }})
    set_species_param(cb, "gambiae", "Indoor_Feeding_Fraction", 0.9)
    set_species_param(cb, "gambiae", "Adult_Life_Expectancy", 20)

    cb.update_params({"Demographics_Filenames": ['Interventions/single_node_demographics.json',
                                   'Interventions/single_node_drugs_overlay.json'],

                      'Antigen_Switch_Rate': pow(10, -9.116590124),
                      'Base_Gametocyte_Production_Rate': 0.06150582,
                      'Base_Gametocyte_Mosquito_Survival_Rate': 0.002011099,
                      'Base_Population_Scale_Factor': 0.1,

                      "Birth_Rate_Dependence": "FIXED_BIRTH_RATE",
                      "Death_Rate_Dependence": "NONDISEASE_MORTALITY_BY_AGE_AND_GENDER",
                      "Disable_IP_Whitelist": 1,

                      "Enable_Vital_Dynamics": 1,
                      "Enable_Birth": 1,
                      'Enable_Default_Reporting': 1,
                      'Enable_Demographics_Other': 1,
                      'Enable_Property_Output': 1,

                      'Falciparum_MSP_Variants': 32,
                      'Falciparum_Nonspecific_Types': 76,
                      'Falciparum_PfEMP1_Variants': 1070,
                      'Gametocyte_Stage_Survival_Rate': 0.588569307,

                      'MSP1_Merozoite_Kill_Fraction': 0.511735322,
                      'Max_Individual_Infections': 3,
                      'Nonspecific_Antigenicity_Factor': 0.415111634,

                      'x_Temporary_Larval_Habitat': 1,
                      'logLevel_default': 'ERROR',

                      "Simulation_Duration": sim_duration,
                      "Serialization_Test_Cycles": 0,
                      'Serialized_Population_Filenames': ['state-21900.dtk'],
                      "Vector_Species_Names": ['gambiae']
                      })

    run_sim_args = {'config_builder': cb,
                    'exp_name': expname,
                    'exp_builder': builder}

    exp_manager = ExperimentManagerFactory.from_setup()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)
