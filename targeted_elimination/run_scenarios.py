import os
import pandas as pd
import numpy as np

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.utils.reports.VectorReport import add_human_migration_report
from malaria.reports.MalariaReport import add_event_counter_report, add_filtered_spatial_report
from simtools.ModBuilder import ModBuilder, ModFn
from dtk.vector.species import update_species_param
from malaria.interventions.health_seeking import add_health_seeking
from malaria.interventions.malaria_drug_campaigns import add_drug_campaign
from dtk.interventions.ivermectin import add_ivermectin
from simtools.Analysis.BaseAnalyzers.SimulationDirectoryMapAnalyzer import SimulationDirectoryMapAnalyzer
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.Analysis.AnalyzeManager import AnalyzeManager

from configure_forest_system import configure_forest_system


# General
exp_name = 'Forest_endectocides'
village_nodes = [1,2]
forest_nodes = [3]
years = 4  # length of simulation, in years
report_migration = False

burnin_expid = 'bb20cb63-2a1f-e911-a2be-c4346bcb1554'

interventions = [
    'ivermectin_to_forest',
    'drug_MDA',
    'drug_to_forest',
    'drug_ivermectin_MDA',
    'forest_HS',
    'ivermectin_MDA'
                 ]

ivm_durations = [14, 30, 60, 90] #range(7, 92, 7)

num_seeds = 100
starting_seed = 0

# Serialization
serialize = False  # If true, save serialized files
pull_from_serialization = True  # requires serialization date and experiment id


def create_sim_map(burnin_id) :

    am = AnalyzeManager(burnin_id, analyzers=SimulationDirectoryMapAnalyzer())
    am.analyze()

    sim_map_dict = am.analyzers[0].results
    df = pd.concat([pd.DataFrame(exp) for exp_id, exp in sim_map_dict.items()])
    return df


def update_serialization_params(cb):

    serialization_date = 100 * 365
    cb.update_params({'Serialized_Population_Filenames':
                     ['state-%05d.dtk' % serialization_date],
    })

    return {}


def configure_health_seeking(cb, pull_from_serialization, village_nodes, set_forest_HS=True) :

    if pull_from_serialization :
        add_health_seeking(cb, start_day=0,
                           drug=['Artemether', 'Lumefantrine'],
                           nodes=village_nodes,
                           targets=[
                               {'trigger': 'NewClinicalCase', 'coverage': 0.8, 'agemin': 0, 'agemax': 10, 'seek': 1,
                                'rate': 0.3},
                               {'trigger': 'NewClinicalCase', 'coverage': 0.7, 'agemin': 10, 'agemax': 100, 'seek': 1,
                                'rate': 0.3},
                               {'trigger': 'NewSevereCase', 'coverage': 0.95, 'agemin': 0, 'agemax': 100, 'seek': 1,
                                'rate': 0.3}],
                           ind_property_restrictions=[{'ForestGoing' : 'HatesForest'}]
                           )
        if set_forest_HS :
            add_health_seeking(cb, start_day=0,
                               drug=['Artemether', 'Lumefantrine'],
                               nodes=village_nodes,
                               targets=[
                                   {'trigger': 'NewClinicalCase', 'coverage': 0.4, 'seek': 1, 'rate': 0.2}],
                               ind_property_restrictions=[{'ForestGoing': 'LovesForest'}]
                               )
    else :
        add_health_seeking(cb, start_day=0,
                           drug=['Artemether', 'Lumefantrine'],
                           nodes=village_nodes,
                           targets=[
                               {'trigger': 'NewClinicalCase', 'coverage': 0.4, 'agemin': 0, 'agemax': 5, 'seek': 1,
                                'rate': 0.3},
                               {'trigger': 'NewClinicalCase', 'coverage': 0.2, 'agemin': 5, 'agemax': 100, 'seek': 1,
                                'rate': 0.3},
                               {'trigger': 'NewSevereCase', 'coverage': 0.5, 'agemin': 0, 'agemax': 100, 'seek': 1,
                                'rate': 0.3}]
                           )


def add_intervention(cb, intervention, coverage, target, ivm_duration) :

    IPR = [{'HomeVillage' : 'Village1'}] if target == 'vill1' else None

    if intervention == 'ivermectin_to_forest' :
        add_ivermectin(cb, box_duration=ivm_duration, coverage=coverage, start_days=[365],
                       trigger_condition_list=['Immigrating'],
                       nodeIDs=forest_nodes,
                       target_group={'agemin': 5, 'agemax': 200, 'gender': 'Male'},
                       target_residents_only=False,
                       ind_property_restrictions=IPR
                       )
    elif intervention == 'drug_to_forest' :
        add_drug_campaign(cb, 'MDA', 'DP', coverage=coverage, start_days=[365],
                          repetitions=1,
                          trigger_condition_list=['Immigrating'],
                          nodeIDs=forest_nodes,
                          target_residents_only=False,
                          ind_property_restrictions=IPR
                          )
    elif intervention == 'drug_MDA' :
        add_drug_campaign(cb, 'MDA', 'DP', coverage=coverage, start_days=[365 + 152],
                          repetitions=3, tsteps_btwn_repetitions=30,
                          nodeIDs=village_nodes,
                          target_residents_only=False,
                          ind_property_restrictions=IPR
                          )
    elif intervention == 'drug_ivermectin_MDA' :
        add_drug_campaign(cb, 'MDA', 'DP', coverage=coverage, start_days=[365 + 152],
                          repetitions=3, tsteps_btwn_repetitions=30,
                          nodeIDs=village_nodes,
                          target_residents_only=False,
                          ind_property_restrictions=IPR
                          )
        add_ivermectin(cb, box_duration=ivm_duration, coverage=1, start_days=[365],
                       trigger_condition_list=['Received_Campaign_Drugs'],
                       target_group={'agemin': 5, 'agemax': 200, 'gender': 'Male'},
                       nodeIDs=village_nodes,
                       target_residents_only=False
                       )
        add_ivermectin(cb, box_duration=ivm_duration, coverage=1, start_days=[365],
                       trigger_condition_list=['Received_Campaign_Drugs'],
                       target_group={'agemin': 51, 'agemax': 200, 'gender': 'Female'},
                       nodeIDs=village_nodes,
                       target_residents_only=False
                       )
        add_ivermectin(cb, box_duration=ivm_duration, coverage=1, start_days=[365],
                       trigger_condition_list=['Received_Campaign_Drugs'],
                       target_group={'agemin': 5, 'agemax': 12, 'gender': 'Female'},
                       nodeIDs=village_nodes,
                       target_residents_only=False
                       )
    elif intervention == 'ivermectin_MDA' :
        add_ivermectin(cb, box_duration=ivm_duration, coverage=coverage,
                       start_days=[365 + 152 + 30 * x for x in range(3)],
                       nodeIDs=village_nodes,
                       target_group={'agemin': 5, 'agemax': 200, 'gender': 'Male'},
                       target_residents_only=False,
                       ind_property_restrictions=IPR
                       )
        add_ivermectin(cb, box_duration=ivm_duration, coverage=coverage,
                       start_days=[365 + 152 + 30 * x for x in range(3)],
                       nodeIDs=village_nodes,
                       target_group={'agemin': 51, 'agemax': 200, 'gender': 'Female'},
                       target_residents_only=False,
                       ind_property_restrictions=IPR
                       )
        add_ivermectin(cb, box_duration=ivm_duration, coverage=coverage,
                       start_days=[365 + 152 + 30 * x for x in range(3)],
                       nodeIDs=village_nodes,
                       target_group={'agemin': 5, 'agemax': 12, 'gender': 'Female'},
                       target_residents_only=False,
                       ind_property_restrictions=IPR
                       )

    forest_HS_coverage = coverage if intervention == 'forest_HS' else 0.4
    add_health_seeking(cb, start_day=0,
                       drug=['Artemether', 'Lumefantrine'],
                       nodes=village_nodes,
                       targets=[
                           {'trigger': 'NewClinicalCase', 'coverage': forest_HS_coverage, 'seek': 1, 'rate': 0.2}],
                       ind_property_restrictions=[{'ForestGoing': 'LovesForest'}]
                       )

    event_trigger_list=['Received_Treatment']
    if 'ivermectin' in intervention :
        event_trigger_list.append('Received_Ivermectin')
    if 'drug' in intervention :
        event_trigger_list.append('Received_Campaign_Drugs')
    add_event_counter_report(cb, event_trigger_list=event_trigger_list)

    return { 'intervention' : intervention,
             'coverage' : coverage,
             'target' : target,
             'IVM duration' : ivm_duration}


def sample_anthrop_and_outdoor(cb, species, anthropophily, outdoor) :

    update_species_param(cb, species, 'Anthropophily', anthropophily)
    update_species_param(cb, species, 'Indoor_Feeding_Fraction', 1-outdoor)

    return { '%s_Anthropophily' % species : anthropophily,
             '%s_Outdoor_Fraction' % species : outdoor}


if __name__=="__main__":

    SetupParser.init('HPC')
    cb = configure_forest_system('', years, add_outbreak=False)
    cb.update_params({
        'Config_Name': exp_name,
    })

    cb.set_exe_collection("520eaf08-cbf1-e811-a2bd-c4346bcb1555")
    cb.set_dll_collection("53887114-1ef4-e811-a2bd-c4346bcb1555")
    cb.set_input_collection("dc81f29a-e773-e811-a2c0-c4346bcb7275")

    if serialize:
        cb.update_params({'Serialization_Time_Steps': [365*years]})

    if pull_from_serialization :
        set_forest_HS = False
        burnin_df = create_sim_map(burnin_expid)
        burnin_df = burnin_df[(burnin_df['dirus_Anthropophily'] == 0.5) &
                              (burnin_df['x_Temporary_Larval_Habitat'] == 2.4)]
        update_serialization_params(cb)
        builder = ModBuilder.from_list([[
            ModFn(add_intervention, intervention, coverage, target, ivm_duration),
            ModFn(sample_anthrop_and_outdoor, 'dirus', row['dirus_Anthropophily'], row['dirus_Outdoor_Fraction']),
            ModFn(DTKConfigBuilder.update_params, {
                'x_Temporary_Larval_Habitat': row['x_Temporary_Larval_Habitat'],
                'Run_Number': y,
                'Serialized_Population_Path': os.path.join(row['path'], 'output')
            })
            ]
            for y in range(starting_seed, starting_seed+num_seeds)
            for coverage in np.linspace(0,1,6)
            for intervention in interventions
            for target in ['bothvill', 'vill1']
            for ivm_duration in ivm_durations
            for r,row in burnin_df.iterrows()
        ])

    else :
        set_forest_HS = True
        builder = ModBuilder.from_list([[
            ModFn(sample_anthrop_and_outdoor, 'dirus', anthropophily, outdoor),
            ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat', x_LH),
            ModFn(DTKConfigBuilder.set_param, 'Run_Number', y)]
            for y in range(num_seeds) for anthropophily in np.linspace(0,1,11) for outdoor in [0.99]
            for x_LH in [0.6, 1.2, 2.4]
        ])

    configure_health_seeking(cb, pull_from_serialization, village_nodes, set_forest_HS=set_forest_HS)

    # ---- CUSTOM REPORTS ----
    add_filtered_spatial_report(cb, start=365*(years-1), end=365*years,
                                channels=['Population', 'Prevalence',
                                          # 'Adult_Vectors',
                                          # 'PfHRP2_Prevalence',
                                          # 'PCR_Prevalence',
                                          'New_Infections'])
    if report_migration:
        cb.update_params({
            "Report_Event_Recorder": 1,
            "Report_Event_Recorder_Events": ["Immigrating", "Emigrating", 'Births'],
            "Report_Event_Recorder_Ignore_Events_In_List": 0
        })

        add_human_migration_report(cb)

    # Run args
    run_sim_args = {'config_builder': cb,
                    'exp_name': exp_name,
                    'exp_builder': builder}

    em = ExperimentManagerFactory.from_cb(cb)
    em.run_simulations(**run_sim_args)
