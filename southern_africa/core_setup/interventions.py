
import math
import os
import pandas as pd
import numpy as np
from dtk.interventions.input_EIR import add_InputEIR
from dtk.interventions.ivermectin import add_ivermectin
from dtk.interventions.migrate_to import add_migration_event

from helpers.windows_filesystem import get_dropbox_location
from helpers.relative_time import convert_to_day_365

from dtk.interventions.irs import add_IRS
from dtk.interventions.itn_age_season import add_ITN_age_season
from malaria.interventions.health_seeking import add_health_seeking
from malaria.interventions.malaria_drug_campaigns import add_drug_campaign

from gridded_sims.run.build_cb import get_project_folder
from gridded_sims.run.site import find_cells_for_this_catchment

project_folder = get_project_folder()
interventions_folder = os.path.join(project_folder, "dtk_simulation_input/mozambique/")


def add_regional_EIR_node(cb):
    regional_EIR_node_label = 100000
    monthly_profile = 0.5*np.array([12.1, 23.9, 33.5, 14.8, 6.8, 4.9, 3.3, 3.8, 3.5, 3.3, 2.8, 4.3])
    add_InputEIR(cb,
                 monthlyEIRs=list(monthly_profile),
                 nodeIDs=[regional_EIR_node_label],
                 start_day=0)





# Input: intervention_df = pd.read_csv("grid_itn_events.csv")
# format: interventions by grid cell.  Has columns for "grid_cell", "fulldate", "simday"
def try_campaign_compression(intervention_df, field_list, bin_fidelity=0.1, simday_fidelity=10):
    def round_nearest(x, a):
        rounded = round(round(x / a) * a, -int(math.floor(math.log10(a))))
        return rounded

    binned_intervention_df = intervention_df.copy(deep=True)

    # Bin these fields with bin width = bin_fidelity
    for field in field_list:
        binned_intervention_df[field] = intervention_df[field].map(lambda x: round_nearest(x, bin_fidelity))

    binned_intervention_df['simday'] = intervention_df['simday'].map(lambda x: round_nearest(x, simday_fidelity))

    # Group by the new binned fields, as well as by date
    data_fields = ['simday'] + field_list
    binned_and_grouped = binned_intervention_df.groupby(data_fields)

    return [binned_intervention_df, binned_and_grouped, data_fields]


# All of the following functions add_X require input of a dataframe with the following format:
# simday: day of simulation to implement the intervention
# grid_cell: node ID to implement intervention on
# other fields, intervention-dependent

def add_hs(cb, events_df):
    hs_field_list = ["cov_newclin_youth", "cov_newclin_adult", "cov_severe_youth", "cov_severe_adult", "duration"]
    binned_intervene_events, binned_and_grouped, data_fields = try_campaign_compression(events_df, hs_field_list)

    for table, group in binned_and_grouped:
        table_dict = dict(zip((data_fields), table))

        node_list = sorted(group['grid_cell'])
        node_dict = {"class": "NodeSetNodeList", "Node_List": node_list}

        add_health_seeking(cb,
                           start_day=table_dict['simday'],
                           targets=[{'trigger': 'NewClinicalCase',
                                     'coverage': float(table_dict['cov_newclin_youth']),
                                     'agemin': 0,
                                     'agemax': 5,
                                     'seek': 1,
                                     'rate': 0.3},
                                    {'trigger': 'NewClinicalCase',
                                     'coverage': float(table_dict['cov_newclin_adult']),
                                     'agemin': 5,
                                     'agemax': 100,
                                     'seek': 1,
                                     'rate': 0.3},
                                    {'trigger': 'NewSevereCase',
                                     'coverage': float(table_dict['cov_severe_youth']),
                                     'agemin': 0,
                                     'agemax': 5,
                                     'seek': 1,
                                     'rate': 0.5},
                                    {'trigger': 'NewSevereCase',
                                     'coverage': float(table_dict['cov_severe_adult']),
                                     'agemin': 5,
                                     'agemax': 100,
                                     'seek': 1,
                                     'rate': 0.5}],
                           drug=['Artemether', 'Lumefantrine'],
                           dosing='FullTreatmentNewDetectionTech',
                           nodes=node_list,
                           duration=table_dict['duration'])



def add_itn(cb, events_df):
    itn_field_list = ["age_cov", "cov_all", "min_season_cov", "fast_fraction"]
    binned_intervene_events, binned_and_grouped, data_fields = try_campaign_compression(events_df, itn_field_list)

    birthnet_df = events_df.copy(deep=True)
    birthnet_df.sort_values(by='simday', inplace=True)
    birthnet_df['duration'] = birthnet_df.groupby('grid_cell')['simday'].shift(-1).sub(birthnet_df['simday'])
    birthnet_df['duration'].fillna(-1, inplace=True)
    birthnet_field_list = itn_field_list + ["duration"]
    BIRTH_binned_intervene_events, BIRTH_binned_and_grouped, BIRTH_data_fields = try_campaign_compression(birthnet_df, birthnet_field_list)

    for table, group in binned_and_grouped:
        table_dict = dict(zip((data_fields), table))
        node_list = sorted(list(set(group['grid_cell']))) # fixme Needed to add this because sometimes there are duplicate nodes in the list, and this breaks things

        start = float(table_dict['simday'])
        if start >= 0:
            # Regular bednet distribution
            add_ITN_age_season(cb,
                               start=float(table_dict['simday']),
                               age_dep={'youth_cov': float(table_dict['age_cov']),
                                        'youth_min_age': 5,
                                        'youth_max_age': 20},
                               coverage_all=float(table_dict['cov_all']),
                               as_birth=False,
                               seasonal_dep={'min_cov': float(table_dict['min_season_cov']),
                                             'max_day': 60},
                               discard={'halflife1': 260,
                                        'halflife2': 2106,
                                        'fraction1': float(table_dict['fast_fraction'])},
                               nodeIDs=node_list)

    # Birthnet distribution
    for table, group in BIRTH_binned_and_grouped:
        table_dict = dict(zip((BIRTH_data_fields), table))
        node_list = sorted(group['grid_cell'])

        start = float(table_dict['simday'])
        if start >= 0:
            add_ITN_age_season(cb,
                               as_birth=True,
                               duration=table_dict['duration'],
                               start=float(table_dict['simday']),
                               age_dep={'youth_cov': float(table_dict['age_cov']),
                                        'youth_min_age': 5,
                                        'youth_max_age': 20},
                               coverage_all=float(table_dict['cov_all']),
                               seasonal_dep={'min_cov': float(table_dict['min_season_cov']),
                                             'max_day': 60},
                               discard={'halflife1': 260,
                                        'halflife2': 2106,
                                        'fraction1': float(table_dict['fast_fraction'])},
                               nodeIDs=node_list)



def add_irs(cb, events_df):
    irs_field_list = ["cov_all", "killing", "exp_duration", "box_duration"]
    binned_intervene_events, binned_and_grouped, data_fields = try_campaign_compression(events_df, irs_field_list)

    for table, group in binned_and_grouped:
        table_dict = dict(zip((data_fields), table))
        node_list = sorted(group['grid_cell'])

        add_IRS(cb, start=int(table_dict['simday']),
                coverage_by_ages=[{'coverage': float(table_dict['cov_all'])}],
                waning={"Killing_Config": {
                    "Initial_Effect": float(table_dict['killing']),
                    "Decay_Time_Constant": float(table_dict['exp_duration']),
                    "Box_Duration": float(table_dict['box_duration']),
                    "class": "WaningEffectBoxExponential"
                }},
                nodeIDs=node_list)


def add_mda(cb, events_df):
    mda_field_list = ["cov_all"]
    binned_intervene_events, binned_and_grouped, data_fields = try_campaign_compression(events_df, mda_field_list)

    for table, group in binned_and_grouped:
        table_dict = dict(zip((data_fields), table))
        node_list = sorted(group['grid_cell'])

        add_drug_campaign(cb,
                          campaign_type='MDA',
                          drug_code='DP',
                          start_days=[float(table_dict['simday'])],
                          coverage=table_dict['cov_all'],
                          repetitions=1,
                          target_residents_only=False,
                          nodeIDs=node_list)


def add_pure_ivm_mda(cb, events_df, ivm_duration, restricted_ivm=False):
    mda_field_list = ["cov_all"]
    binned_intervene_events, binned_and_grouped, data_fields = try_campaign_compression(events_df, mda_field_list)

    for table, group in binned_and_grouped:
        table_dict = dict(zip((data_fields), table))
        node_list = sorted(group['grid_cell'])

        if not restricted_ivm:
            add_ivermectin(cb,
                           box_duration=ivm_duration,
                           coverage=table_dict['cov_all'],
                           start_days = [float(table_dict['simday'])],
                           target_residents_only=False,
                           nodeIDs=node_list)
        else:
            target_groups = [{'agemin': 5, 'agemax': 200, 'gender': 'Male'},
                             {'agemin': 51, 'agemax': 200, 'gender': 'Female'},
                             {'agemin': 5, 'agemax': 12, 'gender': 'Female'}]

            for tg in target_groups:
                add_ivermectin(cb,
                               box_duration=ivm_duration,
                               coverage=table_dict['cov_all'],
                               start_days=[float(table_dict['simday'])],
                               target_residents_only=False,
                               nodeIDs=node_list,
                               target_group=tg)


def add_rcd(cb, events_df):
    rcd_field_list = ["coverage", "trigger_coverage", "interval"]
    binned_intervene_events, binned_and_grouped, data_fields = try_campaign_compression(events_df, rcd_field_list)

    for table, group in binned_and_grouped:
        table_dict = dict(zip((data_fields), table))
        node_list = sorted(group['grid_cell'])

        add_drug_campaign(cb,
                          campaign_type='rfMSAT',
                          drug_code='AL',
                          diagnostic_type='BLOOD_SMEAR_PARASITES',
                          diagnostic_threshold=0,
                          start_days=[float(table_dict['simday'])],
                          coverage=float(table_dict['coverage']),
                          trigger_coverage=float(table_dict['trigger_coverage']),
                          nodeIDs=node_list)


def add_travellers(cb, catch):
    all_nodes = list(find_cells_for_this_catchment(catch))

    add_migration_event(cb,
                        nodeto=100000,
                        coverage=0.5,
                        repetitions=500,
                        tsteps_btwn=30,
                        duration_of_stay=3,
                        duration_before_leaving_distr_type='UNIFORM_DURATION',
                        duration_before_leaving=0,
                        duration_before_leaving_2=30,
                        nodesfrom={'class': 'NodeSetNodeList',
                                   'Node_List': all_nodes},
                        ind_property_restrictions=[{'TravelerStatus': 'IsTraveler'}])



hs_default = os.path.join(interventions_folder, "grid_all_healthseek_events.csv")
itn_default = os.path.join(interventions_folder, "grid_all_itn_events.csv")
irs_default = os.path.join(interventions_folder, "grid_all_irs_events.csv")
mda_default = os.path.join(interventions_folder, "grid_all_mda_events.csv")
rcd_default = os.path.join(interventions_folder, "grid_all_react_events_condensed.csv")
print("Note: using condensed RCD file.")
def preload_intervention_csvs(catch,
                              sim_start_date,
                              hs_file=hs_default,
                              itn_file=itn_default,
                              irs_file=irs_default,
                              mda_file=mda_default,
                              rcd_file=rcd_default):
    hs_df = pd.read_csv(hs_file)
    itn_df = pd.read_csv(itn_file)
    irs_df = pd.read_csv(irs_file)
    mda_df = pd.read_csv(mda_file)
    rcd_df = pd.read_csv(rcd_file)

    df_dict = {"hs": hs_df,
               "itn": itn_df,
               "irs": irs_df,
               "mda": mda_df,
               "rcd": rcd_df}

    # Add simday column for adding to campaign file
    for key in list(df_dict.keys()):
        df = df_dict[key]

        # Add column which specifies which simulation day the intervention of that row is to be implemented
        df['simday'] = [convert_to_day_365(x, sim_start_date, "%Y-%m-%d") for x in df.fulldate]

    # Restrict to catchment of interest
    catch_cells = find_cells_for_this_catchment(catch)
    for key in list(df_dict.keys()):
        df = df_dict[key]
        df = df[np.in1d(df["grid_cell"], catch_cells)]
        df.reset_index(inplace=True)
        df_dict[key] = df

    return df_dict



def add_intervention_combos(cb, intervention_df_dict, catch, itn, irs, mda, rcd, hs=True, regional_EIR_node=True, travellers=True):
    if hs:
        add_hs(cb, intervention_df_dict["hs"])
    if regional_EIR_node:
        add_regional_EIR_node(cb)
    if travellers:
        add_travellers(cb, catch)

    if itn:
        add_itn(cb, intervention_df_dict["itn"])
    if irs:
        add_irs(cb, intervention_df_dict["irs"])
    if mda:
        add_mda(cb, intervention_df_dict["mda"])
    if rcd:
        add_rcd(cb, intervention_df_dict["rcd"])

    return {"itn": itn,
            "irs": irs,
            "mda": mda,
            "rcd": rcd}

def add_all_interventions(cb, catch, sim_start_date):
    intervention_df_dict = preload_intervention_csvs(catch, sim_start_date)
    return add_intervention_combos(cb, intervention_df_dict, catch, True, True, True, True)