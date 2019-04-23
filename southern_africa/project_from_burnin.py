from dtk.interventions.ivermectin import add_ivermectin
from malaria.reports.MalariaReport import add_filtered_report, add_filtered_spatial_report, add_event_counter_report
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment, retrieve_simulation

from core_setup.simplified_ento import catch_1_yr_spline

SetupParser.default_block = "HPC"
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder

from malaria.interventions.malaria_drug_campaigns import add_drug_campaign

from core_setup.build_cb import *
from core_setup.interventions import add_all_interventions, preload_intervention_csvs, add_intervention_combos, \
    add_mda, add_hs, add_itn, add_irs, add_pure_ivm_mda
from core_setup.reports import add_all_reports


from gridded_sims.calib.calib_helpers import set_ento

from COMPS import Client
Client.login("comps.idmod.org")


def draw_from_burnin(cb, burnin_sim_id):
    sim = retrieve_simulation(burnin_sim_id)
    serialized_file_path = sim.get_path()

    cb.update_params({
        'Serialized_Population_Path': os.path.join(serialized_file_path, 'output'),
        'Serialized_Population_Filenames': ['state-19710-000.dtk', 'state-19710-001.dtk']
    })



def iver_sweep(cb, intervention_df_dict, vc_pack, drug_pack, vc_coverage, mda_coverage, ivm_duration, restricted_ivm=False):
    add_hs(cb, intervention_df_dict["hs"])
    event_trigger_list = ["Received_Treatment"]

    if vc_pack == "none":
        pass
    elif vc_pack == "ITN only":
        itn_events = intervention_df_dict["itn"]
        itn_events["cov_all"] = vc_coverage
        add_itn(cb, itn_events)
    elif vc_pack == "IRS only":
        irs_events = intervention_df_dict["irs"]
        irs_events["cov_all"] = vc_coverage
        add_irs(cb, irs_events)
    elif vc_pack == "ITN and IRS":
        itn_events = intervention_df_dict["itn"]
        itn_events["cov_all"] = vc_coverage
        irs_events = intervention_df_dict["irs"]
        irs_events["cov_all"] = vc_coverage
        add_itn(cb, itn_events)
        add_irs(cb, irs_events)

    # Is there a DP component to the MDA?
    if mda_coverage > 0:
        mda_events = intervention_df_dict["mda"]
        mda_events["cov_all"] = mda_coverage

        if drug_pack != "Pure IVM MDA": # MDA includes DP
            event_trigger_list.append('Received_Campaign_Drugs')
            add_mda(cb, mda_events)
        elif drug_pack == "Pure IVM MDA": # MDA with only IVM and no DP
            event_trigger_list.append('Received_Ivermectin')
            add_pure_ivm_mda(cb, mda_events, ivm_duration, restricted_ivm=restricted_ivm)

    if drug_pack == "MDA with IVM":
        event_trigger_list.append('Received_Ivermectin')
        if not restricted_ivm:
            add_ivermectin(cb, box_duration=ivm_duration, coverage=1.0, start_days=[1],
                           trigger_condition_list=['Received_Campaign_Drugs'], target_residents_only=False)
        else:
            target_groups = [{'agemin': 5, 'agemax': 200, 'gender': 'Male'},
                             {'agemin': 51, 'agemax': 200, 'gender': 'Female'},
                             {'agemin': 5, 'agemax': 12, 'gender': 'Female'}]

            for tg in target_groups:
                add_ivermectin(cb, box_duration=ivm_duration, coverage=1.0, start_days=[1],
                               trigger_condition_list=['Received_Campaign_Drugs'], target_residents_only=False,
                               target_group=tg)

    add_event_counter_report(cb, event_trigger_list=event_trigger_list)
    return {"vc_pack": vc_pack,
            "drug_pack": drug_pack,
            "vc_coverage": vc_coverage,
            "mda_coverage": mda_coverage,
            "ivm_duration": ivm_duration}



def add_filtered_reports(cb, catch, start=0, duration=10000000):
    add_filtered_spatial_report(cb, channels=["Population", "True_Prevalence"], start=start, end=(start + duration))

    # Filtered report just for work node, and just for catchment:
    regional_EIR_node_label = 100000
    catch_node_ids = find_cells_for_this_catchment(catch)

    add_filtered_report(cb, nodes=[regional_EIR_node_label], description='Work', start=start,
                        end=(start + duration))
    add_filtered_report(cb, nodes=catch_node_ids, description='Catchment', start=start,
                        end=(start + duration))

def generate_input_variations():
    vc_packs = ["none", "ITN only", "IRS only", "ITN and IRS"]
    drug_packs = ["none","Pure IVM MDA","MDA without IVM", "MDA with IVM"]
    mda_coverages = [0.2,0.4,0.6,0.8,1.0]
    ivm_durations = [14, 30, 60, 90]

    tuples_list = []
    for vc_pack in vc_packs:
        for drug_pack in drug_packs:
            for mda_coverage in mda_coverages:
                for ivm_duration in ivm_durations:
                    if drug_pack == "none" and mda_coverage == 0.2 and ivm_duration == 14:
                        tuples_list.append((vc_pack, drug_pack, -1, -1))
                    elif drug_pack == "MDA without IVM" and ivm_duration == 14:
                        tuples_list.append((vc_pack, drug_pack, mda_coverage, -1))
                    elif drug_pack == "MDA with IVM":
                        tuples_list.append((vc_pack, drug_pack, mda_coverage, ivm_duration))

    return tuples_list


if __name__ == "__main__":
    # Run parameters:
    catch = "Magude-Sede-Facazissa"
    num_cores = 2
    num_seeds = 100
    priority = ""
    coreset = ""
    burnin_sim_id = "" 


    experiment_name = "MSF_iver_all_restricted"
    arab_times, arab_spline = catch_1_yr_spline(catch, "gambiae")
    funest_times, funest_spline = catch_1_yr_spline(catch, "funestus")

    cb = build_project_cb()
    catchment_cb_params(cb, catch)

    # Burnin serialized at 2009.  So current run, which goes to beginning of 2019, has run duration of 10 years:
    # sim_start_date = "2009-01-01"
    sim_start_date = "2014-01-01"
    cb.set_param("Simulation_Duration", 4 * 365)
    draw_from_burnin(cb, burnin_sim_id)

    # NO IMPORTATIONS
    cb.set_param("x_Regional_Migration", 0)

    ivm_interventions_folder = os.path.join(project_folder, "dtk_simulation_input/mozambique/ivermectin_paper/")
    hs_timeshifted = os.path.join(ivm_interventions_folder, "grid_all_healthseek_events.csv")
    itn_simplified = os.path.join(ivm_interventions_folder, "grid_all_itn_events.csv")
    irs_simplified = os.path.join(ivm_interventions_folder, "grid_all_irs_events.csv")
    mda_simplified = os.path.join(ivm_interventions_folder, "grid_all_mda_events.csv")
    intervention_df_dict = preload_intervention_csvs(catch, sim_start_date,
                                                     hs_file=hs_timeshifted,
                                                     itn_file=itn_simplified,
                                                     irs_file=irs_simplified,
                                                     mda_file=mda_simplified)
    # add_intervention_combos(cb, intervention_df_dict, catch, True, False, False, False, travellers=False)
    add_all_reports(cb, catch, start=0 * 365)

    arab_burnin = 9.388
    funest_burnin = 10.2778
    set_ento(cb, arab_burnin, funest_burnin, arab_times, arab_spline, funest_times, funest_spline)



    SetupParser.init()
    SetupParser.set("HPC", "priority", priority)
    SetupParser.set("HPC", "node_group", coreset)

    modlists = []

    new_modlist = [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed) for seed in range(num_seeds)]
    modlists.append(new_modlist)

    # Generate argument list:
    input_variations = generate_input_variations()
    print ("{} variations (100 seeds each)".format(len(input_variations)))

    new_modlist = [ModFn(iver_sweep, intervention_df_dict, x[0], x[1], x[2], x[3], x[4], False) for x in input_variations]

    modlists.append(new_modlist)

    builder = ModBuilder.from_combos(*modlists)

    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(config_builder=cb,
                                exp_name=experiment_name,
                                exp_builder=builder)