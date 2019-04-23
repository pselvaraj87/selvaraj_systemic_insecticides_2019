from dtk.interventions.ivermectin import add_ivermectin
from malaria.reports.MalariaReport import add_filtered_report, add_filtered_spatial_report, add_event_counter_report
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment

SetupParser.default_block = "HPC"
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder

from malaria.interventions.malaria_drug_campaigns import add_drug_campaign

from core_setup.build_cb import *
from core_setup.interventions import add_all_interventions, preload_intervention_csvs, add_intervention_combos, add_mda
from core_setup.reports import add_all_reports


from core_setup.calib_helpers import set_ento

from core_setup.simplified_ento import catch_1_yr_spline

def add_filtered_reports(cb, catch, start=0, duration=10000000):
    add_filtered_spatial_report(cb, channels=["Population", "True_Prevalence"], start=start, end=(start + duration))

    # Filtered report just for work node, and just for catchment:
    regional_EIR_node_label = 100000
    catch_node_ids = find_cells_for_this_catchment(catch)

    add_filtered_report(cb, nodes=[regional_EIR_node_label], description='Work', start=start,
                        end=(start + duration))
    add_filtered_report(cb, nodes=catch_node_ids, description='Catchment', start=start,
                        end=(start + duration))


if __name__ == "__main__":
    # Run parameters:
    catch = "Magude-Sede-Facazissa"
    num_cores = 2
    num_seeds = 1
    # priority = "Normal"
    # coreset = "emod_abcd"
    priority = "Highest"
    coreset = "emod_32cores"
    # ====================================================================================================================

    experiment_name = "MSF_iver_twocore_burnin"
    arab_times, arab_spline = catch_1_yr_spline(catch, "gambiae")
    funest_times, funest_spline = catch_1_yr_spline(catch, "funestus")

    cb = build_project_cb(num_cores=num_cores)
    catchment_cb_params(cb, catch)
    # cb.set_param("x_Regional_Migration", 0.0081)

    # Burnin serialized at 2009.  So current run, which goes to beginning of 2019, has run duration of 10 years:
    # Burnin has run duration of 65 years:
    cb.set_param("Simulation_Duration", 65 * 365)
    cb.set_param("Serialization_Time_Steps", [54 * 365])
    end_year = 2020
    serialize_year = 2009  # Gives 1 year of buffer
    sim_start_date = "1955-01-01"

    itn_simplified = "interventions/grid_all_itn_events.csv"
    irs_simplified = "interventions/grid_all_irs_events.csv"
    mda_simplified = "interventions/grid_all_mda_events.csv"
    intervention_df_dict = preload_intervention_csvs(catch, sim_start_date,
                                                     itn_file=itn_simplified,
                                                     irs_file=irs_simplified,
                                                     mda_file=mda_simplified)
    add_intervention_combos(cb, intervention_df_dict, catch, False, False, False, False, travellers=False)
    add_all_reports(cb, catch, start=0 * 365)

    # NO IMPORTATIONS
    cb.set_param("x_Regional_Migration", 0)

    arab_burnin = 10
    funest_burnin = 9.4
    set_ento(cb, arab_burnin, funest_burnin, arab_times, arab_spline, funest_times, funest_spline)

    SetupParser.init()
    SetupParser.set("HPC", "priority", priority)
    SetupParser.set("HPC", "node_group", coreset)

    modlists = []

    new_modlist = [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed) for seed in range(num_seeds)]
    modlists.append(new_modlist)

    new_modlist = [ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat', 10)] #x) for x in range(11,20,2)]
    modlists.append(new_modlist)


    builder = ModBuilder.from_combos(*modlists)

    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(config_builder=cb,
                                exp_name=experiment_name,
                                exp_builder=builder)