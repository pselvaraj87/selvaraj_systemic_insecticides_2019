
# Add filtered reports to separate this node out
import pandas as pd
from malaria.reports.MalariaReport import add_filtered_spatial_report, add_event_counter_report, add_filtered_report

from gridded_sims.run.site import find_cells_for_this_catchment


def add_all_reports(cb, catch, start=0, duration=10000000):

    add_event_counter_report(cb, event_trigger_list=["Received_Treatment"], start=start, duration=duration)
    add_filtered_spatial_report(cb, channels=["Population", "True_Prevalence"], start=start, end=(start+duration))

    # Filtered report just for work node, and just for catchment:
    regional_EIR_node_label = 100000
    catch_node_ids = find_cells_for_this_catchment(catch)

    add_filtered_report(cb, nodes=[regional_EIR_node_label], description='Work', start=start, end=(start+duration))
    add_filtered_report(cb, nodes=catch_node_ids, description='Catchment', start=start, end=(start+duration))