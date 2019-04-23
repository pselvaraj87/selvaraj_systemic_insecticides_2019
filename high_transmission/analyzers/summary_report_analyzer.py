"""
Prashanth Selvaraj
Feb 2019


"""

import os
import pandas as pd
import numpy as np
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser


class SummaryAnalyzer(BaseAnalyzer):

    def __init__(self, expt_name, report_names=["Daily"], sweep_variables=None, working_dir="."):
        super(SummaryAnalyzer, self).__init__(working_dir=working_dir,
                                        filenames=["output/MalariaSummaryReport_{name}.json".format(name=name)
                                                      for name in report_names]
                                           )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.reports=report_names
        self.expt_name = expt_name
        self.age_bins = [5, 15, 100]

    def select_simulation_data(self, data, simulation):
        simdata = []

        for report in self.reports:

            datatemp = data["output/MalariaSummaryReport_{name}.json".format(name=report)]

            incidence = datatemp['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin']
            population = datatemp['DataByTimeAndAgeBins']['Average Population by Age Bin']
            clinical_cases = (np.sum(np.array(incidence)*np.array(population), axis=0)/365).tolist()

            df = pd.DataFrame(dict(zip([5, 15, 100], clinical_cases)), index=[0])
            simdata.append(df)
        simdata = pd.concat(simdata)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
            else:
                simdata[sweep_var] = 0
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        df = pd.concat(selected).reset_index(drop=True)
        df.to_csv(os.path.join(self.working_dir, '%s_%s.csv' % (self.expt_name, 'clinical_cases')), index=False)


if __name__ == "__main__":

    SetupParser.default_block = 'HPC'
    SetupParser.init()

    out_dir = os.path.join(projectdir, 'sim_data')

    experiments = {
                   "SMC_ivermection_scenarios_uncorrelated": "99b7796b-4d2b-e911-a2c5-c4346bcb7273"
                   }

    for expt_name, exp_id in experiments.items():
        am = AnalyzeManager(exp_list=exp_id,
                            analyzers=[
                                SummaryAnalyzer(working_dir=out_dir,
                                                expt_name=expt_name,
                                                report_names = ['Daily_Report'],
                                                sweep_variables=["Run_Number", "Coverage", "Intervention_type"]),

                                       ],
                            force_analyze=False)

        print(am.experiments)
        am.analyze()

