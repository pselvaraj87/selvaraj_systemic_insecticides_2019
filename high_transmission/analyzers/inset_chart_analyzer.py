"""
Jaline Gerardin
Dec 2018


"""

import os
import pandas as pd
import numpy as np
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser


projectdir = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)', 'Malaria Team Folder', 'projects', 'burkina_indie')


class InsetAnalyzer(BaseAnalyzer):

    def __init__(self, expt_name, report_names=["InsetChart"], sweep_variables=None, working_dir="."):
        super(InsetAnalyzer, self).__init__(working_dir=working_dir,
                                        filenames=["output/{name}.json".format(name=name)
                                                      for name in report_names]
                                           )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.reports=report_names
        self.expt_name = expt_name

    def select_simulation_data(self, data, simulation):
        simdata = []

        for report in self.reports:

            datatemp = data["output/{name}.json".format(name=report)]

            infections_reduced = \
                np.abs(1 - datatemp['Channels']['Infected']['Data'][-1] / datatemp['Channels']['Infected']['Data'][-365])
            annual_EIR = \
                np.sum(datatemp['Channels']['Daily EIR']['Data'][-365:])

            df = pd.DataFrame(dict(zip(['Infections reduced', 'Annual EIR'], [infections_reduced, annual_EIR])), index=[0])
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
        df.to_csv(os.path.join(self.working_dir, '%s_%s.csv' % (self.expt_name, 'transmission')), index=False)


if __name__ == "__main__":

    SetupParser.default_block = 'LOCAL'
    SetupParser.init()

    out_dir = os.path.join(projectdir, 'sim_data')

    experiments = {
                   "SMC_ivermection_scenarios_uncorrelated": "01768676-af24-e911-a2be-c4346bcb1554"
                   }

    for expt_name, exp_id in experiments.items():
        am = AnalyzeManager(exp_list=exp_id,
                            analyzers=[
                                InsetAnalyzer(working_dir=out_dir,
                                                expt_name=expt_name,
                                                report_names = ['InsetChart'],
                                                sweep_variables=["Run_Number", "Coverage", "Intervention_type"])
                                       ],
                            force_analyze=False)

        print(am.experiments)
        am.analyze()

