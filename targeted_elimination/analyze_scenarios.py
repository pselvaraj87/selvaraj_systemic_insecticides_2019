import pandas as pd
import numpy as np
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser

import os

if not SetupParser.initialized:
    SetupParser.init('HPC')

userpath = 'E:/'

wdir = os.path.join(userpath, 'Dropbox (IDM)', 'Malaria Team Folder/projects/cambodia_forest_scenarios/ivermectin')
datadir = os.path.join(wdir, 'data')


class IvermectinAnalyzer(BaseAnalyzer) :

    def __init__(self, expt_name, sweep_variables, inset_channel='Infected', working_dir='.'):
        super(IvermectinAnalyzer, self).__init__(working_dir=working_dir,
                                                 filenames=["output/InsetChart.json",
                                                            'output/ReportEventCounter.json']
                                           )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.channel = inset_channel
        self.event_channel = 'Event_Count'
        self.expt_name = expt_name

    def select_simulation_data(self, data, simulation):

        if 'ivermectin' in simulation.tags['intervention'] :
            event_channel = 'Received_Ivermectin'
        elif 'drug' in simulation.tags['intervention'] :
            event_channel = 'Received_Campaign_Drugs'
        else :
            event_channel = 'Received_Treatment'

        channeldata = data[self.filenames[0]]['Channels'][self.channel]['Data']
        eventdata = data[self.filenames[1]]['Channels'][event_channel]['Data']

        simdata = pd.DataFrame( { self.channel : channeldata,
                                  self.event_channel: eventdata,
                                  'time' : list(range(len(channeldata)))})

        simdata['year'] = simdata['time'].apply(lambda x : int(x/365))

        # annual average
        adf = simdata.groupby('year')[self.channel].agg(np.mean).reset_index()
        adf.rename(columns={self.channel : 'annual mean %s' % self.channel}, inplace=True)

        # end of the year
        df = simdata[(simdata['time'] % 365 == 364)]
        df = df[[self.channel, 'year']]
        df.rename(columns={self.channel : 'year_end_%s' % self.channel}, inplace=True)
        adf = pd.merge(left=adf, right=df, on='year')

        # annual events (different by intervention type)
        df = simdata.groupby('year')[self.event_channel].agg(np.sum).reset_index()
        adf = pd.merge(left=adf, right=df, on='year')

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                adf[sweep_var] = simulation.tags[sweep_var]
            else:
                adf[sweep_var] = 0

        return adf

    def finalize(self, all_data):

        selected = [data for sim,data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)

        adf.sort_values(by=self.sweep_variables + ['year'], inplace=True)
        adf.to_csv(os.path.join(self.working_dir, '%s.csv' % self.expt_name), index=False)


if __name__ == '__main__' :

    out_dir = os.path.join(wdir, 'data')

    experiments = {
        "End60" : "229a4440-8d1f-e911-a2be-c4346bcb1554"
                   }

    for expt_name, exp_id in experiments.items():
        am = AnalyzeManager(exp_list=exp_id,
                            analyzers=[IvermectinAnalyzer(working_dir=out_dir,
                                                          expt_name=expt_name,
                                                          sweep_variables=["Run_Number",
                                                                           "x_Temporary_Larval_Habitat",
                                                                           'dirus_Anthropophily',
                                                                           'target',
                                                                           'coverage',
                                                                           "intervention"
                                                                           ])],
                            force_analyze=True)

        print(am.experiments)
        am.analyze()
