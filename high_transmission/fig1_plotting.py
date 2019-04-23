"""
Prashanth Selvaraj
Feb 2019

"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess
mpl.rcParams['pdf.fonttype'] = 42

from plotting.colors import load_color_palette

projectdir = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)', 'Malaria Team Folder', 'projects', 'endectocides', 'prashanth_sahel')
datadir = os.path.join(projectdir, 'sim_data')
plotdir = os.path.join(projectdir, 'sim_plots')
if not os.path.exists(plotdir):
    os.mkdir(plotdir)


def data_wrangling(sim_df, ref_df, columns_to_keep=None):
    final_df = pd.DataFrame()
    sim_tags = list(sim_df['Intervention_type'].unique())
    for tag in sim_tags:
        df = sim_df[sim_df['Intervention_type'] == tag]
        df[columns_to_keep] = (-df[columns_to_keep] + ref_df[columns_to_keep].values)/ref_df[columns_to_keep].values
        final_df = pd.concat([final_df, df])

    return final_df


def plot_data(df, tag, to_compare=None, ymin=None, ymax=None, savename='None', save=1):

    fig = plt.figure(tag, figsize=(5, 5))
    fig.subplots_adjust(left=0.05, right=0.98)
    palette = load_color_palette()

    for i, ia in enumerate(['Total']):
        ax = fig.add_subplot(1, 1, i+1)
        plt.grid()
        for iarm, (comp, compdf) in enumerate(df.groupby(to_compare)):
            if comp == 'SMC15':
                ax.plot(compdf['Coverage'], compdf[ia + '_mean'], '-', color=palette[iarm], label=comp)
            elif comp == 'SMC10':
                ax.plot(compdf['Coverage'], compdf[ia + '_mean'], '-', color=palette[iarm], label=comp)
            else:
                ax.plot(compdf['Coverage'], compdf[ia + '_mean'], '-', color=palette[iarm], label=comp.split('+')[2])
            ax.fill_between(compdf['Coverage'], compdf[ia + '_mean']-compdf[ia + '_std'], compdf[ia + '_mean'] + compdf[ia + '_std'], color=palette[iarm], linewidth=0, alpha=0.3)
        if ia == 'Annual EIR':
            ax.set_title('%s' % ia)
        else:
            ax.set_title('clinical cases')
        ax.set_ylim(ymin, ymax)
        ax.set_xlim(0.0, 1.0)
        if i == 0:
            ax.legend(bbox_to_anchor=(1.05, 0.5),)
    fig.suptitle('%s' % savename)
    fig.text(0.5, -0.01, 'Coverage', ha='center', va='center')
    fig.text(0.0, 0.5, '% reduction', ha='center', va='center', rotation='vertical')

    if save:
        save_double(fig, savename=savename)

    return fig


def save_double(fig, savename) :

    if savename :
        fig.savefig(os.path.join(plotdir, '%s_total.png' % savename), bbox_inches="tight")
        fig.savefig(os.path.join(plotdir, '%s_total.pdf' % savename), format='PDF', bbox_inches="tight")


if __name__ == '__main__' :

    sim_expt_name = "All_expts"
    sim_expt_tags = [
        'DP5\+Everyone', 'DP5\+Excluding\<5', 'DP5\+Excludingwomenand\<5',
        'DP10\+Everyone\+',
        'DP10\+ExpandedSMC\+', 'DP10\+ExpandedSMCexcludingwomen',
        'DP15\+Everyone\+', 'DP15\+ExpandedSMC\+', 'DP15\+ExpandedSMCexcludingwomen'
    ]
    save_names = [
        'SMC5_Everyone', 'SMC5_nochildren', 'SMC5_nowomenandchildren',
        'eSMC10_Everyone',
        'eSMC10_nochildren', 'eSMC10_nowomenandchildren',
        'eSMC15_Everyone', 'eSMC15_nochildren', 'eSMC15_nowomenandchildren'
                  ]

    smc_types = ['DP5', 'DP15']
    ivm_durations = [7, 14, 30, 90]

    summary_data_fname = os.path.join(datadir, "%s_clinical_cases.csv" % sim_expt_name)
    df = pd.read_csv(summary_data_fname)

    ref_df = df[df['Intervention_type'] == 'SMC5']

    columns = ['Total']
    final_df = data_wrangling(df, ref_df, columns_to_keep=columns)

    rows = sim_expt_tags   # Each row is comparing labels across all simulations
    labels = ivm_durations
    columns_to_groupby = ['Coverage', 'Intervention_type']
    final_df_mean = final_df.groupby(columns_to_groupby)[columns].apply(np.mean).reset_index()
    final_df_std = final_df.groupby(columns_to_groupby)[columns].apply(np.std).reset_index()
    df = pd.merge(final_df_mean, final_df_std, on=['Coverage', 'Intervention_type'], suffixes=('_mean', '_std'))

    if rows == sim_expt_tags:
        for i, tag in enumerate(sim_expt_tags):
            dftemp = df[df['Intervention_type'].str.contains(tag)]
            if 'DP10' in tag:
                dfesmc = df[df['Intervention_type'] == 'SMC10']
                dftemp = pd.concat([dftemp, dfesmc])
            elif 'DP15' in tag:
                dfesmc = df[df['Intervention_type'] == 'SMC15']
                dftemp = pd.concat([dftemp, dfesmc])
            plot_data(dftemp, tag, to_compare='Intervention_type', ymin=0, ymax=1.0, savename=save_names[i])
    exit()

    plt.show()

