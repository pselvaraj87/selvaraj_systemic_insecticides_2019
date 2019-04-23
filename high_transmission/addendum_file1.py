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

    fig = plt.figure(tag, figsize=(20, 5))
    fig.subplots_adjust(left=0.05, right=0.98)
    palette = load_color_palette()

    for i, ia in enumerate(['Annual EIR', '5', '15', '100']):
        ax = fig.add_subplot(1, 4, i+1)
        plt.grid()
        for iarm, (comp, compdf) in enumerate(df.groupby(to_compare)):
            ax.plot(compdf['Coverage'], compdf[ia + '_mean'], '-', color=palette[iarm], label=comp)
            ax.fill_between(compdf['Coverage'], compdf[ia + '_mean']-compdf[ia + '_std'], compdf[ia + '_mean'] + compdf[ia + '_std'], color=palette[iarm], linewidth=0, alpha=0.3)
        if ia == 'Annual EIR':
            ax.set_title('%s' % ia)
        else:
            ax.set_title('clinical cases <%sy' % ia)
        # ax.set_ylim(ymin, max(df[ia + '_mean']) + 1.1*max(df[ia + '_std']))
        ax.set_ylim(0.0, 1.0)
        ax.set_xlim(0.0, 1.0)
        if i == 3:
            ax.legend()
    fig.text(0.5, -0.01, 'Coverage', ha='center', va='center')
    fig.text(0.0, 0.5, '% reduction', ha='center', va='center', rotation='vertical')

    if save:
        save_double(fig, savename=savename)

    return fig


def save_double(fig, savename) :

    if savename :
        fig.savefig(os.path.join(plotdir, '%s.png' % savename), bbox_inches="tight")
        fig.savefig(os.path.join(plotdir, '%s.pdf' % savename), format='PDF', bbox_inches="tight")


if __name__ == '__main__' :

    sim_expt_name = "eSMC10_IVM"
    smc = 'DP10'
    # sim_expt_name = "IVM_final"
    # smc = 'DP15'
    addendum_expt_name = "eSMC_IVM_everyone"
    addendum_expt_name2 = "SMC_10"
    ref_expt_name = "IVM_final"
    sim_expt_tags = ["%s\+ExpandedSMC\+" % smc,
                     "%s\+ExpandedSMCexcludingwomen\+" % smc,
                     "%s\+Everyone\+" % smc
                     ]
    smc_types = ['%s' % smc]
    ivm_durations = [7, 14, 30, 90]

    summary_data_fname = os.path.join(datadir, "%s_clinical_cases.csv" % sim_expt_name)
    summary_df = pd.read_csv(summary_data_fname)
    if addendum_expt_name:
        add_summary_data_fname = os.path.join(datadir, "%s_clinical_cases.csv" % addendum_expt_name)
        add_summary_df = pd.read_csv(add_summary_data_fname)
        summary_df = pd.concat([summary_df, add_summary_df], ignore_index=True)
    if addendum_expt_name2:
        add_summary_data_fname = os.path.join(datadir, "%s_clinical_cases.csv" % addendum_expt_name2)
        add_summary_df = pd.read_csv(add_summary_data_fname)
        if 'Unnamed: 0' in add_summary_df.columns:
            del add_summary_df['Unnamed: 0']
        summary_df = pd.concat([summary_df, add_summary_df], ignore_index=True)

    inset_data_fname = os.path.join(datadir, "%s_transmission.csv" % sim_expt_name)
    inset_df = pd.read_csv(inset_data_fname)
    if addendum_expt_name:
        add_inset_data_fname = os.path.join(datadir, "%s_transmission.csv" % addendum_expt_name)
        add_inset_df = pd.read_csv(add_inset_data_fname)
        inset_df = pd.concat([inset_df, add_inset_df], ignore_index=True)
    if addendum_expt_name2:
        add_inset_data_fname = os.path.join(datadir, "%s_transmission.csv" % addendum_expt_name2)
        add_inset_df = pd.read_csv(add_inset_data_fname)
        if 'Unnamed: 0' in add_inset_df.columns:
            del add_inset_df['Unnamed: 0']
        inset_df = pd.concat([inset_df, add_inset_df], ignore_index=True)

    df = pd.merge(summary_df, inset_df, on=['Coverage', 'Run_Number', 'Intervention_type'])
    df = df.sort_values(by=['Coverage', 'Run_Number', 'Intervention_type'])

    summary_data_fname = os.path.join(datadir, "%s_clinical_cases.csv" % ref_expt_name)
    summary_df = pd.read_csv(summary_data_fname)
    inset_data_fname = os.path.join(datadir, "%s_transmission.csv" % ref_expt_name)
    inset_df = pd.read_csv(inset_data_fname)
    ref_df = pd.merge(summary_df, inset_df, on=['Coverage', 'Run_Number', 'Intervention_type'])
    ref_df = ref_df.sort_values(by=['Coverage', 'Run_Number', 'Intervention_type'])
    ref_df = ref_df[ref_df['Intervention_type'] == 'SMC5']

    columns = ['5', '15', '100', 'Annual EIR']
    final_df = data_wrangling(df, ref_df, columns_to_keep=columns)

    rows = sim_expt_tags   # Each row is comparing labels across all simulations
    labels = ivm_durations
    columns_to_groupby = ['Coverage', 'Intervention_type']
    final_df_mean = final_df.groupby(columns_to_groupby)[columns].apply(np.mean).reset_index()
    final_df_std = final_df.groupby(columns_to_groupby)[columns].apply(np.std).reset_index()
    df = pd.merge(final_df_mean, final_df_std, on=['Coverage', 'Intervention_type'], suffixes=('_mean', '_std'))

    if rows == sim_expt_tags:
        for tag in sim_expt_tags:
            dftemp = df[df['Intervention_type'].str.contains(tag)]
            if 'DP10' in tag:
                dfesmc = df[df['Intervention_type'] == 'SMC10']
                dftemp = pd.concat([dftemp, dfesmc])
            plot_data(dftemp, tag, to_compare='Intervention_type', ymin=0, ymax=100, savename='%s_EIR_clinical' %tag.replace('\\', ''))

    exit()

    plt.show()

