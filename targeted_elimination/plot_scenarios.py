import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import os
from plotting.colors import load_color_palette

mpl.rcParams['pdf.fonttype'] = 42

userpath = 'E:/'

wdir = os.path.join(userpath, 'Dropbox (IDM)', 'Malaria Team Folder/projects/cambodia_forest_scenarios/ivermectin')
plotdir = os.path.join(wdir, 'plots')
datadir = os.path.join(wdir, 'data')


def plot_by_intervention(adf, plotname) :

    num_years = len(adf['year'].unique())
    num_targets = len(adf['target'].unique())

    ann_prev_channel = 'annual mean Infected'
    day_prev_channel = 'year_end_Infected'

    sweep_variables = ["x_Temporary_Larval_Habitat",
                       'dirus_Anthropophily',
                       ]

    sns.set_style('whitegrid', {'axes.linewidth' : 0.5})
    palette = load_color_palette('wes')
    for (x_temp, anth), df in adf.groupby(sweep_variables) :
        fig1 = plt.figure('%s Anth%.2f LH%.1f prevalence' % (plotname, anth, x_temp), figsize=(10,6))
        fig1.subplots_adjust(left=0.07, right=0.95)
        axes1 = [fig1.add_subplot(num_targets, num_years, i + 1) for i in range(num_years * num_targets)]
        fig2 = plt.figure('%s Anth%.2f LH%.1f fraction eliminated' % (plotname, anth, x_temp), figsize=(10,6))
        fig2.subplots_adjust(left=0.07, right=0.95)
        axes2 = [fig2.add_subplot(num_targets, num_years, i + 1) for i in range(num_years * num_targets)]
        for i, ((target, year), gdf) in enumerate(df.groupby(['target', 'year'])) :
            ax1 = axes1[i]
            ax2 = axes2[i]
            for j, ((intervention, ivm_duration), ddf) in enumerate(gdf.groupby(['intervention','IVM duration'])) :
                if 'drug' in intervention :
                    if 'ivermectin' not in intervention :
                        label = 'DP'
                    else :
                        label = 'DP + End%d' % ivm_duration
                elif 'ivermectin' in intervention :
                    label = 'End%d' % ivm_duration
                else :
                    label = ''

                sdf = ddf.groupby('coverage')[ann_prev_channel].agg(np.mean).reset_index()
                ax1.plot(sdf['coverage'], sdf[ann_prev_channel], color=palette[j], label=label)
                sdf = ddf.groupby('coverage')[day_prev_channel].agg(count_zero).reset_index()
                ax2.plot(sdf['coverage'], sdf[day_prev_channel], color=palette[j], label=label)

            for ax in [ax1, ax2] :
                ax.set_xlim(0,1)
                if i%num_years > 0 :
                    ax.set_yticklabels([])
                else :
                    ax.set_ylabel(target)
                if i < num_years :
                    ax.set_title('year %d' % year)
                if i < len(axes1)-num_years :
                    ax.set_xticklabels([])
                else :
                    ax.set_xlabel('coverage')
            ax1.set_ylim(0, 0.2)
            ax2.set_ylim(0, 1)
        axes1[-1].legend()
        axes2[-1].legend()
        fig2.savefig(os.path.join(plotdir, '%s Anth%.2f LH%.1f fraction eliminated.png' % (plotname, anth, x_temp)))
        fig2.savefig(os.path.join(plotdir, '%s Anth%.2f LH%.1f fraction eliminated.pdf' % (plotname, anth, x_temp)),
                     format='PDF')


def plot_by_anth(adf, plotname) :

    num_years = len(adf['year'].unique())
    num_targets = len(adf['target'].unique())
    day_prev_channel = 'year_end_Infected'

    palette = load_color_palette('wes')

    sns.set_style('whitegrid', {'axes.linewidth' : 0.5})
    fig = plt.figure(plotname, figsize=(10, 6))
    fig.subplots_adjust(left=0.07, right=0.95)
    axes = [fig.add_subplot(num_targets, num_years, i + 1) for i in range(num_years * num_targets)]
    for i, ((target, year), gdf) in enumerate(adf.groupby(['target', 'year'])):
        ax = axes[i]

        for j, ((intervention, ivm_duration), ddf) in enumerate(gdf.groupby(['intervention', 'IVM duration'])):
            if 'drug' in intervention:
                if 'ivermectin' not in intervention:
                    label = 'DP'
                else:
                    label = 'DP + End%d' % ivm_duration
            elif 'ivermectin' in intervention:
                label = 'End%d' % ivm_duration
            else:
                label = ''

            sdf = ddf.groupby('dirus_Anthropophily')[day_prev_channel].agg(count_zero).reset_index()
            ax.plot(sdf['dirus_Anthropophily'], sdf[day_prev_channel], color=palette[j], label=label)

        if i % num_years > 0:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel(target)
        if i < num_years:
            ax.set_title('year %d' % year)
        if i < len(axes) - num_years:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('dirus anthropophily')
        ax.set_ylim(0, 1)
        ax.set_xlim(0, 1)
    axes[-1].legend()

    fig.savefig(os.path.join(plotdir, '%s fraction eliminated.png' % plotname))
    fig.savefig(os.path.join(plotdir, '%s fraction eliminated.pdf' % plotname), format='PDF')


def plot_by_event_count(adf, plotname) :

    num_years = len(adf['year'].unique())
    num_targets = len(adf['target'].unique())
    day_prev_channel = 'year_end_Infected'

    palette = load_color_palette('wes')

    sns.set_style('whitegrid', {'axes.linewidth' : 0.5})
    fig = plt.figure(plotname, figsize=(10, 6))
    fig.subplots_adjust(left=0.07, right=0.95)
    axes = [fig.add_subplot(num_targets, num_years, i + 1) for i in range(num_years * num_targets)]
    for i, ((target, year), gdf) in enumerate(adf.groupby(['target', 'year'])):
        ax = axes[i]

        for j, ((intervention, ivm_duration), ddf) in enumerate(gdf.groupby(['intervention', 'IVM duration'])):
            if 'drug' in intervention:
                if 'ivermectin' not in intervention:
                    label = 'DP'
                else:
                    label = 'DP + End%d' % ivm_duration
            elif 'ivermectin' in intervention:
                label = 'End%d' % ivm_duration
            else:
                label = ''

            cdf = adf[(adf['target'] == target) & (adf['intervention'] == intervention) &
                      (adf['IVM duration'] == ivm_duration) & (adf['year'] <= year)]
            num_seeds = np.max(ddf['Run_Number'])

            zeroes = ddf.groupby('coverage')[day_prev_channel].agg(count_zero).reset_index()
            events = cdf.groupby('coverage')['Event_Count'].agg(np.sum).reset_index()
            ax.plot(events['Event_Count']/num_seeds, zeroes[day_prev_channel], color=palette[j], label=label)

        if i % num_years > 0:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel(target)
        if i < num_years:
            ax.set_title('year %d' % year)
        if i < len(axes) - num_years:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('doses distributed')
        ax.set_ylim(0, 1)
        ax.set_xlim(0, 1700)
    axes[-1].legend()
    fig.savefig(os.path.join(plotdir, '%s fraction eliminated v doses.png' % plotname))
    fig.savefig(os.path.join(plotdir, '%s fraction eliminated v doses.pdf' % plotname), format='PDF')


def plot_by_ivm_duration(df, plotname) :

    num_years = len(df['year'].unique())
    num_targets = len(df['target'].unique())
    min_ivm_duration = np.min(df[df['IVM duration'] > 0]['IVM duration'])
    max_ivm_duration = np.max(df['IVM duration'])

    day_prev_channel = 'year_end_Infected'

    sns.set_style('whitegrid', {'axes.linewidth' : 0.5})
    palette = load_color_palette('wes')
    fig = plt.figure(plotname, figsize=(10,6))
    fig.subplots_adjust(left=0.07, right=0.95)
    axes = [fig.add_subplot(num_targets, num_years, i + 1) for i in range(num_years * num_targets)]
    for i, ((target, year), gdf) in enumerate(df.groupby(['target', 'year'])) :
        ax = axes[i]
        for j, (intervention, ddf) in enumerate(gdf.groupby(['intervention'])) :

            sdf = ddf.groupby('IVM duration')[day_prev_channel].agg(count_zero).reset_index()
            if 'ivermectin' not in intervention :
                ax.plot([min_ivm_duration, max_ivm_duration], [sdf[day_prev_channel], sdf[day_prev_channel]],
                        color=palette[j], label=intervention)
            else :
                ax.plot(sdf['IVM duration'], sdf[day_prev_channel], color=palette[j], label=intervention)

        if i%num_years > 0 :
            ax.set_yticklabels([])
        else :
            ax.set_ylabel(target)
        if i < num_years :
            ax.set_title('year %d' % year)
        if i < len(axes)-num_years :
            ax.set_xticklabels([])
        else :
            ax.set_xlabel('IVM duration')
        ax.set_ylim(0, 1)
    axes[-1].legend()
    fig.savefig(os.path.join(plotdir, '%s fraction eliminated v IVM duration.png' % plotname))
    fig.savefig(os.path.join(plotdir, '%s fraction eliminated v IVM duration.pdf' % plotname), format='PDF')


def plot_scenarios_by_coverage(adf, anth, LH) :

    scenarios = ['MDA', 'to_forest']
    df = adf[(adf['dirus_Anthropophily'] == anth) & (adf['x_Temporary_Larval_Habitat'] == LH)]
    for scenario_name in scenarios :
        plot_by_intervention(df[df['intervention'].str.contains(scenario_name)], scenario_name)


def plot_scenarios_by_anth(adf, LH, coverage) :

    scenarios = ['to_forest', 'MDA', 'HS']

    df = adf[adf['x_Temporary_Larval_Habitat'] == LH]
    df = df[df['coverage'] == coverage]

    for scenario_name in scenarios :
        plot_by_anth(df[df['intervention'].str.contains(scenario_name)], '%s cov%.1f' % (scenario_name, coverage))


def plot_scenarios_by_event_count(adf, anth, LH) :

    scenarios = ['MDA', 'to_forest']
    df = adf[(adf['dirus_Anthropophily'] == anth) & (adf['x_Temporary_Larval_Habitat'] == LH)]
    for scenario_name in scenarios :
        plot_by_event_count(df[df['intervention'].str.contains(scenario_name)],
                            '%s Anth%.2f LH%.1f' % (scenario_name, anth, LH))


def plot_scenarios_by_ivm_duration(adf, anth, LH) :

    scenarios = ['MDA', 'to_forest']
    df = adf[(adf['dirus_Anthropophily'] == anth) & (adf['x_Temporary_Larval_Habitat'] == LH)]
    for scenario_name in scenarios :
        plot_by_ivm_duration(df[df['intervention'].str.contains(scenario_name)],
                             scenario_name)


def count_zero(x) :
    return (len(x) - np.count_nonzero(x))/len(x)


if __name__ == '__main__' :

    adf = pd.read_csv(os.path.join(datadir, 'End_all.csv'))
    adf = adf[adf['year'] > 0]

    plot_scenarios_by_coverage(adf, anth=0.5, LH=2.4)
    plot_scenarios_by_anth(adf, LH=2.4, coverage=0.6)
    plot_scenarios_by_event_count(adf, anth=0.5, LH=2.4)
    plot_scenarios_by_ivm_duration(adf, anth=0.5, LH=2.4)
