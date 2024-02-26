from __future__ import (division, print_function)

import sys
import numpy as np
import scipy.stats as st
import scipy.linalg as la
import matplotlib
import seaborn as sns
import pandas as pd
from seaborn.categorical import boxplot

custom_params={"axes.spines.right":False,
                "axes.spines.top":False,
                "axes.grid":True,
                "axes.facecolor":'#E6E6E6',
                "grid.color":'white',
                "font.family":['serif']
                }

sns.set_theme(context='talk',style="ticks",font_scale=1.2,palette="pastel",rc=custom_params)
#sns.set_style("ticks")

matplotlib.use('Agg') #WSL compatibility 
from matplotlib import pyplot as plt
from matplotlib import ticker as mticker

def boxPlot(output, filename):
    # setup plot
    fig, ax = plt.subplots(figsize=(12,8))

    # Count number of run
    number_runs = np.count_nonzero(output[:,2] == 1)
    
    # Noramlize runtimes
    sum = np.sum(output[:,1])
    normalized_time = np.divide(output[:,1],sum)
    output[:,1] = normalized_time * number_runs * 100

    data = pd.DataFrame(output, columns=['Sequence Length','Time','Process Rank'])

    sns.boxplot(data=data, x='Process Rank',y='Time',hue='Sequence Length',palette=['#c27ba0'])

    # legend = ax.get_legend()
    # handles = legend.legendHandles
    # legend.remove()
    # ax.legend(handles,[sys.argv[3],sys.argv[2]])

    plt.title('Process Occupation in [%]', loc='left',c='grey',fontsize='large')
    plt.suptitle('Computation Time per Process',x=0.125, y=0.98, ha='left',fontsize='xx-large',fontweight=500)
    plt.xlabel('Process Rank')
    plt.ticklabel_format(axis='y',style='sci', scilimits=(-2,2))
    ax.set(ylabel=None)

    # Format yaxis as percentages
    fmt = '%.2f%%'
    yticks = matplotlib.ticker.FormatStrFormatter(fmt)
    ax.yaxis.set_major_formatter(yticks)

    #ax.xaxis.set_major_locator(plt.MaxNLocator(18))

    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2))
    ax.xaxis.set_major_formatter('{x:.0f}')

    # Set y-axis limit to 100%
    #ax.set_ybound(0,100)

    fig.align_ylabels()
    fig.savefig('plots/'+str(filename),dpi=500)

def violinPlot(output, filename):
    # setup plot
    fig, ax = plt.subplots(figsize=(12,8))

    # Count number of run
    number_runs = np.count_nonzero(output[:,2] == 1)
    
    # Noramlize runtimes
    sum = np.sum(output[:,1])
    normalized_time = np.divide(output[:,1],sum)
    output[:,1] = normalized_time * number_runs * 100

    data = pd.DataFrame(output, columns=['Sequence Length','Time','Process Rank'])

    sns.violinplot(data=data, x='Process Rank',y='Time',hue='Sequence Length',palette="tab10")

    # legend = ax.get_legend()
    # handles = legend.legendHandles
    # legend.remove()
    # ax.legend(handles,[sys.argv[3],sys.argv[2]])

    plt.title('Process Occupation in [%]', loc='left',c='grey',fontsize='large')
    plt.suptitle('Computation Time per Process',x=0.125, y=0.96, ha='left',fontsize='xx-large',fontweight=500)
    plt.xlabel('Process Rank')
    plt.ticklabel_format(axis='y',style='sci', scilimits=(-2,2))
    ax.set(ylabel=None)

    # Format yaxis as percentages
    fmt = '%.0f%%'
    yticks = matplotlib.ticker.FormatStrFormatter(fmt)
    ax.yaxis.set_major_formatter(yticks)

    # Set y-axis limit to 100%
    #ax.set_ybound(0,100)

    fig.align_ylabels()
    fig.savefig('plots/'+str(filename),dpi=500)

def main():
    # retrieve data from output
    output_3 = np.loadtxt(sys.argv[1],delimiter='\t')
    #output_4 = np.loadtxt(sys.argv[2],delimiter='\t')

    boxPlot(output_3, 'boxPlot_procs_occupation.png')
    #violinPlot(output_3, 'violinPlot_procs_occupation.png')
    #boxPlot(output_4, 'boxPlot_threads_occupation.png')


    
if __name__ == '__main__':
    main()