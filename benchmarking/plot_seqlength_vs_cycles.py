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
#context can be: paper, notebook, talk, poster
sns.set_theme(context='talk',style="ticks",font_scale=1.2,palette="pastel",rc=custom_params)
#sns.set_style("ticks")

matplotlib.use('Agg') #WSL compatibility 
from matplotlib import pyplot as plt
from matplotlib import ticker as mticker

# Helper class for formatting
class MathTextSciFormatter(mticker.Formatter):
        def __init__(self, fmt="%1.2e"):
            self.fmt = fmt
        def __call__(self, x, pos=None):
            s = self.fmt % x
            decimal_point = '.'
            positive_sign = '+'
            tup = s.split('e')
            significand = tup[0].rstrip(decimal_point)
            sign = tup[1][0].replace(positive_sign, '')
            exponent = tup[1][1:].lstrip('0')
            if exponent:
                exponent = '10^{%s%s}' % (sign, exponent)
            if significand and exponent:
                s =  r'%s{\times}%s' % (significand, exponent)
            else:
                s =  r'%s%s' % (significand, exponent)
            return "${}$".format(s)

# absolute time
def plotAbs_Seq(output_00, output_01, output_02, output_03, output_04,output_05, output_06,y_ax,filename):
    # setup plot
    fig, ax = plt.subplots(figsize=(12,8))

    # create suitable np array
    output_00f = np.insert(output_00, 2, np.zeros(np.size(output_00,0)), axis=1)
    output_01f = np.insert(output_01, 2, np.ones(np.size(output_01,0)), axis=1)
    output_02f = np.insert(output_02, 2, np.full(np.size(output_02,0),2), axis=1)
    output_03f = np.insert(output_03, 2, np.full(np.size(output_03,0),3), axis=1)
    output_04f = np.insert(output_04, 2, np.full(np.size(output_04,0),4), axis=1)
    output_05f = np.insert(output_05, 2, np.full(np.size(output_05,0),5), axis=1)
    output_06f = np.insert(output_06, 2, np.full(np.size(output_06,0),6), axis=1)

    #output = np.concatenate((output_00f,output_01f,output_02f,output_03f,output_05f,output_04f,output_06f,output_10f,output_12f,output_13f,output_14f))
    output = np.concatenate((output_00f,output_01f,output_02f,output_03f,output_04f,output_05f,output_06f))
    
    # convert to pd.DataFrame 
    data = pd.DataFrame(output, columns=['N','Count','Algorithm'])
    g = sns.lineplot(data=data, x='N', y='Count', hue='Algorithm', palette=['#c90076','#6a329f','#2986cc','#16537e','#6aa84f','#38761d','#f44336'], linewidth=.5, markers=True, dashes=False, markersize=1, markeredgewidth=0.2, ci='sd', style='Algorithm')
    
    # allow recreation of legend
    legend = ax.get_legend()
    handles = legend.legendHandles
    legend.remove()
    ax.legend(handles, ['Cell Matrix', 'Integer Matrices', 'Diagonal Index Matrix','Diagonal Integer Vectors','Diagonal Integer Vectors SIMD','Diagonal Integer Vectors Nested','Diagonal Integer Vectors Pragma SIMD'],loc='lower right',prop={'size': 15}) # 'Sequential Diagonal' 'Sequential Diagonal Simd', 'Sequential Diagonal Simd Pragma'
    #ax.legend(handles, ['Sequential Baseline', 'Sequential Optimized' ,'Sequential Diagonal', 'Sequential Diagonal Optimized','Sequential Diagonal Simd','Sequential Diagonal Simd Pragma','Sequential Diagonal Optimized Nested','Parallel Baseline','Parallel Shared','Parallel OpenMp','Parallel OpenMpOpt']) # 'Sequential Diagonal' 'Sequential Diagonal Simd', 'Sequential Diagonal Simd Pragma'
    #ax.legend(handles, [sys.argv[1], sys.argv[2], sys.argv[3]])
    

    # Get number max seq
    max_seq = int(np.amax(output_00[:,0]))

    # Get number of runs
    number_runs = np.count_nonzero(output_00[:,0] == max_seq)

    # Get number of total runs
    tot_runs = output_00.shape[0]

    runtime_00 = np.mean(output_00[tot_runs-number_runs:tot_runs,1])
    runtime_01 = np.mean(output_01[tot_runs-number_runs:tot_runs,1])
    runtime_02 = np.mean(output_02[tot_runs-number_runs:tot_runs,1])
    runtime_03 = np.mean(output_03[tot_runs-number_runs:tot_runs,1])
    runtime_04 = np.mean(output_04[tot_runs-number_runs:tot_runs,1])
    runtime_05 = np.mean(output_05[tot_runs-number_runs:tot_runs,1])
    runtime_06 = np.mean(output_06[tot_runs-number_runs:tot_runs,1])

    ax.annotate(round(runtime_00,2), xy=(max_seq, runtime_00), xycoords='data',
            xytext=(50, 25), textcoords='offset points', ha='center',size=12, color='#c90076',
            bbox=dict(boxstyle="round4", fc="0.9",color='#c90076'),arrowprops=dict(arrowstyle="->",color='#c90076'))
    ax.annotate(round(runtime_01,2), xy=(max_seq, runtime_01), xycoords='data',
            xytext=(50, -25), textcoords='offset points', ha='center',size=12, color='#6a329f',
            bbox=dict(boxstyle="round4", fc="0.9",color='#6a329f'),arrowprops=dict(arrowstyle="->",color='#6a329f'))
    ax.annotate(round(runtime_02,2), xy=(max_seq, runtime_02), xycoords='data',
            xytext=(50, 10), textcoords='offset points', ha='center',size=12, color='#2986cc',
            bbox=dict(boxstyle="round4", fc="0.9",color='#2986cc'),arrowprops=dict(arrowstyle="->",color='#2986cc'))
    ax.annotate(round(runtime_03,2), xy=(max_seq, runtime_03), xycoords='data',
            xytext=(50, 15), textcoords='offset points', ha='center',size=12, color='#16537e',
            bbox=dict(boxstyle="round4", fc="0.9",color='#16537e'),arrowprops=dict(arrowstyle="->",color='#16537e'))
    ax.annotate(round(runtime_04,2), xy=(max_seq, runtime_04), xycoords='data',
            xytext=(50, -5), textcoords='offset points', ha='center',size=12, color='#6aa84f',
            bbox=dict(boxstyle="round4", fc="0.9",color='#6aa84f'),arrowprops=dict(arrowstyle="->",color='#6aa84f'))
    ax.annotate(round(runtime_05,2), xy=(max_seq, runtime_05), xycoords='data',
            xytext=(50, -15), textcoords='offset points', ha='center',size=12, color='#38761d',
            bbox=dict(boxstyle="round4", fc="0.9",color='#38761d'),arrowprops=dict(arrowstyle="->",color='#38761d'))
    ax.annotate(round(runtime_06,2), xy=(max_seq, runtime_06), xycoords='data',
            xytext=(50, 0), textcoords='offset points', ha='center',size=12, color='#f44336',
            bbox=dict(boxstyle="round4", fc="0.9",color='#f44336'),arrowprops=dict(arrowstyle="->",color='#f44336'))


    #set log-yaxis
    if(y_ax == 1):
        ax.set_yscale('log')
    
    # set title and format
    plt.title('[seconds]', loc='left',c='grey', fontsize='large')
    plt.suptitle('Sequential Absolute Runtimes',x=0.125, y=0.980, ha='left',fontsize='xx-large', fontweight=10)
    plt.xlabel('Sequence Length')
    ax.set(ylabel=None)

    xlabels = ['{:,.0f}'.format(x) + 'K' for x in g.get_xticks()/1000]
    g.set_xticklabels(xlabels)
    
    # Format with 2 decimal places
    #ax.get_yaxis().set_major_formatter(MathTextSciFormatter("%1.2e"))

    # Save figure
    fig.savefig('plots/'+str(filename),dpi=500)

def plotAbs_Par(output_01,output_10,output_11,output_12,output_13,output_14,y_ax, filename):
    # setup plot
    fig, ax = plt.subplots(figsize=(12,8))

    # create suitable np array
    output_10f = np.insert(output_10, 2, np.zeros(np.size(output_10,0)), axis=1)
    output_11f = np.insert(output_11, 2, np.ones(np.size(output_11,0)), axis=1)
    output_12f = np.insert(output_12, 2, np.full(np.size(output_12,0),2), axis=1)
    output_13f = np.insert(output_13, 2, np.full(np.size(output_13,0),3), axis=1)
    output_14f = np.insert(output_14, 2, np.full(np.size(output_14,0),4), axis=1)
    output_01f = np.insert(output_01, 2, np.full(np.size(output_01,0),5), axis=1)

    #output = np.concatenate((output_10f,output_12f,output_13f,output_14f))
    output = np.concatenate((output_10f,output_11f,output_12f,output_13f,output_14f))
    
    # convert to pd.DataFrame 
    data = pd.DataFrame(output, columns=['N','Count','Algorithm'])
    sequentialData = pd.DataFrame(output_01f, columns=['N','Count','Algorithm'])
    g=sns.lineplot(data=data, x='N', y='Count', hue='Algorithm', palette=['#c90076','#6a329f','#2986cc','#16537e','#6aa84f'], linewidth=.5, markers=True, dashes=False, markersize=1, markeredgewidth=0.2, ci='sd', style='Algorithm') #palette=['#38761d','#6a329f','#f1c232','#990000','#16537e']
    sns.lineplot(data=sequentialData, x='N', y='Count', hue='Algorithm', palette=['#f44336'], linewidth=.5, markers=False, dashes=[(6, 3)], markersize=1, markeredgewidth=0.2, ci='sd', style='Algorithm')

    # allow recreation of legend
    legend = ax.get_legend()
    handles = legend.legendHandles
    legend.remove()
    #ax.legend(handles, ['MPI Cell Matrix','Shared MPI Cell Matrix','OpenMP Cell Vector','OpenMP Integer Vectors'],loc='upper left',prop={'size': 15})
    ax.legend(handles, ['MPI Cell Matrix','MPI Integer Matrices','Shared MPI Integer Matrices','OpenMP Cell Vector','OpenMP Integer Vectors'],loc='upper left',prop={'size': 15})

    ax.annotate("Best Sequential Implementation", xy=(20000, output_01[21,1]), xycoords='data',
            xytext=(-40, 30), textcoords='offset points', ha='center',size=15, color='#351c75',
            bbox=dict(boxstyle="round4", fc="0.9",color='#b4a7d6'),
            arrowprops=dict(arrowstyle="->",color='#351c75'))
    
    # Get number max seq
    max_seq = int(np.amax(output_10[:,0]))

    # Get number of runs
    number_runs = np.count_nonzero(output_10[:,0] == max_seq)

    # Get number of total runs
    tot_runs = output_10.shape[0]

    runtime_10 = np.mean(output_10[tot_runs-number_runs:tot_runs,1])
    runtime_11 = np.mean(output_11[tot_runs-number_runs:tot_runs,1])
    runtime_12 = np.mean(output_12[tot_runs-number_runs:tot_runs,1])
    runtime_13 = np.mean(output_13[tot_runs-number_runs:tot_runs,1])
    runtime_14 = np.mean(output_14[tot_runs-number_runs:tot_runs,1])
    runtime_01 = np.mean(output_01[tot_runs-number_runs:tot_runs,1])

    ax.annotate(round(runtime_10,2), xy=(max_seq, runtime_10), xycoords='data',
            xytext=(10, 0), textcoords='offset points', ha='center',size=12, color='#38761d',
            bbox=dict(boxstyle="round4", fc="0.9",color='#38761d'))

    ax.annotate(round(runtime_11,2), xy=(max_seq, runtime_11), xycoords='data',
            xytext=(10, -10), textcoords='offset points', ha='center',size=12, color='#b45f06',
            bbox=dict(boxstyle="round4", fc="0.9",color='#b45f06'),
            arrowprops=dict(arrowstyle="->",color='#b45f06'))

    ax.annotate(round(runtime_12,2), xy=(max_seq, runtime_12), xycoords='data',
            xytext=(10, 0), textcoords='offset points', ha='center',size=12, color='#bf9000',
            bbox=dict(boxstyle="round4", fc="0.9",color='#bf9000'))

    ax.annotate(round(runtime_13,2), xy=(max_seq, runtime_13), xycoords='data',
            xytext=(10, -5), textcoords='offset points', ha='center',size=12, color='#990000',
            bbox=dict(boxstyle="round4", fc="0.9",color='#990000'))

    ax.annotate(round(runtime_14,2), xy=(max_seq, runtime_14), xycoords='data',
            xytext=(10, 5), textcoords='offset points', ha='center',size=12, color='#16537e',
            bbox=dict(boxstyle="round4", fc="0.9",color='#16537e'))

    ax.annotate(round(runtime_01,2), xy=(max_seq, runtime_01), xycoords='data',
            xytext=(10, 0), textcoords='offset points', ha='center',size=12, color='#351c75',
            bbox=dict(boxstyle="round4", fc="0.9",color='#351c75'))

    #set log-yaxis
    if(y_ax == 1):
        ax.set_yscale('log')

    # set title and format
    plt.title('[seconds]', loc='left',c='grey', fontsize='large')
    plt.suptitle('Parallel Absolute Runtimes',x=0.125, y=0.980, ha='left',fontsize='xx-large', fontweight=10)
    plt.xlabel('Sequence Length')
    ax.set(ylabel=None)

    xlabels = ['{:,.0f}'.format(x) + 'K' for x in g.get_xticks()/1000]
    g.set_xticklabels(xlabels)
    
    # Format with 2 decimal places
    #ax.get_yaxis().set_major_formatter(MathTextSciFormatter("%1.2e"))

    # Save figure
    fig.savefig('plots/'+str(filename),dpi=500)

def plotAbs_Par_long(output_10,output_11,output_12,y_ax, filename):
    # setup plot
    fig, ax = plt.subplots(figsize=(12,8))

    # create suitable np array
    output_10f = np.insert(output_10, 2, np.zeros(np.size(output_10,0)), axis=1)
    output_11f = np.insert(output_11, 2, np.ones(np.size(output_11,0)), axis=1)
    output_12f = np.insert(output_12, 2, np.full(np.size(output_12,0),2), axis=1)

    #output = np.concatenate((output_10f,output_12f,output_13f,output_14f))
    output = np.concatenate((output_10f,output_11f,output_12f))
    
    # convert to pd.DataFrame 
    data = pd.DataFrame(output, columns=['N','Count','Algorithm'])
    g=sns.lineplot(data=data, x='N', y='Count', hue='Algorithm', palette=['#c90076','#2986cc','#6aa84f'], linewidth=.5, markers=True, dashes=False, markersize=1, markeredgewidth=0.2, ci='sd', style='Algorithm') #palette=['#38761d','#6a329f','#f1c232','#990000','#16537e']

    # allow recreation of legend
    legend = ax.get_legend()
    handles = legend.legendHandles
    legend.remove()
    ax.legend(handles, ['MPI Cell Matrix','MPI Integer Matrices','Shared MPI Intger Matrices'],loc='upper left',prop={'size': 15})
    
    # Get number max seq
    max_seq = int(np.amax(output_10[:,0]))

    # Get number of runs
    number_runs = np.count_nonzero(output_10[:,0] == max_seq)

    # Get number of total runs
    tot_runs = output_10.shape[0]

    runtime_10 = np.mean(output_10[tot_runs-number_runs:tot_runs,1])
    runtime_11 = np.mean(output_11[tot_runs-number_runs:tot_runs,1])
    runtime_12 = np.mean(output_12[tot_runs-number_runs:tot_runs,1])

    ax.annotate(round(runtime_11,2), xy=(max_seq, runtime_11), xycoords='data',
            xytext=(10, -10), textcoords='offset points', ha='center',size=12, color='#b45f06',
            bbox=dict(boxstyle="round4", fc="0.9",color='#b45f06'))

    ax.annotate(round(runtime_12,2), xy=(max_seq, runtime_12), xycoords='data',
            xytext=(10, 0), textcoords='offset points', ha='center',size=12, color='#bf9000',
            bbox=dict(boxstyle="round4", fc="0.9",color='#bf9000'))


    ax.annotate(round(runtime_10,2), xy=(max_seq, runtime_10), xycoords='data',
            xytext=(10, 0), textcoords='offset points', ha='center',size=12, color='#38761d',
            bbox=dict(boxstyle="round4", fc="0.9",color='#38761d'))

    #set log-yaxis
    if(y_ax == 1):
        ax.set_yscale('log')

    # set title and format
    plt.title('[seconds]', loc='left',c='grey', fontsize='large')
    plt.suptitle('Parallel Absolute Runtimes',x=0.125, y=0.980, ha='left',fontsize='xx-large', fontweight=10)
    plt.xlabel('Sequence Length')
    ax.set(ylabel=None)


    xlabels = ['{:,.0f}'.format(x) + 'K' for x in g.get_xticks()/1000]
    g.set_xticklabels(xlabels)
    
    # Format with 2 decimal places
    #ax.get_yaxis().set_major_formatter(MathTextSciFormatter("%1.2e"))

    # Save figure
    fig.savefig('plots/'+str(filename),dpi=500)

def main():

        # retrieve data from output
        #output_00 = np.loadtxt(sys.argv[1],delimiter='\t')
        #output_01 = np.loadtxt(sys.argv[2],delimiter='\t')
        #output_02 = np.loadtxt(sys.argv[3],delimiter='\t') 
        #output_03 = np.loadtxt(sys.argv[4],delimiter='\t')
        #output_04 = np.loadtxt(sys.argv[5],delimiter='\t')
        #output_05 = np.loadtxt(sys.argv[6],delimiter='\t')
        #output_06 = np.loadtxt(sys.argv[7],delimiter='\t')
        output_10 = np.loadtxt(sys.argv[8],delimiter='\t')
        output_11 = np.loadtxt(sys.argv[9],delimiter='\t')
        output_12 = np.loadtxt(sys.argv[10],delimiter='\t')
        #output_13 = np.loadtxt(sys.argv[11],delimiter='\t')
        #output_14 = np.loadtxt(sys.argv[12],delimiter='\t')

        # linear y axis
        #plotAbs_Seq(output_00,output_01,output_02,output_03,output_04,output_05,output_06,0,'linAbsoluteRuntimeSeq')
        #plotAbs_Par(output_01,output_10[:,:2],output_11[:,:2],output_12[:,:2],output_13[:,:2],output_14[:,:2],0,'linAbsoluteRuntimePar')
        plotAbs_Par_long(output_10[:,:2],output_11[:,:2],output_12[:,:2],0,'linAbsoluteRuntimePar_long')

        # log y axis
        #plotAbs_Seq(output_00,output_01,output_02,output_03,output_04,output_05,output_06,1,'logAbsoluteRuntimeSeq')
        #plotAbs_Par(output_01,output_10[:,:2],output_11[:,:2],output_12[:,:2],output_13[:,:2],output_14[:,:2],1,'logAbsoluteRuntimePar')
        plotAbs_Par_long(output_10[:,:2],output_11[:,:2],output_12[:,:2],1,'logAbsoluteRuntimePar_long')

    
if __name__ == '__main__':
    main()
