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

def plotLinePlot_all(output_10,output_13,output_14,scale,filename):
    # setup plot
    fig, ax = plt.subplots(figsize=(12,8))

    # Get number of runs
    number_runs = np.count_nonzero(output_10[:,2] == 1)

    # Get number of processors
    number_proc = int(np.amax(output_10[:,2]))

    # linear speed up for all individual algorithms
    lin_speed_10 = np.copy(output_10)
    lin_speed_10[:number_runs,1] = np.mean(lin_speed_10[:number_runs,1])
    # lin_speed_11 = np.copy(output_11)
    # lin_speed_11[:number_runs,1] = np.mean(lin_speed_11[:number_runs,1])
    # lin_speed_12 = np.copy(output_12)
    # lin_speed_12[:number_runs,1] = np.mean(lin_speed_12[:number_runs,1])
    lin_speed_13 = np.copy(output_13)
    lin_speed_13[:number_runs,1] = np.mean(lin_speed_13[:number_runs,1])
    lin_speed_14 = np.copy(output_14)
    lin_speed_14[:number_runs,1] = np.mean(lin_speed_14[:number_runs,1])

    start = number_runs

    for i in range (2,number_proc+1):
        lin_speed_10[start:start+number_runs,1] = np.divide(lin_speed_10[:number_runs,1],i)
        # lin_speed_11[start:start+number_runs,1] = np.divide(lin_speed_11[:number_runs,1],i)
        # lin_speed_12[start:start+number_runs,1] = np.divide(lin_speed_12[:number_runs,1],i)
        lin_speed_13[start:start+number_runs,1] = np.divide(lin_speed_13[:number_runs,1],i)
        lin_speed_14[start:start+number_runs,1] = np.divide(lin_speed_14[:number_runs,1],i)
        start = start + number_runs
    

    
    output_10f = np.insert(output_10, 3, np.zeros(np.size(output_10,0)), axis=1)
    # output_11f = np.insert(output_11, 3, np.ones(np.size(output_11,0)), axis=1)
    # output_12f = np.insert(output_12, 3, np.full(np.size(output_12,0),2), axis=1)
    output_13f = np.insert(output_13, 3, np.full(np.size(output_13,0),3), axis=1)
    output_14f = np.insert(output_14, 3, np.full(np.size(output_14,0),4), axis=1)
    output_linspeed_10f = np.insert(lin_speed_10, 3, np.full(np.size(lin_speed_10,0),5), axis=1)
    #output_linspeed_11f = np.insert(lin_speed_11, 3, np.full(np.size(lin_speed_11,0),6), axis=1)
    # output_linspeed_12f = np.insert(lin_speed_12, 3, np.full(np.size(lin_speed_12,0),7), axis=1)
    output_linspeed_13f = np.insert(lin_speed_13, 3, np.full(np.size(lin_speed_13,0),8), axis=1)
    output_linspeed_14f = np.insert(lin_speed_14, 3, np.full(np.size(lin_speed_14,0),9), axis=1)
    

    output = np.concatenate((output_10f,output_13f,output_14f, output_linspeed_10f,output_linspeed_13f,output_linspeed_14f))
    #output = np.concatenate((output_10f,output_11f,output_13f,output_14f, output_linspeed_10f,output_linspeed_11f,output_linspeed_13f,output_linspeed_14f))
    #output = np.concatenate((output_10f,output_11f,output_12f,output_13f,output_14f, output_linspeed_10f,output_linspeed_11f,output_linspeed_12f,output_linspeed_13f,output_linspeed_14f))

    # convert to pd.DataFrame 
    data = pd.DataFrame(output, columns=['Sequence Length','Count','Processors','Algorithm'])

    sns.lineplot(data=data, x='Processors',y='Count',hue='Algorithm',palette=['#38761d','#990000','#16537e','#81d65c','#f44336','#2986cc'], linewidth=.5, markers=True, dashes=False, markersize=1, markeredgewidth=0.2, ci='sd') #palette=['#38761d','#6a329f','#741b47','#990000','#16537e','#81d65c','#674ea7','#c90076','#f44336','#2986cc']

    
    # allow recreation of legend
    legend = ax.get_legend()
    handles = legend.legendHandles
    legend.remove()
    ax.legend(handles, ['MPI Cell Matrix','OpenMP Cell Vector','OpenMP Integer Vectors','Lin Speedup for MPI Cell Matrix','Lin Speedup for OpenMP Cell Vector','Lin Speedup for OpenMP Integer Vectors'],loc='upper right',prop={'size': 15})
    #ax.legend(handles, ['MPI Cell Matrix','MPI Integer Matrcies','OpenMP Cell Vector','OpenMP Integer Vectors','Lin Speedup for MPI Cell Matrix','Lin Speedup for MPI Integer Matrices','Lin Speedup for OpenMP Cell Vector','Lin Speedup for OpenMP Integer Vectors'],prop={'size': 15})
    #ax.legend(handles, ['MPI Cell Matrix','MPI Integer Matrcies','Shared MPI Integer Matrcies','OpenMP Cell Vector','OpenMP Integer Vectors','Lin Speedup for MPI Cell Matrix','Lin Speedup for MPI Integer Matrices','Lin Speedup Shared MPI Integer Matrcies','Lin Speedup for OpenMP Cell Vector','Lin Speedup for OpenMP Integer Vectors'],prop={'size': 15})

    #set log-yaxis
    #ax.set_yscale('log')
    
    # set title and format
    plt.title('[seconds]', loc='left',c='grey',fontsize='large')
    plt.suptitle('Processors vs. Time',x=0.125, y=0.98, ha='left',fontsize='xx-large',fontweight=500)
    plt.xlabel('Number Processors')
    plt.ticklabel_format(axis='y',style='sci', scilimits=(-2,2))
    ax.set(ylabel=None)
    

    if(scale == 1):
        ax.set_yscale('log')
    else:
        plt.ticklabel_format(style='plain', axis='y')
    # Format with 2 decimal places
    #ax.get_yaxis().set_major_formatter(MathTextSciFormatter("%1.2e"))
    #plt.ticklabel_format(style='plain', axis='y')

    # Format x-axis
    ax.xaxis.set_ticks(np.arange(1, number_proc+1, 1.))

    plt.xlim([0, 36])

    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2))
    ax.xaxis.set_major_formatter('{x:.0f}')
    

    # Save figure
    fig.savefig('plots/'+str(filename),dpi=500)

def plotLinePlot_MPI(output_10,output_11,scale,filename):
    # setup plot
    fig, ax = plt.subplots(figsize=(12,8))

    # Get number of runs
    number_runs = np.count_nonzero(output_10[:,2] == 1)

    # Get number of processors
    number_proc = int(np.amax(output_10[:,2]))

    # linear speed up for all individual algorithms
    lin_speed_10 = np.copy(output_10)
    lin_speed_10[:number_runs,1] = np.mean(lin_speed_10[:number_runs,1])
    lin_speed_11 = np.copy(output_11)
    lin_speed_11[:number_runs,1] = np.mean(lin_speed_11[:number_runs,1])
    # lin_speed_12 = np.copy(output_12)
    # lin_speed_12[:number_runs,1] = np.mean(lin_speed_12[:number_runs,1])


    start = number_runs

    for i in range (2,number_proc+1):
        lin_speed_10[start:start+number_runs,1] = np.divide(lin_speed_10[:number_runs,1],i)
        lin_speed_11[start:start+number_runs,1] = np.divide(lin_speed_11[:number_runs,1],i)
        # lin_speed_12[start:start+number_runs,1] = np.divide(lin_speed_12[:number_runs,1],i)
        start = start + number_runs
    

    
    output_10f = np.insert(output_10, 3, np.zeros(np.size(output_10,0)), axis=1)
    output_11f = np.insert(output_11, 3, np.ones(np.size(output_11,0)), axis=1)
    # output_12f = np.insert(output_12, 3, np.full(np.size(output_12,0),2), axis=1)

    output_linspeed_10f = np.insert(lin_speed_10, 3, np.full(np.size(lin_speed_10,0),5), axis=1)
    output_linspeed_11f = np.insert(lin_speed_11, 3, np.full(np.size(lin_speed_11,0),6), axis=1)
    # output_linspeed_12f = np.insert(lin_speed_12, 3, np.full(np.size(lin_speed_12,0),7), axis=1)

    output = np.concatenate((output_10f,output_11f, output_linspeed_10f,output_linspeed_11f))
    #output = np.concatenate((output_10f,output_11f,output_13f,output_14f, output_linspeed_10f,output_linspeed_11f,output_linspeed_13f,output_linspeed_14f))
    #output = np.concatenate((output_10f,output_11f,output_12f,output_13f,output_14f, output_linspeed_10f,output_linspeed_11f,output_linspeed_12f,output_linspeed_13f,output_linspeed_14f))

    # convert to pd.DataFrame 
    data = pd.DataFrame(output, columns=['Sequence Length','Count','Processors','Algorithm'])

    sns.lineplot(data=data, x='Processors',y='Count',hue='Algorithm',palette=['#38761d','#6a329f','#81d65c','#674ea7'], linewidth=.5, markers=True, dashes=False, markersize=1, markeredgewidth=0.2, ci='sd') #palette=['#38761d','#6a329f','#741b47','#990000','#16537e','#81d65c','#674ea7','#c90076','#f44336','#2986cc']

    
    # allow recreation of legend
    legend = ax.get_legend()
    handles = legend.legendHandles
    legend.remove()
    ax.legend(handles, ['MPI Cell Matrix','MPI Integer Matrcies','Lin Speedup for MPI Cell Matrix','Lin Speedup for MPI Integer Matrcies'],loc='upper right',prop={'size': 15})
    #ax.legend(handles, ['MPI Cell Matrix','MPI Integer Matrcies','OpenMP Cell Vector','OpenMP Integer Vectors','Lin Speedup for MPI Cell Matrix','Lin Speedup for MPI Integer Matrices','Lin Speedup for OpenMP Cell Vector','Lin Speedup for OpenMP Integer Vectors'],prop={'size': 15})
    #ax.legend(handles, ['MPI Cell Matrix','MPI Integer Matrcies','Shared MPI Integer Matrcies','OpenMP Cell Vector','OpenMP Integer Vectors','Lin Speedup for MPI Cell Matrix','Lin Speedup for MPI Integer Matrices','Lin Speedup Shared MPI Integer Matrcies','Lin Speedup for OpenMP Cell Vector','Lin Speedup for OpenMP Integer Vectors'],prop={'size': 15})

    #set log-yaxis
    #ax.set_yscale('log')
    
    # set title and format
    plt.title('[seconds]', loc='left',c='grey',fontsize='large')
    plt.suptitle('Processors vs. Time for MPI',x=0.125, y=0.98, ha='left',fontsize='xx-large',fontweight=500)
    plt.xlabel('Number Processors')
    plt.ticklabel_format(axis='y',style='sci', scilimits=(-2,2))
    ax.set(ylabel=None)
    

    if(scale == 1):
        ax.set_yscale('log')
    else:
        plt.ticklabel_format(style='plain', axis='y')
    # Format with 2 decimal places
    #ax.get_yaxis().set_major_formatter(MathTextSciFormatter("%1.2e"))
    #plt.ticklabel_format(style='plain', axis='y')

    # Format x-axis
    ax.xaxis.set_ticks(np.arange(1, number_proc+1, 1.))

    plt.xlim([0, 36])

    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2))
    ax.xaxis.set_major_formatter('{x:.0f}')
    

    # Save figure
    fig.savefig('plots/'+str(filename),dpi=500)

def plotLinePlot_OpenMP(output_13,output_14,scale,filename):
    # setup plot
    fig, ax = plt.subplots(figsize=(12,8))

    # Get number of runs
    number_runs = np.count_nonzero(output_13[:,2] == 1)

    # Get number of processors
    number_proc = int(np.amax(output_13[:,2]))

    # linear speed up for all individual algorithms
    lin_speed_13 = np.copy(output_13)
    lin_speed_13[:number_runs,1] = np.mean(lin_speed_13[:number_runs,1])
    lin_speed_14 = np.copy(output_14)
    lin_speed_14[:number_runs,1] = np.mean(lin_speed_14[:number_runs,1])

    start = number_runs

    for i in range (2,number_proc+1):
        lin_speed_13[start:start+number_runs,1] = np.divide(lin_speed_13[:number_runs,1],i)
        lin_speed_14[start:start+number_runs,1] = np.divide(lin_speed_14[:number_runs,1],i)
        start = start + number_runs
    

    output_13f = np.insert(output_13, 3, np.full(np.size(output_13,0),3), axis=1)
    output_14f = np.insert(output_14, 3, np.full(np.size(output_14,0),4), axis=1)
    output_linspeed_13f = np.insert(lin_speed_13, 3, np.full(np.size(lin_speed_13,0),8), axis=1)
    output_linspeed_14f = np.insert(lin_speed_14, 3, np.full(np.size(lin_speed_14,0),9), axis=1)
    

    output = np.concatenate((output_13f,output_14f,output_linspeed_13f,output_linspeed_14f))
    #output = np.concatenate((output_10f,output_11f,output_13f,output_14f, output_linspeed_10f,output_linspeed_11f,output_linspeed_13f,output_linspeed_14f))
    #output = np.concatenate((output_10f,output_11f,output_12f,output_13f,output_14f, output_linspeed_10f,output_linspeed_11f,output_linspeed_12f,output_linspeed_13f,output_linspeed_14f))

    # convert to pd.DataFrame 
    data = pd.DataFrame(output, columns=['Sequence Length','Count','Processors','Algorithm'])

    sns.lineplot(data=data, x='Processors',y='Count',hue='Algorithm',palette=['#990000','#16537e','#f44336','#2986cc'], linewidth=.5, markers=True, dashes=False, markersize=1, markeredgewidth=0.2, ci='sd') #palette=['#38761d','#6a329f','#741b47','#990000','#16537e','#81d65c','#674ea7','#c90076','#f44336','#2986cc']

    
    # allow recreation of legend
    legend = ax.get_legend()
    handles = legend.legendHandles
    legend.remove()
    ax.legend(handles, ['OpenMP Cell Vector','OpenMP Integer Vectors','Lin Speedup for OpenMP Cell Vector','Lin Speedup for OpenMP Integer Vectors'],loc='upper right',prop={'size': 15})
    #ax.legend(handles, ['MPI Cell Matrix','MPI Integer Matrcies','OpenMP Cell Vector','OpenMP Integer Vectors','Lin Speedup for MPI Cell Matrix','Lin Speedup for MPI Integer Matrices','Lin Speedup for OpenMP Cell Vector','Lin Speedup for OpenMP Integer Vectors'],prop={'size': 15})
    #ax.legend(handles, ['MPI Cell Matrix','MPI Integer Matrcies','Shared MPI Integer Matrcies','OpenMP Cell Vector','OpenMP Integer Vectors','Lin Speedup for MPI Cell Matrix','Lin Speedup for MPI Integer Matrices','Lin Speedup Shared MPI Integer Matrcies','Lin Speedup for OpenMP Cell Vector','Lin Speedup for OpenMP Integer Vectors'],prop={'size': 15})

    #set log-yaxis
    #ax.set_yscale('log')
    
    # set title and format
    plt.title('[seconds]', loc='left',c='grey',fontsize='large')
    plt.suptitle('Processors vs. Time for OpenMP',x=0.125, y=0.98, ha='left',fontsize='xx-large',fontweight=500)
    plt.xlabel('Number Processors')
    plt.ticklabel_format(axis='y',style='sci', scilimits=(-2,2))
    ax.set(ylabel=None)
    

    if(scale == 1):
        ax.set_yscale('log')
    else:
        plt.ticklabel_format(style='plain', axis='y')
    # Format with 2 decimal places
    #ax.get_yaxis().set_major_formatter(MathTextSciFormatter("%1.2e"))
    #plt.ticklabel_format(style='plain', axis='y')

    # Format x-axis
    ax.xaxis.set_ticks(np.arange(1, number_proc+1, 1.))
    
    plt.ylim([0,110])
    plt.xlim([0, 36])

    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2))
    ax.xaxis.set_major_formatter('{x:.0f}')
    

    # Save figure
    fig.savefig('plots/'+str(filename),dpi=500)

def boxPlot(output, filename):
    # setup plot
    fig, ax = plt.subplots(figsize=(12,8))
    
    data = pd.DataFrame(output, columns=['Sequence Length','Count','Processors'])

    sns.boxplot(data=data, x='Processors',y='Count',hue='Sequence Length',palette="tab10")


    plt.title('[seconds]', loc='left',c='grey',fontsize='large')
    plt.suptitle('Processors vs. Time',x=0.125, y=0.96, ha='left',fontsize='xx-large',fontweight=500)
    plt.xlabel('Number Processors')
    plt.ticklabel_format(axis='y',style='sci', scilimits=(-2,2))
    ax.set(ylabel=None)

    # Format with 2 decimal places
    ax.get_yaxis().set_major_formatter(MathTextSciFormatter("%1.2e"))

    fig.align_ylabels()
    fig.savefig('plots/'+str(filename),dpi=500)
  

def main():
    # retrieve data from output
    output_10 = np.loadtxt(sys.argv[1],delimiter='\t')
    output_11 = np.loadtxt(sys.argv[2],delimiter='\t')
    #output_12 = np.loadtxt(sys.argv[1],delimiter='\t')
    output_13 = np.loadtxt(sys.argv[3],delimiter='\t')
    output_14 = np.loadtxt(sys.argv[4],delimiter='\t')

    #boxPlot(output_2,'boxPlot_procs.png')
    #plotLinePlot_all(output_10[:,:3],output_13[:,:3],output_14[:,:3],0, 'linLinePlot_procs.png')
    #plotLinePlot_all(output_10[:,:3],output_13[:,:3],output_14[:,:3],1, 'logLinePlot_procs.png')

    plotLinePlot_OpenMP(output_13[:,:3],output_14[:,:3],0, 'linLinePlotOpenMP_procs.png')
    plotLinePlot_OpenMP(output_13[:,:3],output_14[:,:3],1, 'logLinePlotOpenMP_procs.png')
    #plotLinePlot(output_13[:,:3], 'linePlot_threads.png','Parallel OpenMP')

    #plotLinePlot_MPI(output_10[:,:3],output_11[:,:3],0, 'linLinePlotMPI_procs.png')
    #plotLinePlot_MPI(output_10[:,:3],output_11[:,:3],1, 'logLinePlotMPI_procs.png')


    
if __name__ == '__main__':
    main()
