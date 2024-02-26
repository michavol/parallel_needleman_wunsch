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


def linePlot(output, filename):
    # setup plot
    fig, ax = plt.subplots(figsize=(12,8))
    
    # convert to pd.DataFrame 
    data = pd.DataFrame(output, columns=['N','Count','Processors','BlockLength'])
    sns.lineplot(data=data, x='BlockLength', y='Count', hue='Processors', palette="tab10", linewidth=.5, markers=True, dashes=False, markersize=1, markeredgewidth=0.2, ci=99)
    
    # allow recreation of legend
    # legend = ax.get_legend()
    # handles = legend.legendHandles
    # legend.remove()
    # ax.legend(handles, ['Sequential Baseline', 'Sequential Optimized', 'Sequential Diagonal' ,'Parallel Baseline']) #'Sequential Diagonal'
    #ax.legend(handles, [sys.argv[1], sys.argv[2], sys.argv[3]])
    
    #set log-yaxis
    #ax.set_yscale('log')
    
    # set title and format
    plt.title('[seconds]', loc='left',c='grey', fontsize='large')
    plt.suptitle('Absolute Runtime',x=0.125, y=0.980, ha='left',fontsize='xx-large', fontweight=10)
    plt.xlabel('Block Size')
    ax.set(ylabel=None)
    
    # Format with 2 decimal places
    #ax.get_yaxis().set_major_formatter(MathTextSciFormatter("%1.2e"))

    # Save figure
    fig.savefig('plots/'+str(filename),dpi=500)
  

def main():
    # retrieve data from output
    output_2 = np.loadtxt(sys.argv[1],delimiter='\t')

    linePlot(output_2,'blockLength_vs_cycles.png')


    
if __name__ == '__main__':
    main()
