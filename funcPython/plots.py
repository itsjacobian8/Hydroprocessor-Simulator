#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib import rcParams, cycler

rcParams['legend.frameon'] = False
rcParams['axes.linewidth'] = 1.5
rcParams['axes.grid'] = True
rcParams['axes.xmargin'] = 0.01
rcParams['axes.ymargin'] = 0.01
rcParams['grid.alpha'] = 1.0
rcParams['grid.linestyle'] = 'dotted'
rcParams['grid.color'] = 'gray'
rcParams['grid.linewidth'] = 0.5
rcParams['xtick.major.size'] = 5.0
rcParams['ytick.major.size'] = 5.0
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['font.family'] = 'sans-serif'
rcParams['font.monospace'] = 'libris'
rcParams['xtick.top'] = True
rcParams['ytick.right'] = True 
rcParams['lines.linewidth'] = 3
rcParams['axes.prop_cycle'] = cycler(color=['#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D', '#666666'])

def plot_summary(df, system, recycleInfo, fname):
    """
    plots summary results
    
    args:
        df: pandas dataframe
        system: specifies whether the data belongs to bubble column or ebullated bed (boolean)
        recycleInfo: specifies internal recycle (boolean)
    """
    if(system.startswith('bubble') and system.endswith('column')):
        if(recycleInfo):
            plt.figure(figsize=[14,10])
            plt.subplot(2,2,1)
            plt.plot(df['UgIn'], df['eg'])
            plt.xlabel(r'inlet gas velocity (m/s)')
            plt.ylabel(r'overall gas holdup')
            
            plt.subplot(2,2,2)
            plt.plot(df['UgIn'], df['eta'])
            plt.xlabel(r'inlet gas velocity (m/s)')
            plt.ylabel(r'separation efficiency')
            
            plt.subplot(2,2,3)
            plt.plot(df['UgIn'], df['qgr'], label=r'gas')
            plt.plot(df['UgIn'], df['qlr'], label=r'liquid')
            plt.xlabel(r'inlet gas velocity (m/s)')
            plt.ylabel(r'recycle flow rates ($m^3/s$)')
            plt.legend()
            
            plt.subplot(2,2,4)
            plt.plot(df['UgIn'], df['Rgas'], label=r'gas')
            plt.plot(df['UgIn'], df['R'], label=r'liquid')
            plt.xlabel(r'inlet gas velocity (m/s)')
            plt.ylabel(r'recycled fractions')
            plt.legend()
            
            plt.savefig(fname + ".png", dpi=300)
            
        else:
            plt.figure(figsize=[7,5])
            plt.plot(df['Ug'], df['eg'])
            plt.xlabel(r'superficial gas velocity')
            plt.ylabel(r'overall gas holdup')
            plt.savefig(fname + ".png", dpi=300)
            
    elif(system.startswith('ebullated') and system.endswith('bed')):
        if(recycleInfo):
            plt.figure(figsize = [14,15])
            plt.subplot(3,2,1)
            plt.plot(df['Ugt'], df['egbed'], label=r'gas')
            plt.plot(df['Ugt'], df['epsilonl'], label=r'liquid')
            plt.plot(df['Ugt'], df['epsilonp'], label=r'solids')
            plt.xlabel(r'inlet gas velocity (m/s)')
            plt.ylabel(r'bed holdup')
            plt.legend()
            
            plt.subplot(3,2,2)
            plt.plot(df['Ugt'], df['egfb'])
            plt.xlabel(r'inlet gas velocity (m/s)')
            plt.ylabel(r'freeboard gas holdup')
            
            plt.subplot(3,2,3)
            plt.plot(df['Ugt'], df['eta'])
            plt.xlabel(r'inlet gas velocity (m/s)')
            plt.ylabel(r'separation efficiency')
            
            plt.subplot(3,2,4)
            plt.plot(df['Ugt'], df['qgr'], label=r'gas')
            plt.plot(df['Ugt'], df['qlr'], label=r'liquid')
            plt.xlabel(r'inlet gas velocity (m/s)')
            plt.ylabel(r'recycle flow rates ($m^3/s$)')
            plt.legend()
            
            plt.subplot(3,2,5)
            plt.plot(df['Ugt'], df['Rgas'], label=r'gas')
            plt.plot(df['Ugt'], df['R'], label=r'liquid')
            plt.xlabel(r'inlet gas velocity (m/s)')
            plt.ylabel(r'recycle fractions')
            plt.legend()
            
            plt.subplot(3,2,6)
            plt.plot(df['Ugt'], df['Umf'])
            plt.xlabel(r'inlet gas velocity (m/s)')
            plt.ylabel(r'minimum fluidization velocity (m/s)')
            
            plt.savefig(fname + ".png", dpi=300)
        else:
            plt.figure(figsize=[14,10])
            plt.subplot(2,2,1)
            plt.plot(df['Ug'], df['egbed'], label=r'gas')
            plt.plot(df['Ug'], df['epsilonl'], label=r'liquid')
            plt.plot(df['Ug'], df['epsilonp'], label=r'solid')
            plt.xlabel(r'superficial gas velocity (m/s)')
            plt.ylabel(r'bed holdups')
            plt.legend()
            
            plt.subplot(2,2,2)
            plt.plot(df['Ug'], df['egfb'])
            plt.xlabel(r'superficial gas velocity (m/s)')
            plt.ylabel(r'freeboard gas holdup')
            
            plt.subplot(2,2,3)
            plt.plot(df['Ug'], df['hbed'])
            plt.xlabel(r'superficial gas velocity (m/s)')
            plt.ylabel(r'bed height')
            
            plt.subplot(2,2,4)
            plt.plot(df['Ug'], df['Umf'])
            plt.xlabel(r'superficial gas velocity (m/s)')
            plt.ylabel(r'minimum fluidization velocity (m/s)')
            
            plt.savefig(fname + ".png", dpi=300)
            
def plot_bsvd(df, fname):
    """
    plots bubble size distribution
    
    args:
        df: pandas data frame
        fname: file name for storing the plots
    """
    plt.figure(figsize=[7,10])
    plt.subplot(2,1,1)
    plt.plot(df['dbi'], df['pdf'])
    plt.xlabel(r'volume equivalent diameter (mm)')
    plt.ylabel(r'probability density')
    
    plt.subplot(2,1,2)
    plt.plot(df['usi'], df['pdf'])
    plt.xlabel(r'bubble slip velocity (m/s)')
    plt.ylabel(r'probability density')
    plt.savefig(fname + ".png", dpi=300)  