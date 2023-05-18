#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['legend.frameon'] = False
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['axes.grid'] = False
mpl.rcParams['axes.xmargin'] = 0.01
mpl.rcParams['axes.ymargin'] = 0.01
mpl.rcParams['grid.alpha'] = 1.0
mpl.rcParams['grid.linestyle'] = 'dotted'
mpl.rcParams['grid.color'] = 'gray'
mpl.rcParams['grid.linewidth'] = 0.5
mpl.rcParams['xtick.major.size'] = 5.0
mpl.rcParams['ytick.major.size'] = 5.0
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.monospace'] = 'Helvetica'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True 
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D', '#666666'])

# current working directory
cwd = os.getcwd()

# directory containing data
output_dir = cwd + "/OUTPUTS/"

# bubble size and slip velocity distribution
plt.figure(figsize=[14,15])
for i in range(10,51,10):
    fname = output_dir + "QGT-" + str(i) + ".dat"
    df = pd.read_csv(fname, delimiter='\t')
    x = np.asarray(df['dbi[mm]'])
    y = np.asarray(df['pdf'])
    plt.subplot(3,2,1)
    plt.plot(x, y, label=r"$Q_{gt}$ = " + str((i)/100.0) + r" $m^3/s$")
    plt.xlabel(r'Volume equivalent bubble diameter, $d_{bi}$ [ mm ]')
    plt.ylabel(r'Probability density')
    plt.legend()
    plt.subplot(3,2,2)
    x = np.asarray(df['usi[m/s]'])
    plt.plot(x, y, label=r"$Q_{gt}$ = " + str((i)/100.0) + r" $m^3/s$")
    plt.xlabel(r'Bubble slip velocity, $u_{si}$ [ m/s ]')
    plt.ylabel(r'Probability density')
    plt.legend()

# summary plots
fname = output_dir + "/summary_results.dat"
df = pd.read_csv(fname, delimiter='\t')

plt.subplot(3,2,3)
x = df['Ugt[m/s]']
y = df['egfb']
y2 = df['egbed']
plt.plot(x, y, label = 'freeboard')
plt.plot(x,y2, label = 'bed')
plt.xlabel(r'Inlet gas velocity [ m/s ]')
plt.ylabel(r'Gas holdup')
plt.legend()

plt.subplot(3,2,4)
y = df['epsilonl']
plt.plot(x,y)
plt.xlabel(r'Inlet gas velocity [ m/s ]')
plt.ylabel(r'Bed liquid holdup')

plt.subplot(3,2,5)
y = df['R']
plt.plot(x, y)
plt.xlabel('Inlet gas velocity [ m/s ]')
plt.ylabel('Liquid Recycle Ratio [ - ]')

plt.subplot(3,2,6)
y = df['Rgas']
plt.plot(x, y)
plt.xlabel('Inlet gas velocity [ m/s ]')
plt.ylabel('Recycled Gas Fraction')
plt.savefig('FIGURES/all_plots.png')