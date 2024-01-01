#!/usr/bin/env python3
from getInputs import GetInputs
from solvers import Solvers
from time import ctime, time
from termcolor import colored
import pandas as pd
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
rcParams['font.monospace'] = 'Alegreya'
rcParams['xtick.top'] = True
rcParams['ytick.right'] = True 
rcParams['lines.linewidth'] = 3
rcParams['axes.prop_cycle'] = cycler(color=['#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D', '#666666'])

class Callers(GetInputs):
    """
    This is a collection of miscellaneous methods for:
        - obtaining input parameters from the user
        - making call to main solver for each inlet gas flow rate
        - writing results to file
        - plotting the results
    
    class inherits from GetInputs class
    
    class attributes
    -----------------
    this class has no attributes and can be instantiated without any inputs
    
    class methods
    --------------
    callToInputs: None
        makes call to main methods for obtaining input parameters from user
    
    callToSolvers: None
        makes call to main solver
    
    plotNoRecycleBubble: None
        plots summary results for bubble column with no internal recycle
    
    plotRecycleBubble: None
        plots summary results for bubble column with constant and variable internal recycle
    
    plotNoRecycleEbullated: None
        plots summary results for ebullated bed with no internal recycle
    
    plotVariableRecycleEbullated: None
        plots summary results for ebullated bed with variable internal recycle
    
    plot_bsvd: None
        plots bubble size and slip velocity distribution
    
    for more information, consult docstrings for individual methods
    """
    def __init__(self):
        super().__init__()
    
    def callToInputs(self):
        """
        makes call to main functions for obtaining input parameters
        """
        # obtain system information
        self.getSystem()
        
        # obtain recycle information
        self.getRecycleInfo()
        
        # obtain inlet flow rates
        self.getInletFlowRates()
        
        # obtain gas-liquid properties
        self.getGasLiquidProperties()
        
        # obtain solid properties
        if(self.system.startswith("ebullated") and self.system.endswith("bed")):
            self.getSolidProperties()
        else:
            self.inputs['Vp'] = 0.0
            self.inputs['Ap'] = 0.0
            self.inputs['rhop'] = 0.0
            self.inputs['mcat'] = 0.0
            self.inputs['hset'] = 0.0
        
        # obtain grid parameters
        self.getGridParameters()
        
        # obtain column parameters
        self.getColumnParameters()
        
        # obtain convergence criteria
        self.getConvergenceCriteria()
        
    def callToSolvers(self):
        """
        makes call to main solver for each inlet gas flow rate
        """
        
        # obtain input parameters
        self.callToInputs()
        sysInfo = self.system
        recInfo = self.recycleInfo
        
        # number of inlet gas flow rates
        nGasFlows = int(round((self.inputs['qgtMax'] - self.inputs['qgtMin'])/self.inputs['qgtIncrement'] + 1))
        
        summary_list = []
        for i in range(nGasFlows):
            # next gas flow rate
            self.inputs['qgt'] = self.inputs['qgtMin'] + float(i)*self.inputs['qgtIncrement']
            print(colored(f"[{ctime(time())}] SIMULATING INLET FLOW RATE = {self.inputs['qgt']:.2f} m^3/s", "cyan", "on_black"))
            # initialize solver
            solve = Solvers(self.inputs)
            
            # make specific solver calls
            if(sysInfo.startswith("bubble") and sysInfo.endswith("column")):
                if(recInfo.startswith("no") and recInfo.endswith("recycle")):
                    summary, bsvd = solve.noRecycleBubbleColumn()
                elif(recInfo.startswith("constant") and recInfo.endswith("recycle")):
                    summary, bsvd = solve.constantRecycleBubbleColumn(self.inputs['R'])
            elif(sysInfo.startswith("ebullated") and sysInfo.endswith("bed")):
                if(recInfo.startswith("no") and recInfo.endswith("recycle")):
                    summary, bsvd = solve.noRecycleEbullatedBed()
                elif(recInfo.startswith("variable") and recInfo.endswith("recycle")):
                    summary, bsvd = solve.variableRecycleEbullatedBed()
            
            # accumulate summary results
            summary_list.append(summary)
            
            # write bubble size and slip velocity distribution to file
            df = pd.DataFrame(bsvd)
            fname = f"bsvd-qgt-{self.inputs['qgt']:.2f}"
            df.to_csv(fname + ".txt", sep='\t')
            self.plot_bsvd(df, fname)
            
        # write summary results to file
        summary_df = pd.DataFrame(summary_list)
        fname = f"summary-qlvr-{self.inputs['qlvr']:.4f}"
        summary_df.to_csv(fname + ".txt", sep='\t')
        
        # plot summary results
        if(sysInfo.startswith("bubble") and sysInfo.endswith("column")):
            if(recInfo.startswith("no") and recInfo.endswith("recycle")):
                self.plotNoRecycleBubble(summary_df, fname)
            else:
                self.plotRecycleBubble(summary_df, fname)
        elif(sysInfo.startswith("ebullated") and sysInfo.endswith("bed")):
            if(recInfo.startswith("no") and recInfo.endswith("recycle")):
                self.plotNoRecycleEbullated(summary_df, fname)
            elif(recInfo.startswith("variable") and recInfo.endswith("recycle")):
                self.plotVariableRecycleEbullated(summary_df, fname)
            
    def plotNoRecycleBubble(self, df, fname):
        """
        plots summary results for no recycle bubble column
        
        args:
            df: data frame containing summary results
            fname: file name without extension
        """
        plt.figure(figsize=[7,5])
        plt.plot(df['Ugt'], df['egfb'])
        plt.xlabel(r"superficial velocity [ m/s ]")
        plt.ylabel(r"overall gas holdup [ - ]")
        plt.savefig(fname + ".png", dpi=300)
        plt.close()
    def plotRecycleBubble(self, df, fname):
        """
        plots summary results for constant and variable recycle bubble column
        
        args:
            df: data frame containing summary results
            fname: file name without extension
        """
        plt.figure(figsize=[14, 10])
        plt.subplot(2,2,1)
        plt.plot(df['Ugt'], df['egfb'])
        plt.xlabel(r"inlet gas velocity [ m/s ]")
        plt.ylabel(r"overall gas holdup [ - ]")
        
        plt.subplot(2,2,2)
        plt.plot(df['Ugt'], df['eta'])
        plt.xlabel(r"inlet gas velocity [ m/s ]")
        plt.ylabel(r"gas-liquid separation efficiency [ - ]")
        
        plt.subplot(2,2,3)
        plt.plot(df['Ugt'], df['qgr'], label='gas')
        plt.plot(df['Ugt'], df['qlr'], label='liquid')
        plt.xlabel(r"inlet gas velocity [ m/s ]")
        plt.ylabel(r"recycled flow rates [ $m^3/s$]")
        plt.legend()
        
        plt.subplot(2,2,4)
        plt.plot(df['Ugt'], df['Rgas'], label='gas')
        plt.plot(df['Ugt'], df['R'], label='liquid')
        plt.xlabel(r"inlet gas velocity [ m/s ]")
        plt.ylabel(r"recycled fractions [ - ]")
        plt.legend()
        plt.savefig(fname + ".png", dpi=300)
        plt.close()
    def plotNoRecycleEbullated(self, df, fname):
        """
        plots summary results for no recycle ebullated bed
        
        args:
            df: data frame containing summary results
            fname: file name without extension
        """
        plt.figure(figsize=[14,10])
        plt.subplot(2,2,1)
        plt.plot(df['Ugt'], df['egbed'], label='gas')
        plt.plot(df['Ugt'], df['epsilonl'], label='liquid')
        plt.plot(df['Ugt'], df['epsilonp'], label='solid')
        plt.xlabel(r"superficial gas velocity [ m/s ]")
        plt.ylabel(r"overall holdups [ - ]")
        plt.legend()
        
        plt.subplot(2,2,2)
        plt.plot(df['Ugt'], df['egfb'])
        plt.xlabel(r"superficial gas velocity [ m/s ]")
        plt.ylabel(r"freeboard gas holdup [ - ]")
        
        plt.subplot(2,2,3)
        plt.plot(df['Ugt'], df['hbed'])
        plt.xlabel(r"superficial gas velocity [ m/s ]")
        plt.ylabel(r"bed height [ m ]")
        
        plt.subplot(2,2,4)
        plt.plot(df['Ugt'], df['Umf'])
        plt.xlabel(r"superficial gas velocity [ m/s ]")
        plt.ylabel(r"minimum fluidization velocity [ m/s ]")
        plt.savefig(fname + ".png", dpi=300)
        plt.close()
        
    def plotVariableRecycleEbullated(self, df, fname):
        """
        plots summary results for variable recycle ebullated bed
        
        args:
            df: data frame containing summary results
            fname: file name without extension
        """
        plt.figure(figsize = [14,15])
        plt.subplot(3,2,1)
        plt.plot(df['Ugt'], df['egbed'], label=r'gas')
        plt.plot(df['Ugt'], df['epsilonl'], label=r'liquid')
        plt.plot(df['Ugt'], df['epsilonp'], label=r'solids')
        plt.xlabel(r'inlet gas velocity [ m/s ]')
        plt.ylabel(r'bed holdup [ - ]')
        plt.legend()
        
        plt.subplot(3,2,2)
        plt.plot(df['Ugt'], df['egfb'])
        plt.xlabel(r'inlet gas velocity [ m/s ]')
        plt.ylabel(r'freeboard gas holdup')
        
        plt.subplot(3,2,3)
        plt.plot(df['Ugt'], df['eta'])
        plt.xlabel(r'inlet gas velocity [ m/s ]')
        plt.ylabel(r'separation efficiency')
        
        plt.subplot(3,2,4)
        plt.plot(df['Ugt'], df['qgr'], label=r'gas')
        plt.plot(df['Ugt'], df['qlr'], label=r'liquid')
        plt.xlabel(r'inlet gas velocity [ m/s ]')
        plt.ylabel(r'recycle flow rates [ $m^3/s$ ]')
        plt.legend()
        
        plt.subplot(3,2,5)
        plt.plot(df['Ugt'], df['Rgas'], label=r'gas')
        plt.plot(df['Ugt'], df['R'], label=r'liquid')
        plt.xlabel(r'inlet gas velocity [ m/s ]')
        plt.ylabel(r'recycle fractions [ - ]')
        plt.legend()
        
        plt.subplot(3,2,6)
        plt.plot(df['Ugt'], df['Umf'])
        plt.xlabel(r'inlet gas velocity [ m/s ]')
        plt.ylabel(r'minimum fluidization velocity [ m/s ]')
        
        plt.savefig(fname + ".png", dpi=300)
        plt.close()
    def plot_bsvd(self, df, fname):
        """
        plots bubble size and slip velocity distribution
        
        args:
            df: data frame containing summary results
            fname: file name without extension
        """
        plt.figure(figsize=[14,5])
        plt.subplot(1,2,1)
        plt.plot(df['dbi'], df['pdf'])
        plt.xlabel(r"volume equivalent diameter [ mm ]")
        plt.ylabel(r"probability density")
        
        plt.subplot(1,2,2)
        plt.plot(df['usi'], df['pdf'])
        plt.xlabel(r"bubble slip velocity [ m/s ]")
        plt.ylabel(r"probability density")
        plt.savefig(fname + ".png", dpi=300)
        plt.close()