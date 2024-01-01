#!/usr/bin/env python3
from inputs import *
import ebullatedBed
import bubbleColumn
from scipy.constants import pi
import pandas as pd
from plots import *

if __name__ == '__main__':
    # specify system to be simulated (bubble column vs. ebullated bed)
    system = getSystem()
    
    if(system.startswith("bubble") and system.endswith("column")):
        # get recycle information
        recycleInfo = getRecycleInfo()
        
        # get input parameters
        inputs = getBubbleColumnParameters(recycleInfo)
        
        # number of gas flows to be simulated
        nGasFlows = int(round((inputs["qgMax"] - inputs["qgMin"])/inputs["qgIncrement"] + 1))
        
        summary_list = []
        for i in range(nGasFlows):
            # next gas flow rate
            qg = inputs['qgMin'] + float(i)*inputs['qgIncrement']
        
            # call to solver
            if(recycleInfo):
                summary, bsvd = bubbleColumn.solve(qg, inputs['qlIn'], inputs['rhog'], inputs['rhol'], inputs['mul'], inputs['sigma'], inputs['P'], inputs['Dc'], inputs['dor'], inputs['nor'], inputs['nr'], inputs['classWidth'], inputs['alpha'], inputs['absTol'], inputs['relTol'], recycleInfo, inputs['R'], inputs['vsep'], inputs['drec'])
            else:
                summary, bsvd = bubbleColumn.solve(qg, inputs['qlIn'], inputs['rhog'], inputs['rhol'], inputs['mul'], inputs['sigma'], inputs['P'], inputs['Dc'], inputs['dor'], inputs['nor'], inputs['nr'], inputs['classWidth'], inputs['alpha'], inputs['absTol'], inputs['relTol'], recycleInfo)
            summary_list.append(summary)
            
            # write bubble size distribution to file
            df = pd.DataFrame(bsvd)
            fname = f"bubbleColumn-bsvd-QL{inputs['qlIn']}-QG{qg}"
            df.to_csv(fname + ".txt", sep='\t')
            plot_bsvd(df, fname)
            
        # write summary results to file
        summary_df = pd.DataFrame(summary_list)
        fname = f"bubbleColumn-summary-QL-{inputs['qlIn']}"
        summary_df.to_csv(fname + ".txt", sep='\t')
        plot_summary(summary_df, system, recycleInfo, fname)
        
    elif(system.startswith("ebullated") and system.endswith("bed")):
        # get recycle information
        recycleInfo = getRecycleInfo()
        
        # get input parameters
        inputs = getEbullatedBedParameters(recycleInfo)
        
        # number of gas flows to be simulated
        nGasFlows = int(round((inputs['qgtMax'] - inputs['qgtMin'])/inputs['qgtIncrement'] + 1))
        
        summary_list =[]
        for i in range(nGasFlows):
            # next gas flow rate
            qgt = inputs['qgtMin'] + float(i)*inputs['qgtIncrement']
            
            # call to solver
            if(recycleInfo):
                summary, bsvd = ebullatedBed.solve(qgt, inputs['qlvr'], inputs['rhog'], inputs['rhol'], inputs['mul'], inputs['mcat'], inputs['Vp'], inputs['Ap'], inputs['rhop'], inputs['sigma'], inputs['P'], inputs['Dc'], inputs['hset'], inputs['dor'], inputs['nor'], inputs['nr'], inputs['classWidth'], inputs['alpha'], inputs['absTol'], inputs['relTol'], recycleInfo, inputs['drec'], inputs['vsep'], inputs['Rstep'], inputs['maxIter'])
            else:
                summary, bsvd = ebullatedBed.solve(qgt, inputs['qlvr'], inputs['rhog'], inputs['rhol'], inputs['mul'], inputs['mcat'], inputs['Vp'], inputs['Ap'], inputs['rhop'], inputs['sigma'], inputs['P'], inputs['Dc'], inputs['hset'], inputs['dor'], inputs['nor'], inputs['nr'], inputs['classWidth'], inputs['alpha'], inputs['absTol'], inputs['relTol'], recycleInfo)
                
            summary_list.append(summary)
            
            # write bubble size distribution to file
            df = pd.DataFrame(bsvd)
            fname = f"ebullatedBed-bsvd-QL{inputs['qlvr']}-QG{qgt}"
            df.to_csv(fname + ".txt", sep='\t')
            plot_bsvd(df, fname)
            
        # write summary results to file
        summary_df = pd.DataFrame(summary_list)
        fname = f"ebullated-bed-summary-QL-{inputs['qlvr']}"
        summary_df.to_csv(fname + ".txt", sep='\t')
        plot_summary(summary_df, system, recycleInfo, fname)
