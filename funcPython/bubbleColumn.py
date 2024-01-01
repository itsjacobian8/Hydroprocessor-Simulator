#!/usr/bin/env python3
from freeboard import *
from time import ctime, time
from termcolor import colored

def solve(qgIn, qlIn, rhog, rhol, mul, sigma, P, Dc, dor, nor, nr, classWidth, alpha, absTol, relTol, recycleInfo, R=0.0, vsep=0.0, drec=0.0):
    """
    solver for bubble column with or without internal recycle
    
    args:
        qgIn: inlet gas flow rate (float)
        qlIn: inlet liquid flow rate (float)
        rhog: gas density (float)
        rhol: liquid density (float)
        mul: liquid viscosity (float)
        sigma: gas-liquid interfacial tension (float)
        P: operating pressure (float)
        Dc: column diameter (float)
        dor: outlet orifice diameter (float)
        nor: number of outlet orifices per riser (float)
        nr: number of risers (float)
        classWidth: bubble class width (float)
        alpha: confidence interval for bubble size distribution (float)
        absTol: absolute tolerance (float)
        relTol: relative tolerance (float)
        recycleInfo: boolean to indicate absence or presence of internal recycle (float)
        R: liquid recycle ratio (float)
        vsep: gas-liquid separator volume (float)
        drec: recycle line diameter (float)
    
    return:
        summary: summary results (dictionary)
        bsvd: bubble size distribution (dictionary)
    """
    print(colored(f"[{ctime(time())}] SIMULATING QG = {qgIn:.4f}","green"))
    p = P/0.1
    if(recycleInfo):
        # bubbble column with internal recycle
        # this is different from air-lift reactor
        
        # initialize recycled gas and liquid flow rates
        qgr = 0.0
        qlr = 0.0
        
        # initialize column superficial velocities
        UgIn = qgIn/(0.25*pi*(Dc**2.0))
        UlIn = qlIn/(0.25*pi*(Dc**2.0))
        
        Ug = UgIn
        Ul = UlIn
        
        relError = 1.0
        
        while(relError > relTol):
            # assign current relative error to previous value
            relError_prev = relError
            
            # assign current column velocities to previous values
            Ug_prev = Ug
            Ul_prev = Ul
            
            # update gas-liquid separation efficiency
            eta = GLSeparator(rhog, rhol, mul, sigma, Ug_prev, Ul_prev, dor, nor, nr, Dc, drec, R, vsep, p, classWidth, alpha, absTol, relTol)
        
            # update recycle flow rates
            qgr, qlr = recycleFlowRates(Ug_prev, Ul_prev, Dc, drec, R, eta)
            
            # update column flow rates
            Ug, Ul = bedFlowRates(qgIn, qlIn, Dc, drec, qgr, qlr)

            # update relative error
            relError = abs(Ul - Ul_prev)/Ul + abs(Ug - Ug_prev)/Ug
            
            # absolute error on relError to break from loop
            absError = abs(relError - relError_prev)
            if(absError < absTol):
                print(colored(f"[{ctime(time())}] SOLUTION CONVERGED. REL. ERROR = {relError}", "green"))
                break
            print(colored(f"[{ctime(time())}] relative error = {relError}","green"))
        
        print("\n")
        # solve for bubble size and slip velocity distribution and gas holdup
        eg, dbi, usi, pdf = freeboardGasHoldup(rhog, rhol, sigma, mul, p, Ug, Ul, dor, nor, nr, Dc, drec, absTol, relTol, alpha, classWidth)
        
        # recycled gas fraction
        Rgas = R*(1.0 - eta)
        summary = {"UgIn":UgIn, "UlIn":UlIn, "Ug":Ug, "Ul":Ul, "eg":eg, "eta":eta, "qgr":qgr, "qlr":qlr, 'R':R, 'Rgas':Rgas}
    else:
        # superficial velocities
        Ug = qgIn/(0.25*pi*(Dc**2.0))
        Ul = qlIn/(0.25*pi*(Dc**2.0))
        
        # solution to bubble column with no internal recycle
        eg, dbi, usi, pdf = freeboardGasHoldup(rhog, rhol, sigma, mul, p, Ug, Ul, dor, nor, nr, Dc, drec, absTol, relTol, alpha, classWidth)
        
        summary = {"Ug":Ug, "Ul":Ul, "eg":eg}
    
    # bubble size distribution
    bsvd = {"dbi":dbi, "usi":usi, "pdf":pdf}
    
    return summary, bsvd