#!/usr/bin/env python3
from fluidization import *
from freeboard import *
from time import time, ctime
from scipy.constants import pi
from termcolor import colored
import pandas as pd

def solve(qgt, qlvr, rhog, rhol, mul, mcat, Vp, Ap, rhop, sigma, P, Dc, hset, dor, nor, nr, classWidth, alpha, absTol, relTol, recycleInfo, drec=0.0, vsep=0.0, Rstep=0.0, maxIter=0):
    """
    Iterative solver for an ebullated bed reactor with internal recycle
    
    args:
        qgt: inlet gas flow rate (float)
        qlvr: inlet liquid flow rate (float)
        rhog: gas density (float)
        rhol: liquid density (float)
        mul: liquid viscosity (float)
        mcat: mass of catalyst inventory (float)
        Vp: particle volume (float)
        Ap: particle surface area (float)
        rhop: particle density (float)
        sigma: gas-liquid interfacial tension (float)
        P: operating pressure (float)
        Dc: column diameter (float)
        hset: desired bed height (float)
        dor: outlet orifice diameter (float)
        nor: number of outlet orifices per riser (int but received as float)
        nr: number of risers on the distributor grid (int but received as float)
        classWidth: bubble class width (float)
        alpha: confidence interval for generating bubble size distribution (float)
        absTol: absolute tolerance (float)
        relTol: relative tolerance (float)
        recycleInfo: switch for internal recycle (boolean)
        drec: recycle line diameter (float)
        vsep: separator volume (float)
        Rstep: liquid recycle ratio step (float)
        maxIter: maximum number of iterations (int)
    
    return:
        summary: scalar output (dictionary)
        bsvd: bubble size and velocity distribution array organized into a (dictionary)
    """
    print(colored(f"[{ctime(time())}] SIMULATING QGT = {qgt:.4f} m^3/s","green"))
    # dimensionless operating pressure
    p = P/0.1
    if(recycleInfo):
        # initialize recycle flow rates
        qgr = 0.0
        qlr = 0.0
    
        # initialize liquid recycle ratio
        R = 0.00
    
        # initialize loop counter
        loopCounter = 0
        
        while(loopCounter < maxIter):
            # update bed flow rates
            Ugbed, Ulbed = bedFlowRates(qgt, qlvr, Dc, drec, qgr, qlr)
        
            # calculate minimum fluidization velocity
            dv = volumeEquivalentDiameter(Vp)
            phi = sphericity(Vp, Ap)
            Arp = archimedesNumber(rhol, rhop, mul, dv)
            emfo = minimumVoid(phi)
            Umfo = LSMinimumFluidization(rhol, mul, dv, phi, Arp, emfo)
            Umf = GLSMinimumFluidization(rhog, rhol, rhop, mul, Ugbed, dv, phi, emfo, Umfo, relTol)
        
            # ensure bed liquid velocity is greater than minimum fluidization velocity
            while(Ulbed < Umf):
                # increase recycle ratio since bed liquid velocity is still short of threshold
                R = max(min(R + Rstep, 1.0), 0.0)
            
                # update gas-liquid separation efficiency
                eta = GLSeparator(rhog, rhol, mul, sigma, Ugbed, Ulbed, dor, nor, nr, Dc, drec, R, vsep, p, classWidth, alpha, absTol, relTol)
            
                # update recycle flow rates
                qgr, qlr = recycleFlowRates(Ugbed, Ulbed, Dc, drec, R, eta)
            
                # update bed flow rates
                Ugbed, Ulbed = bedFlowRates(qgt, qlvr, Dc, drec, qgr, qlr)
            
                # update minimum fluidization velocity
                Umf = GLSMinimumFluidization(rhog, rhol, rhop, mul, Ugbed, dv, phi, emfo, Umfo, relTol)
        
            # calculate particles holdup
            k = wallEffectParameter(dv, Dc)
            n = richardsonZakiCoefficient(dv, Dc, Arp)
            upinf = terminalVelocity(rhol, mul, dv, phi, Arp)
            epsilonp = particlesHoldup(Ugbed, Ulbed, k, upinf, n)
        
            # calculate bed height
            hbed = bedHeight(rhop, drec, Dc, mcat, epsilonp)
        
            # deviation of calculated bed height from desired value
            error = (hbed - hset)
            
            # break from loop if the absolute error on desired bed height is smaller than specified tolerance
            if(abs(error) < absTol):
                print(colored(f"[{ctime(time())}] SOLUTION CONVERGED! HEIGHT = {hbed:.5f}, ERROR = {error:.5f}","green"))
                break
        
            # otherwise update liquid recycle ratio appropriately
            R = max(min(R - Rstep*error, 1.0), 0.0)
        
            # update gas-liquid separation efficiency
            eta = GLSeparator(rhog, rhol, mul, sigma, Ugbed, Ulbed, dor, nor, nr, Dc, drec, R, vsep, p, classWidth, alpha, absTol, relTol)
        
            # update recycle flow rates
            qgr, qlr = recycleFlowRates(Ugbed, Ulbed, Dc, drec, R, eta)
        
            # display simulation progress
            print(colored(f"[{ctime(time())}] iteration = {loopCounter}, height = {hbed:.5f}, error = {error:.5f}", "green"))
        
            # update loop counter
            loopCounter += 1
    
        # calculate freeboard gas holdup
        egfb, dbi, usi, pdf = freeboardGasHoldup(rhog, rhol, sigma, mul, p, Ugbed, Ulbed, dor, nor, nr, Dc, drec, absTol, relTol, alpha, classWidth)
    
        # bed gas holdup
        egbed = egfb/1.3
    
        # bed liquid holdup
        epsilonl = 1.0 - egbed - epsilonp
    
        # recycle gas fraction
        Rgas = R * (1.0 - eta)
    
        # inlet gas velocity
        Ugt = qgt/(0.25*pi*(Dc**2.0 - drec**2.0))
    
        # inlet liquid velocity
        Ulvr = qlvr/(0.25*pi*(Dc**2.0 - drec**2.0))
    
        # collect scalar results in a dictionary
        summary = {"Ugt":Ugt, "Ulvr":Ulvr, "Umf":Umf, "Ugbed":Ugbed, "Ulbed":Ulbed, "egbed":egbed, "epsilonl":epsilonl, "epsilonp":epsilonp, "egfb":egfb, "hbed":hbed, "eta":eta, "qgr":qgr, "qlr":qlr, "Rgas":Rgas, "R":R}

    else:
        # superficial velocities
        Ug = qgt/(0.25*pi*Dc**2.0)
        Ul = qlvr/(0.25*pi*Dc**2.0)
        
        # minimum fluidization velocity
        dv = volumeEquivalentDiameter(Vp)
        phi = sphericity(Vp, Ap)
        Arp = archimedesNumber(rhol, rhop, mul, dv)
        emfo = minimumVoid(phi)
        Umfo = LSMinimumFluidization(rhol, mul, dv, phi, Arp, emfo)
        Umf = GLSMinimumFluidization(rhog, rhol, rhop, mul, Ug, dv, phi, emfo, Umfo, relTol)
        
        # check if inlet flow rates are sufficient to fluidize the bed
        if(Ul < Umf):
            print(colored(f"[ctime(time())] Supplied inlet liquid flow rate insufficient to fluidize the bed", "red"))
            print(colored(f"[{ctime(time())}] Setting liquid velocity to a value 15% higher than minimum fluidization velocity", "green"))
            Ul = 1.15*Umf
            
        # solids holdup
        k = wallEffectParameter(dv, Dc)
        n = richardsonZakiCoefficient(dv, Dc, Arp)
        upinf = terminalVelocity(rhol, mul, dv, phi, Arp)
        epsilonp = particlesHoldup(Ug, Ul, k, upinf, n)
        
        # bed height
        hbed = bedHeight(rhop, drec, Dc, mcat, epsilonp)
        
        # check if bed height exceeds maximum column height
        if(hbed > hset):
            raise Exception("Bed height exceeds maximum bed height. OVERFLOW!")
            
        # determine freeboard gas holdups
        eg, dbi, usi, pdf = freeboardGasHoldup(rhog, rhol, sigma, mul, p, Ug, Ul, dor, nor, nr, Dc, drec, absTol, relTol, alpha, classWidth)
        
        # bed gas holdup
        egbed = eg/1.30
        
        # bed liquid holdup
        epsilonl = 1.0 - epsilonp - egbed
        
        # summary scalar results
        summary = {"Ug":Ug, "Ul":Ul, "Umf":Umf, "egbed":egbed, "epsilonl":epsilonl, "epsilonp":epsilonp, "egfb":eg, "hbed":hbed}
    
    # collect bubble size and velocity distribution data in a separate dictionary
    bsvd = {"dbi":dbi, "usi":usi, "pdf":pdf}
    
    print("\n")
    return summary, bsvd
        