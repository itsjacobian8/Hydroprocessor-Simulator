#!/usr/bin/env python3
from time import ctime, time
from termcolor import colored
from freeboard import Freeboard

class Solvers(Freeboard):
    """
    collection of solvers for different two phase and three phase systems
    
    class inherits from Freeboard class.
    
    class attributes
    ------------------
    inputs: dictionary input that includes the following input parameters
        Rstep: liquid recycle ratio step size (float)
        maxIter: maximum number of iterations (int)
        hset: desired bed/column height (float)
    
    class methods:
    ----------------
    noRecycleBubbleColumn: dict, dict
        solver for a bubble column with no internal recycle
    
    constantRecycleBubbleColumn: dict, dict
        solver for a bubble column with constant liquid recycle ratio
    
    noRecycleEbullatedBed: dict, dict
        solver for ebullated bed with no internal recycle
    
    variableRecycleEbullatedBed: dict, dict
        solver for ebullated bed with variable liquid recycle ratio
    
    for more information, consult docstrings for individual methods.
    """
    def __init__(self, inputs: dict):
        super().__init__(inputs)
        self.Rstep = inputs['Rstep']
        self.maxIter = inputs['maxIter']
        self.hset = inputs['hset']
        
    def noRecycleBubbleColumn(self):
        """
        solver for bubble column with no internal recycle
        
        returns:
            summary: summary of results (dict)
            bsvd: bubble size and slip velocity distribution (dict)
        """
        Ugbed = self.qgt/self.Ac
        Ulbed = self.qlvr/self.Ac
        
        # There is only one function call for this use case
        egfb, dbi, usi, pdf = self.freeboardGasHoldup(Ugbed, Ulbed)
        
        # summary results
        summary = {"Ugt":Ugbed, "Ulvr":Ulbed, "egfb":egfb}
        
        # bubble size and velocity distribution
        bsvd = {"dbi":dbi, "usi":usi, "pdf":pdf}
        print("\n")
        return summary, bsvd
        
    def constantRecycleBubbleColumn(self, R):
        """
        solver for a bubble column with constant internal liquid recycle ratio
        
        args:
            R: constant liquid recycle ratio
        
        return:
            summary: summary results (dict)
            bsvd: bubble size and slip velocity distribution (dict)
        """
        # initialize recycled gas and liquid flow rates
        qgr = 0.0
        qlr = 0.0
        
        # initialize column superficial velocities
        Ugbed = (self.qgt + qgr)/self.Ac
        Ulbed = (self.qlvr + qlr)/self.Ac
        
        # initialize relative error
        relError = 1.0
        
        # loop counter
        count = int(0)
        
        while(relError > self.relTol):
            # assign current relative error to previous value
            relError_prev = relError
            
            # assign current column velocities to previous values
            Ugbed_prev = Ugbed
            Ulbed_prev = Ulbed
            
            # update gas-liquid separation efficiency
            eta = self.GLSeparator(Ugbed_prev, Ulbed_prev, R)
            
            # update recycle flow rates
            qgr, qlr = self.recycleFlowRates(Ugbed_prev, Ulbed_prev, R, eta)
            
            # update bed flow rates
            Ugbed, Ulbed = self.bedFlowRates(qgr, qlr)
            
            # update relative error
            relError = abs(Ugbed - Ugbed_prev)/Ugbed + abs(Ulbed - Ulbed_prev)/Ulbed
            
            # update absolute error
            absError = abs(relError - relError_prev)
            
            # break from loop if absolute error is below threshold and maximum number of iterations
            # has been exceeded
            if(count > self.maxIter and absError < self.absTol):
                print(colored(f"[{ctime(time())}] MAXIMUM NUMBER OF ITERATIONS EXCEEDED.", "green", "on_black"))
                break
            
            print(colored(f"[{ctime(time())}] iteration = {count}, relative error = {relError:.4f}, absolute error = {absError:.4f}", "green", "on_black"))
            count = count + 1
        print("\n")
        # solve for gas holdup, bubble size and slip velocity distribution
        egfb, dbi, usi, pdf = self.freeboardGasHoldup(Ugbed, Ulbed)
        
        # recycled gas fraction
        Rgas = R*(1.0 - eta)
        
        # inlet velocities
        Ugt = self.qgt/self.Ac
        Ulvr = self.qlvr/self.Ac
        
        summary = {"Ugt":Ugt, "Ulvr":Ulvr, "Ugbed":Ugbed, "Ulbed":Ulbed, "egfb":egfb, "eta":eta, "qgr":qgr, "qlr":qlr, "R":R, "Rgas":Rgas}
        bsvd = {"dbi":dbi, "usi":usi, "pdf":pdf}
        
        return summary, bsvd
        
    def noRecycleEbullatedBed(self):
        """
        solver for an ebullated bed reactor with no internal recycle
        
        return:
            summary: summary results (dict)
            bsvd: bubble size and slip velocity distribution (dict)
        """
        # superficial velocities
        Ugbed = self.qgt/self.Ac
        Ulbed = self.qlvr/self.Ac
        
        # compute minimum fluidization velocity
        dv = self.volumeEquivalentDiameter()
        phi = self.sphericity()
        Arp = self.archimedesNumber(dv)
        emfo = self.minimumVoid(phi)
        Umfo = self.LSMinimumFluidization(dv, phi, Arp, emfo)
        Umf = self.GLSMinimumFluidization(Ugbed, dv, phi, emfo, Umfo)
        
        # check if inlet flow rates are sufficient to fluidize the bed
        if(Ulbed < Umf):
            print(colored(f"[{ctime(time())}] Supplied inlet liquid flow rate insufficient to fluidize the bed", "red", "on_black"))
            print(colored(f"[{ctime(time())}] Setting liquid velocity to a value 15% higher than minimum fluidization velocity", "green", "on_black"))
            Ulbed = 1.15*Umf
        
        # solids holdup
        k = self.wallEffects(dv)
        n = self.RichardsonZaki(dv, Arp)
        upinf = self.terminalVelocity(dv, phi, Arp)
        epsilonp = self.particlesHoldup(Ugbed, Ulbed, k, upinf, n)
        
        # bed height
        hbed = self.bedHeight(epsilonp)
        
        # check if bed height exceeds maximum column height
        if(hbed > self.hset):
            raise Exception("MAXIMUM BED HEIGHT REACHED!")
        
        # determine freeboard gas holdup
        egfb, dbi, usi, pdf = self.freeboardGasHoldup(Ugbed, Ulbed)
        
        # bed gas holdup
        egbed = egfb/1.30
        
        # bed liquid holdup
        epsilonl = 1.0 - epsilonp - egbed
        
        # summary results
        summary = {"Ugt":Ugbed, "Ulvr":Ulbed, "Umf":Umf, "egbed":egbed, "epsilonl":epsilonl, "epsilonp":epsilonp, "egfb":egfb, "hbed":hbed}
        bsvd = {"dbi":dbi, "usi":usi, "pdf":pdf}
        
        print("\n")
        return summary, bsvd  
        
    def variableRecycleEbullatedBed(self):
        """
        solver for an ebullated bed reactor with variable internal recycle
        
        return:
            summary: summary results (dict)
            bsvd: bubble size and slip velocity distribution (dict)
        """
        # initialize recycle flow rates
        qgr = 0.0
        qlr = 0.0
        
        # initialize liquid recycle ratio
        R = 0.0
        
        # initialize loop counter
        count = int(0)
        
        while(count < self.maxIter):
            # update bed superficial velocities
            Ugbed, Ulbed = self.bedFlowRates(qgr, qlr)
            
            # calculate minimum fluidization velocity
            dv = self.volumeEquivalentDiameter()
            phi = self.sphericity()
            Arp = self.archimedesNumber(dv)
            emfo = self.minimumVoid(phi)
            Umfo = self.LSMinimumFluidization(dv, phi, Arp, emfo)
            Umf = self.GLSMinimumFluidization(Ugbed, dv, phi, emfo, Umfo)
            
            # ensure bed liquid velocity is greater than minimum fluidization velocity
            while(Ulbed < Umf):
                # increase recycle ratio since bed liquid velocity is still short of threshold
                R = max(min(R + self.Rstep, 1.0), 0.0)
                
                # update gas-liquid separation efficiency
                eta = self.GLSeparator(Ugbed, Ulbed, R)
                
                # update recycle flow rates
                qgr, qlr = self.recycleFlowRates(Ugbed, Ulbed, R, eta)
                
                # update bed superficial velocities
                Ugbed, Ulbed = self.bedFlowRates(qgr, qlr)
                
                # update minimum fluidization velocity
                Umf = self.GLSMinimumFluidization(Ugbed, dv, phi, emfo, Umfo)
                
            # calculate particles holdup
            k = self.wallEffects(dv)
            n = self.RichardsonZaki(dv, Arp)
            upinf = self.terminalVelocity(dv, phi, Arp)
            epsilonp = self.particlesHoldup(Ugbed, Ulbed, k, upinf, n)
            
            # calculate bed height
            hbed = self.bedHeight(epsilonp)
            
            # deviation of calculated bed height from desired value
            error = (hbed - self.hset)
            
            # break from loop if absolute error is smaller than tolerance
            if(abs(error) < self.absTol):
                print(colored(f"[{ctime(time())}] SOLUTION CONVERGED!", "green", "on_black"))
                break
            
            # otherwise update liquid recycle ratio appropriately
            R = max(min(R - self.Rstep*error, 1.0), 0.0)
            
            # update gas-liquid separation efficiency
            eta  = self.GLSeparator(Ugbed, Ulbed, eta)
            
            # update recycle flow rates
            qgr, qlr = self.recycleFlowRates(Ugbed, Ulbed, R, eta)
            
            # display simulation progress
            print(colored(f"[{ctime(time())}] iteration = {count}, height = {hbed:.5f}, error = {error:.5f}", "green", "on_black"))
            
            # update loop counter
            count = count + 1
        
        # calculate freeboard gas holdup
        egfb, dbi, usi, pdf = self.freeboardGasHoldup(Ugbed, Ulbed)
        
        # bed gas holdup
        egbed = egfb/1.30
        
        # bed liquid holdup
        epsilonl = 1.0 - egbed - epsilonp
        
        # recycle gas fraction
        Rgas = R*(1.0 - eta)
        
        # inlet velocities
        Ugt = self.qgt/self.Ac
        Ulvr = self.qlvr/self.Ac
        
        # summary results
        summary = {"Ugt":Ugt, "Ulvr":Ulvr, "Umf":Umf, "Ugbed":Ugbed, "Ulbed":Ulbed, "egbed":egbed, "epsilonl":epsilonl, 
        "epsilonp":epsilonp, "egfb":egfb, "hbed":hbed, "eta":eta, "qgr":qgr, "qlr":qlr, "Rgas":Rgas, "R":R}
        
        # bubble size and slip velocity distribution
        bsvd = {"dbi":dbi, "usi":usi, "pdf":pdf}
        
        print("\n")
        return summary, bsvd