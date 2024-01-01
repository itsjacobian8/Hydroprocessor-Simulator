#!/usr/bin/env python3
from fluidization import Fluidization
from scipy.stats import lognorm
from scipy.constants import g, pi
from numpy import log, linspace, exp, tanh, zeros

class Freeboard(Fluidization):
    """
    Collection of methods for computing various parameters associated with the freeboard zone
    
    class inherits attributes from the Fluidization class. 
    
    class attributes:
    ---------------------
    inputs: dictionary input that includes the following input parameters (dict)
        rhog: gas density (float)
        sigma: gas-liquid interfacial tension (float)
        absTol: absolute tolerance (float)
        dor: outlet orifice diameter (float)
        nor: number of outlet orifice per riser (int but received as float)
        nr: number of risers on the distributor grid (int but received as float)
        vsep: gas-liquid separator volume (float)
        p: dimensionless operating pressure (float)
        classwidth: bubble class width (float)
        alpha: confidence interval for generating bubble chord lengths (float)
        qgt: inlet gas flow rate (float)
        qlvr: inlet liquid flow rate (float)
    
    class methods:
    -----------------
    chordToDiameter: float
        converts bubble chord length to volume equivalent diameter
    
    MachDrag: float
        returns bubble drag coefficient
    
    slipVelocity: float
        returns bubble slip velocity
    
    objective: float
        objective function for iteratively computing gas holdup
    
    bisect: float
        custom bisection method for determining gas holdup
    
    lognormParams: float 
        returns lognormal distribution parameters
    
    GLSeparator: float
        returns gas-liquid separation efficiency
    
    recycleFlowRates: float, float
        updates recycled gas and liquid flow rates
    
    bedFlowRates: float, float
        updates bed superficial velocities
    
    freeboardGasHoldup: float, numpy array, numpy array, numpy array
        returns freeboard gas holdup and bubble size distribution
    
    for more information, consult docstrings for individual methods.
    """
    def __init__(self, inputs: dict):
        super().__init__(inputs)
        self.rhog = inputs['rhog']
        self.sigma = inputs['sigma']
        self.absTol = inputs['absTol']
        self.dor = inputs['dor']
        self.nor = inputs['nor']
        self.nr = inputs['nr']
        self.vsep = inputs['vsep']
        self.p = 10.0*inputs['P']
        self.classwidth = inputs['classWidth']
        self.alpha = inputs['alpha']
        self.qgt = inputs['qgt']
        self.qlvr = inputs['qlvr']
        self.Aor = 0.25*pi*inputs['dor']**2.0
    
    def freeboardGasHoldup(self, Ugbed, Ulbed):
        """
        computes freeboard gas holdup
        
        args:
            Ugbed: superficial gas velocity (float)
            Ulbed: superficial liquid velocity (float)
        
        return:
            egfb: freeboard gas holdup (float)
            dbi: volume equivalent bubble diameters (numpy array)
            usi: bubble slip velocities (numpy array)
            pdf: probability densities (numpy array)
        """
        # compute lognormal distribution parameters
        mu, s = self.lognormParams(Ugbed, Ulbed)
        
        # compute minimum and maximum bubble chord lengths
        cbiMin, cbiMax = lognorm.interval(self.alpha, s, loc=0.0, scale=exp(mu))
        
        # number of bubble classes
        nClasses = int(round((cbiMax - cbiMin)/self.classwidth + 1))
        
        # bubble chord lengths
        cbi = linspace(cbiMin, cbiMax, nClasses)
        
        # probability densities
        pdf = lognorm.pdf(cbi, s, loc=0.0, scale=exp(mu))
        
        # initialize volume equivalent bubble diameters
        dbi = zeros(nClasses)
        
        # initialize bubble slip velocities
        usi = zeros(nClasses)
        
        # running sums for computing freeboard gas holdup
        egnum = 0.0
        egden = 0.0
        
        for i in range(nClasses):
            # convert chord length to volume equivalent diameter
            dbi[i] = self.chordToDiameter(cbi[i])
            
            # compute bubble class holdup
            egi = self.bisect(dbi[i], Ugbed, Ulbed)
            
            # compute bubble slip velocity
            usi[i] = self.slipVelocity(dbi[i], egi)
            
            # accumulate numerator and denominator
            egnum = egnum + egi*pdf[i]
            egden = egden + pdf[i]
            
        # freeboard gas holdup
        egfb = egnum/egden
        
        return egfb, dbi, usi, pdf
    
    def bedFlowRates(self, qgr, qlr):
        """
        computes bed superficial velocities
        
        args:
            qgr: recycled gas flow rates (float)
            qlr: recycled liquid flow rates (float)
        
        return:
            Ugbed: superficial bed gas velocity (float)
            Ulbed: superficial liquid velocity (float)
        """
        Ugbed = (self.qgt + qgr)/self.Ac
        Ulbed = (self.qlvr + qlr)/self.Ac
        
        return Ugbed, Ulbed
    
    def recycleFlowRates(self, Ugbed, Ulbed, R, eta):
        """
        computes recycle gas and liquid flow rates
        
        args:
            Ugbed: superficial bed gas velocity (float)
            Ulbed: superficial bed liquid velocity (float)
            R: liquid recycle ratio (float)
            eta: gas-liquid separation efficiency (float)
        
        return:
            qgr: recycled gas flow rate (float)
            qlr: recycled liquid flow rate (float)
        """
        qgr = R*(1.0 - eta)*Ugbed*self.Ac
        qlr = R*Ulbed*self.Ac
        return qgr, qlr
    
    def GLSeparator(self, Ugbed, Ulbed, R):
        """
        computes gas-liquid separation efficiency of the recycle pan
        
        args:
            Ugbed: superficial gas velocity (float)
            Ulbed: superficial liquid velocity (float)
            R: liquid recycle ratio (float)
        
        return:
            eta: gas-liquid separation efficiency (float)
        """
        # compute lognormal distribution parameters
        mu, s = self.lognormParams(Ugbed, Ulbed)
        
        # compute minimum and maximum bubble chord length
        cbiMin, cbiMax = lognorm.interval(self.alpha, s, loc=0.0, scale=exp(mu))
        
        # determine number of bubble classes
        nClasses = int(round((cbiMax - cbiMin)/self.classwidth + 1))
        
        # generate bubble chord lengths
        cbi = linspace(cbiMin, cbiMax, nClasses)
        
        # generate probability densities
        pdf = lognorm.pdf(cbi, s, loc=0.0, scale=exp(mu))
        
        # running sums for computing separation efficiency
        sumpdfr = 0.0
        sumpdfb = 0.0
        
        for i in range(nClasses):
            # convert chord length to volume equivalent diameter
            dbi = self.chordToDiameter(cbi[i])
            
            # compute bubble class holdup
            egi = self.bisect(dbi, Ugbed, Ulbed)
            
            # compute bubble slip velocity
            usi = self.slipVelocity(dbi, egi)
            
            # compute bubble Reynolds number
            Reb = (self.rhol*usi*(0.001*dbi))/self.mul
            
            # compute bubble Eotvos number
            Eo = ((self.rhol - self.rhog)*g*(0.001*dbi)**2.0)/self.sigma
            
            # constant
            Bi = (Eo**0.5)/(-0.0037*Reb*tanh(dbi))
            
            # liquid residence time through the recycle pan
            kprime = self.vsep/(R*Ugbed*self.Ac)
            
            # individual gas-liquid separation efficiency
            etai = 0.29*log(kprime) + Bi
            
            # probability of entraining the bubble class in the liquid recycle
            pdfr = (1.0 - max(min(etai, 1.0), 0.0))*R*pdf[i]
            
            # accumulate bed and recycle probabilities
            sumpdfr = sumpdfr + pdfr
            sumpdfb = sumpdfb + pdf[i]
            
        # separation efficiency
        eta = 1.0 - sumpdfr/(R*sumpdfb)
        
        return eta
    
    def lognormParams(self, Ugbed, Ulbed):
        """
        computes lognormal distribution parameters
        
        the bubble size distribution generated at the grid is assumed
        to be preserved through the bed and the freeboard
        
        assumption is valid when bubble coalescence is significantly inhibited 
        and particles holdup in the bed is low (<25%)
        
        args:
            Ugbed: superficial gas velocity (float)
            Ulbed: superficial liquid velocity (float)
        
        return:
            mu: lognormal mean bubble size (float)
            s: lognormal standard deviation (float)
        """
        # superficial gas velocity through the outlet orifice
        Ugor = Ugbed/(self.nr*self.nor) * (self.Ac/self.Aor)
        
        # superficial liquid velocity through the outlet orifice
        Ulor = Ulbed/(self.nr*self.nor)*(self.Ac/self.Aor)
        
        # outlet orifice based Reynolds number
        Re = (self.rhol*Ulor*self.dor)/self.mul
        
        # outlet orifice based gas froude number
        Fr = Ugor/((g*self.dor)**0.5)
        
        # gas-liquid density ratio
        rhoRatio = self.rhog/self.rhol
        
        # mean bubble size
        muStar = -0.5305 - 0.1469*log(rhoRatio) - 0.4229*log(Re) + 0.345*log(Fr)
        
        # dimensionless mean bubble size
        mu = muStar + log(1000.0*self.dor)
        
        # standard deviation
        s = 0.7594 - 0.0426*log(rhoRatio) - 0.0201*log(Re) + 0.0596*log(Fr)
        
        return mu, s
        
    def chordToDiameter(self, cbi):
        """
        converts bubble chord length to volume equivalent diameter
        
        args:
            cbi: bubble chord length (float)
        
        return:
            dbi: volume equivalent diameter (float)
        """
        # Eotvos number
        Eo = ((self.rhol - self.rhog)*g*(0.0015*cbi)**2.0)/self.sigma
        
        # bubble aspect ratio
        E = 1.0/(1.0 + 0.163*(Eo**0.757))
        
        # initial estimate of volume equivalent diameter
        dbi = 1.5*cbi*E**(-2.0/3.0)
        
        # initialize relative error
        relError = 1.0
        
        while(relError > self.relTol):
            # assign current diameter to previous value
            dbi_prev = dbi
            
            # update Eotvos number
            Eo = ((self.rhol - self.rhog)*g*(0.001*dbi_prev)**2.0)/self.sigma
            
            # update bubble aspect ratio
            E = 1.0/(1.0 + 0.163*(Eo**0.757))
            
            # update bubble diameter
            dbi = 1.5*cbi*E**(-2.0/3.0)
            
            # update relative error
            relError = abs(dbi - dbi_prev)/dbi
        
        return dbi
    
    def bisect(self, dbi, Ugbed, Ulbed):
        """
        custom bisection method for iteratively computing the gas holdup
        
        args:
            dbi: volume equivalent bubble diameter (float)
            Ugbed: superficial gas velocity (float)
            Ulbed: superficial liquid velocity (float)
        
        return:
            egi: gas holdup (float)
        """
        # solution bracket -- gas holdup must be between 0 and 1
        egi_low = 0.0
        egi_high = 1.0
        
        while((egi_high - egi_low) > self.absTol):
            # determine midpoint
            egi = 0.5*(egi_high + egi_low)
            
            # evaluate the objective at the midpoint of the solution bracket
            midEval = self.objective(egi, dbi, Ugbed, Ulbed)
            
            # check if the root has been located
            if(midEval == 0.0):
                break
            
            # evaluate the objective at the lower end of the solution bracket
            lowEval = self.objective(egi_low, dbi, Ugbed, Ulbed)
            
            # adjust solution bracket
            if(lowEval*midEval < 0.0):
                # lower and midpoint evaluations have opposite signs
                # solution must be between the lower end and midpoint of the solution bracket
                # assign midpoint to the upper end of the solution bracket
                egi_high = egi
            elif(lowEval*midEval > 0.0):
                # lower and midpoint evaluations have opposite signs
                # solution must be between the midpoint and the upper end of the solution bracket
                # assign midpoint to lower end of the solution bracket
                egi_low = egi
        return egi
    
    def objective(self, egi, dbi, Ugbed, Ulbed):
        """
        objective function for iteratively determining the gas holdup
        by minimizing the absolute error on bubble slip velocity
        
        slip velocity computed using a drag model is compared to the one 
        computed from superficial velocities and individual holdups
        
        The latter is claimed to be applicable in the bubbly flow regime
        
        args:
            egi: initial guess of gas holdup (float)
            dbi: volume equivalent bubble diameter (float)
            Ugbed: superficial gas velocity (float)
            Ulbed: superficial liquid velocity (float)
        
        return:
            error: error on bubble slip velocity (float)
        """
        # determine bubble slip velocity
        usi = self.slipVelocity(dbi, egi)
        
        # absolute error on bubble slip velocity
        error = Ugbed/max(egi, 1e-8) - Ulbed/max(1.0 - egi, 1e-8) - usi
        
        return error
    
    def slipVelocity(self, dbi, egi):
        """
        computes bubble slip velocity
        
        args:
            dbi: volume equivalent diameter (float)
            egi: individual gas holdup (float)
        
        return:
            usi: bubble slip velocity (float)
        """
        # initialize bubble slip velocity
        usi = 0.001
        
        # initialize relative error
        relError = 1.0
        
        while(relError > self.relTol):
            # assign current slip velocity to previous value
            usi_prev = usi
            
            # compute bubble drag coefficient
            Cdi, Eo = self.MachDrag(dbi, egi, usi_prev)
            
            # compute bubble aspect ratio
            E = 1.0/(1.0 + 0.163*(Eo**0.757))
            
            # update bubble slip velocity
            usi = (4.0/3.0 * (g*0.001*dbi)/Cdi * E**(2.0/3.0) * (self.rhol - self.rhog)/self.rhol)**0.5
            
            # update relative error
            relError = abs(usi - usi_prev)/usi
            
        return usi
        
    def MachDrag(self, dbi, egi, usi):
        """
        computes bubble drag coefficient.
        The bubble is assumed to rise with other bubbles in a liquid
        medium contaminated with surface active compounds
        
        args:
            dbi: volume equivalent diameter (float)
            egi: individual gas holdup (float)
            usi: bubble slip velocity (float)
        
        return:
            Cdi: individual bubble drag coefficient (float)
            Eo: Eotvos number
        """
        # bubble Reynolds number
        Reb = (self.rhol*usi*(0.001*dbi))/self.mul
        
        # bubble Eotvos number
        Eo = ((self.rhol - self.rhog)*g*(0.001*dbi)**2.0)/self.sigma
        
        # bubble drag coefficient
        Cdi = 24.0/Reb * (1.0 + (Reb*Eo)**0.47)*(1.0 - egi)**(2.46*(1.0 - exp(-0.13*self.p/Eo)))
        
        return Cdi, Eo
        
    
        
    
        
    
        
    
        
    
        
    
        
    
        
    
        