#!/usr/bin/env python3
from scipy.constants import pi, g

class Fluidization:
    """
    Collection of methods for simulating the bed region of an ebullated bed reactor
    
    class attributes:
    -----------------
    inputs: dictionary input that includes the following input parameters (dict)
        Vp: particle volume (float)
        Ap: particle surface area (float)
        rhop: particle density (float)
        mcat: mass of catalyst inventory (float)
        rhol: liquid density (float)
        mul: liquid viscosity (float)
        Dc: column diameter (float)
        drec: recycle line diameter (float)
        relTol: relative tolerance (float)
    
    class methods:
    -------------
    volumeEquivalentDiameter: float
        returns volume equivalent diameter, dv
    
    sphericity: float
        returns particle sphericity, phi
    
    archimedesNumber: float
        returns particle-liquid Archimedes number, Arp
    
    minimumVoid: float
        returns void fraction at minimum fluidization, emfo
    
    LSMinimumFluidization: float
        returns gas free minimum fluidization velocity, Umfo
    
    GLSMinimumFluidization: float
        returns minimum fluidization velocity, Umfo 
    
    wallEffects: float
        returns wall effects parameter, k
    
    terminalVelocity: float
        returns particle terminal velocity, upinf
    
    RichardsonZaki: float
        returns Richardson-Zaki coefficient, n 
    
    particlesHoldup: float
        returns particles holdup, epsilonp
    
    bedHeight: float
        returns bed height, hbed
    
    for more information, consult docstrings for individual methods
    """
    def __init__(self, inputs: dict):
        self.Vp = inputs['Vp']
        self.Ap = inputs['Ap']
        self.rhop = inputs['rhop'] 
        self.mcat = inputs['mcat']
        self.rhol = inputs['rhol'] 
        self.mul = inputs['mul']
        self.Dc = inputs['Dc']
        self.drec = inputs['drec']
        self.relTol = inputs['relTol']
        self.Ac = 0.25*pi*(inputs['Dc']**2.0 - inputs['drec']**2.0)
        
    def bedHeight(self, epsilonp):
        """
        computes bed height
        
        args:
            epsilonp: particles holdup (float)
        
        return:
            hbed: bed height (float)
        """
        hbed = self.mcat/(self.rhop * self.Ac * epsilonp)
        return hbed
        
    def particlesHoldup(self, Ugbed, Ulbed, k, upinf, n):
        """
        computes particles holdup
        
        args:
            Ugbed: superficial bed gas velocity (float)
            Ulbed: superficial bed liquid velocity (float)
            k: wall effect parameter (float)
            upinf: particle terminal velocity (float)
            n: Richardson-Zaki coefficient (float)
        
        return:
            epsilonp: particles holdup (float)
        """
        epsilonp = 1.0 - (Ulbed/(k*upinf))**(1.0/n)*(1.0 + 0.22*(Ugbed/Ulbed)**0.92)
        return epsilonp
        
    def wallEffects(self, dv):
        """
        computes wall effect parameter
        
        args:
            dv: volume equivalent diameter (float)
            Dc: column diameter (float)
        
        return:
            k: wall effect parameter (float)
        """
        k = 1.0 - 1.33*(dv/self.Dc)
        return k
        
    def terminalVelocity(self, dv, phi, Arp):
        """
        computes particle terminal velocity
        
        args:
            dv: volume equivalent diameter (float)
            phi: particle sphericity (float)
            Arp: particle-liquid Archimedes numnber
        
        return:
            upinf: particle terminal velocity (float)
        """
        # Reynolds number at terminal conditions
        Repinf = (Arp**(1.0/3.0))/(18.0/(Arp**(2.0/3.0)) + (2.335 - 1.744*phi)/(Arp**(1.0/6.0)))
        
        # terminal velocity
        upinf = (self.mul*Repinf)/(self.rhol*dv)
        
        return upinf
        
    def RichardsonZaki(self, dv, Arp):
        """
        computes Richardson-Zaki coefficient
        
        args:
            dv: volume equivalent diameter (float)
            Arp: particle-liquid Archimedes number (float)
        
        return:
            n: Richardson-Zaki coefficient (float)
        """
        A = 0.043*(Arp**0.57)*(1.0 - 1.24*(dv/self.Dc)**0.27)
        n = (4.8 + 2.4*A)/(1.0 + A)
        return n
        
    def volumeEquivalentDiameter(self):
        """
        computes volume equivalent particle diameter
        
        return:
            dv: volume equivalent diameter (float)
        """
        dv = ((6.0*self.Vp)/pi)**(1.0/3.0)
        return dv
        
    def sphericity(self):
        """
        computes particle sphericity
        
        return:
            phi: particle sphericity (float)
        """
        phi = ((((6.0*self.Vp)**2.0)*pi)**(1.0/3.0))/self.Ap
        return phi
        
    def archimedesNumber(self, dv):
        """
        computes Archimedes number
        
        args:
            dv: volume equivalent diameter (float)
        
        return:
            Arp: particle-liquid Archimedes number (float)
        """
        Arp = (self.rhol*(self.rhop - self.rhol)*g*(dv**3.0))/(self.mul**2.0)
        return Arp
    
    def GLSMinimumFluidization(self, Ugbed, dv, phi, emfo, Umfo):
        """
        Iteratively computes gas-liquid-solid minimum fluidization velocity
        
        args:
            Ugbed: superficial bed gas velocity (float)
            dv: volume equivalent particle diameter (float)
            phi: particle sphericity (float)
            emfo: minimum void fraction (float)
            Umfo: gas free minimum fluidization velocity (float)
        
        return:
            Umf: minimum fluidization velocity (float)
        """
        # initialize minimum fluidization velocity with gas free value
        Umf = Umfo
        
        # initialize relative error
        relError = 1.0
        
        while(relError > self.relTol):
            # assign current minimum fluidization velocity to previous value
            Umf_prev = Umf
            
            # adjust minimum void fraction
            emf = emfo*(1.0 - 0.34*(1.0 - Umf_prev/Umfo) + 0.22*(1.0 - Umf_prev/Umfo)**2.0)
            
            # compute solids free gas holdup at minimum fluidization
            alphamf = (0.16*Ugbed)/(emf*(Ugbed + Umf_prev))
            
            # gas-liquid mixture density at minimum fluidization 
            rhom = self.rhog*alphamf + self.rhol*(1.0 - alphamf)
            
            # gas-liquid-solid Archimedes number
            Argls = (self.rhol*(self.rhop - rhom)*g*(dv**3.0))/(self.mul**2.0)
            
            # constants
            C1 = (150.0*(1.0 - emf))/(3.5*phi)
            C2 = ((emf**3.0)*(1.0 - alphamf)**3.0)/1.75
            
            # Reynolds number at minimum fluidization velocity
            Remf = (C1**2.0 + C2*Argls)**0.5 - C1
            
            # minimum fluidization velocity
            Umf = (self.mul*Remf)/(self.rhol*dv)
            
            # update relative error
            relError = abs(Umf - Umf_prev)/Umf
            
        return Umf
    
    def LSMinimumFluidization(self, dv, phi, Arp, emfo):
        """
        computes two phase liquid-solid minimum fluidization velocity 
        
        args:
            dv: volume equivalent particle diameter (float)
            phi: particle sphericity (float)
            Arp: particle-liquid Archimedes number (float)
            emfo: minimum void fraction (float)
        
        return:
            Umfo: minimum fluidization velocity (float)
        """
        # constants
        C1 = (150.0*(1.0 - emfo))/(3.50*phi)
        C2 = (emfo**3.0)/1.75
        
        # Reynolds number at minimum fluidization
        Remfo = (C1**2.0 + C2*Arp)**0.5 - C1
        
        # minimum fluidization velocity
        Umfo = (self.mul*Remfo)/(self.rhol*dv)
        
        return Umfo
        
    def minimumVoid(self, phi):
        """
        computes void fraction at minimum fluidization
        
        args:
            phi: particle sphericity (float)
        
        return:
            emfo: minimum void fraction (float)
        """
        emfo = 0.415/(phi**(1.0/3.0))
        return emfo
    
        
   
        
    
        
    
        
    