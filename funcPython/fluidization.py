#!/usr/bin/env python3
from scipy.constants import g, pi

def volumeEquivalentDiameter(Vp):
    """
    returns volume equivalent particle diameter

    args:
        Vp: particle volume (float)
    
    return:
        dv: volume equivalent diameter (float)
    """
    dv = ((6.0*Vp)/pi)**(1.0/3.0)
    
    return dv

def sphericity(Vp, Ap):
    """
    returns particle sphericity 

    args:
	Vp: particle volume (float)
	Ap: particle surface area (float)
    
    return:
         phi: particle sphericity (float)
    """
    phi = ((((6.0*Vp)**2.0) * pi)**(1.0/3.0))/Ap
    
    return phi

def archimedesNumber(rhol, rhop, mul, dv):
    """
    returns particle-liquid Archimedes number

    args:
        rhol: liquid density (float)
	    rhop: particle density (float)
	    mul: liquid viscosity (float)
	    dv: volume equivalent diameter (float)
    
    return:
        Arp: particle-liquid Archimedes number (float)
    """
    Arp = (rhol * (rhop - rhol) * g * (dv**3.0))/(mul**2.0)
    return Arp

def minimumVoid(phi):
    """
    returns void fraction at minimum fluidization (liquid-solid)

    args:
	phi: particle sphericity (float)
    
    return:
        emfo: minimum void fraction (float)
    """
    emfo = 0.415/(phi**(1.0/3.0))
    return emfo

def LSMinimumFluidization(rhol, mul, dv, phi, Arp, emfo):
    """
    returns minimum fluidization velocity (liquid-solid system)

    args:
	rhol: liquid density (float)
	mul: liquid viscosity (float)
	dv: volume equivalent particle diameter (float)
	phi: particle sphericity (float)
	Arp: particle-liquid Archimedes number (float)
	emfo: void fraction at minimum fluidization (float)
    
     return:
         Umfo: gas free minimum fluidization velocity (float)
    """
    # constants
    C1 = (150.0 * (1.0 - emfo))/(3.50 * phi)
    C2 = (emfo**3.0)/1.75

    # Reynolds number at minimum fluidization
    Remfo = (C1**2.0 + C2*Arp)**(1.0/2.0) - C1
    
    # minimum fluidization velocity
    Umfo = (mul*Remfo)/(rhol*dv)
    return Umfo

def GLSMinimumFluidization(rhog, rhol, rhop, mul, Ugbed, dv, phi, emfo, Umfo, relTol):
    """
    returns minimum fluidization velocity (gas-liquid-solid system)

    args:
	rhog: gas density (float)
	rhol: liquid density (float)
	rhop: particle density (float)
	mul: liquid viscosity (float)
	Ugbed: superficial gas velocity through the bed (float)
	dv: volume equivalent particle diameter (float)
	phi: particle sphericity (float)
	emfo: gas free minimum void fraction (float)
	Umfo: gas free minimum fluidization velocity (float)
	relTol: relative tolerance (float)
    
    return:
        Umf: minimum fluidization velocity (float)
    """
    # initialize minimum fluidization velocity with gas free value
    Umf = Umfo 

    # relative error
    relError = 1.0

    while(relError > relTol):
        # minimum fluidization at previous iteration
        Umf_prev = Umf

        # minimum void fraction
        emf = emfo * (1.0 - 0.34*(1.0 - (Umf_prev/Umfo)) + 0.22*(1.0 - (Umf_prev/Umfo))**2.0)

        # solids free gas holdup at minimum fluidization
        alphamf = (0.16*Ugbed)/(emf*(Ugbed + Umf_prev))

        # gas-liquid mixture density
        rhom = rhog*alphamf + rhol*(1.0 - alphamf)

        # gas-liquid-solid Archimedes number
        Argls = (rhol*(rhop - rhom)*g*(dv**3.0))/(mul**2.0)

        # constants
        C1 = (150.0*(1.0 - emf))/(3.5*phi)
        C2 = ((emf**3.0)*(1.0 - alphamf)**3.0)/1.75 

        # Reynolds number at minimum fluidization
        Remf = (C1**2.0 + C2*Argls)**(1.0/2.0) - C1
        
        # minimum fluidization velocity. 
        Umf = (mul*Remf)/(rhol*dv)

        # update relative error
        relError = abs(Umf - Umf_prev)/Umf_prev

    return Umf

def wallEffectParameter(dv, Dc):
    """
    returns wall effect parameter, k

    args:
	    dv: volume equivalent particle diameter (float)
	    Dc: column diameter (float)
    
    return:
        k: wall effect parameter (float)
    """
    k = 1.0 - 1.33*(dv/Dc)
    return k

def richardsonZakiCoefficient(dv, Dc, Arp):
    """
    returns Richardson-Zaki Coefficient

    args:
        dv: volume equivalent particle diameter (float)
	    Dc: column diameter (float)
	    Arp: particle-liquid Archimedes number (float)
    
    return:
        n: Richardson-Zaki coefficient (float)
    """
    A = 0.043*(Arp**0.57)*(1.0 - 1.24*(dv/Dc)**0.27)
    n = (4.8 + 2.4*A)/(1.0 + A)
    return n

def terminalVelocity(rhol, mul, dv, phi, Arp):
    """
    returns particle-liquid terminal velocity

    args:
        rhol: liquid density (float)
	    mul: liquid viscosity (float)
	    dv: volume equivalent particle diameter (float)
	    phi: particle sphericity (float)
	    Arp: particle-liquid Archimedes number (float)
    
    return:
        upinf: particle terminal velocity (float)
    """
    # Reynolds number at terminal conditions
    Repinf = (Arp**(1.0/3.0))/((18.0/(Arp**(2.0/3.0))) + (2.335 - 1.744*phi)/(Arp**(1.0/6.0)))
    
    # Temrinal velocity
    upinf = (mul*Repinf)/(rhol*dv)
    return upinf

def particlesHoldup(Ugbed, Ulbed, k, upinf, n):
    """
    returns particles holdup

    args:
        Ugbed: superficial gas velocity through the bed (float)
	    Ulbed: superficial liquid velocity through the bed (float)
	    k: wall effect parameter (float)
	    upinf: particle-liquid terminal velocity (float)
	    n: Richardson-Zaki coefficient (float)
    
    return:
        epsilonp: particles holdup (float)
    """
    epsilonp = 1.0 - (Ulbed/(k*upinf))**(1.0/n) * (1.0 + 0.22*(Ugbed/Ulbed)**0.92)
    return epsilonp

def bedHeight(rhop, drec, Dc, mcat, epsilonp):
    """
    returns bed height

    args:
        rhop: particle density (float)
	    drec: recycle line diameter (float)
	    Dc: column diameter (float)
	    mcat: mass of catalyst inventory (float)
	    epsilonp: particles holdup (float)
    
    return:
        hbed: bed height (float)
    """
    # bed cross-sectional area
    Ac = 0.25*pi*(Dc**2.0 - drec**2.0)
    
    # bed height
    hbed = mcat/(rhop*Ac*epsilonp)
    return hbed
