#!/usr/bin/env python3
from scipy.stats import lognorm
from scipy.constants import g, pi
from numpy import log, linspace, exp, tanh, zeros

def chordToDiameter(cbi, rhog, rhol, sigma, relTol):
    """
    converts bubble chord length to volume equivalent diameter
    
    args:
        cbi: bubble chord length (float)
        rhog: gas density (float)
        rhol: liquid density (float)
        sigma: gas-liquid interfacial tension (float)
        relTol: relative tolerance (float)
    
    return:
        dbi: volume equivalent diameter (float)
    """
    # Eotvos number
    Eo = ((rhol - rhog)*g*((0.0015*cbi)**2.0))/sigma
    
    # bubble aspect ratio
    aspect_ratio = 1.0/(1.0 + 0.163*(Eo**0.757))
    
    # initial estimate of volume equivalent diameter
    dbi = 1.5*cbi*aspect_ratio**(-2.0/3.0)
    
    relError = 1.0
    
    while(relError > relTol):
        # assign current diameter to previous value
        dbi_prev = dbi
        
        # update Eotvos number
        Eo = ((rhol - rhog)*g*((0.001*dbi_prev)**2.0))/sigma
        
        # update bubble aspect ratio
        aspect_ratio = 1.0/(1.0 + 0.163*(Eo**0.757))
        
        # update bubble diameter
        dbi = 1.5*cbi*aspect_ratio**(-2.0/3.0)
        
        # compute relative error
        relError = abs(dbi - dbi_prev)/dbi
        
    return dbi
    
    
def MachDrag(dbi, egi, usi, rhog, rhol, sigma, mul, p):
    """
    returns bubble drag coefficient
    
    args:
        dbi: volume equivalent diameter (float)
        egi: gas holdup (float)
        usi: bubble slip velocity (float)
        rhog: gas density (float)
        rhol: liquid density (float)
        sigma: gas-liquid interfacial tension (float)
        mul: liquid viscosity (float)
        p: dimensionless operating pressure (float)
    
    return:
        Cdi: bubble drag coefficient (float)
    """
    # bubble Reynolds number
    Reb = (rhol*usi*(0.001*dbi))/mul
    
    # bubble Eotvos number
    Eo = ((rhol - rhog)*g*(0.001*dbi)**2.0)/sigma
    
    # bubble drag coefficient
    Cdi = 24.0/Reb * (1.0 + (Reb*Eo)**0.47)*(1.0-egi)**(2.46*(1.0 - exp(-0.13*p/Eo)))
    return Cdi
    
    
def slipVelocity(dbi, egi, rhog, rhol, sigma, mul, p, relTol):
    """
    returns bubble slip velocity
    
    args:
        dbi: volume equivalent bubble diameter (float)
        egi: gas holdup (float)
        rhog: gas density (float)
        rhol: liquid density (float)
        sigma: gas-liquid interfacial tension (float)
        mul: liquid velocity (float)
        p: dimensionless operating pressure (float)
        relTol: relative tolerance (float)
    
    return:
        usi: bubble slip velocity (float)
    """
    # initialize bubble slip velocity
    usi = 0.001
    
    relError = 1.0
    while(relError > relTol):
        # assign current slip velocity to previous value
        usi_prev = usi
        
        # compute bubble drag coefficient
        Cdi = MachDrag(dbi, egi, usi_prev, rhog, rhol, sigma, mul, p)
        
        # bubble Eotvos number
        Eo = ((rhol - rhog)*g*(0.001*dbi)**2.0)/sigma
        
        # bubble aspect ratio
        aspect_ratio = 1.0/(1.0 + 0.163*(Eo**0.757))
        
        # update slip velocity
        usi = (4.0/3.0 * (g*0.001*dbi)/Cdi * aspect_ratio**(2.0/3.0) * (rhol - rhog)/rhol)**0.5
        
        # update relative error
        relError = abs(usi - usi_prev)/usi
    
    return usi

def objectiveFunction(egi, dbi, rhog, rhol, sigma, mul, p, Ugbed, Ulbed, relTol):
    """
    objective function for iteratively determining the gas holdup by minimizing
    the absolute error on slip velocity
    
    args:
        egi: gas holdup (float)
        dbi: bubble diameter (float)
        rhog: gas density (float)
        rhol: liquid density (float)
        sigma: gas-liquid interfacial tension (float)
        mul: liquid viscosity (float)
        p: dimensionless operating pressure (float)
        Ugbed: superficial gas velocity (float)
        Ulbed: superficial liquid velocity (float)
        relTol: relative tolerance (float)
    
    return:
        error: error on bubble slip velocity (float)
    """
    # compute slip velocity
    usi = slipVelocity(dbi, egi, rhog, rhol, sigma, mul, p, relTol)
    
    error = Ugbed/max(egi, 1e-8) - Ulbed/max(1.0 - egi, 1e-8) - usi
    return error

def bisect(dbi, rhog, rhol, sigma, mul, p, Ugbed, Ulbed, absTol, relTol):
    """
    custom bisection method for iteratively computing the gas holdup
    
    args:
        dbi: volume equivalent bubble diameter (float)
        rhog: gas density (float)
        rhol: liquid density (float)
        sigma: gas-liquid interfacial tension (float)
        mul: liquid viscosity (float)
        p: dimensionless operating pressure (float)
        Ugbed: superificial gas velocity (float)
        Ulbed: superficial liquid velocity (float)
        absTol: absolute tolerance (float)
        relTol: relative tolerance (float)
    
    return:
        egi: gas holdup
    """
    # solution bracket -- gas holdup must be between 0 and 1
    egi_low = 0.0
    egi_high = 1.0
    
    while((egi_high - egi_low) > absTol):
        # mid point
        egi = 0.5*(egi_high + egi_low)
        
        # evaluate objective function at the midpoint
        midEval = objectiveFunction(egi, dbi, rhog, rhol, sigma, mul, p, Ugbed, Ulbed, relTol)
        
        # check if the root has been located
        if(midEval == 0.0):
            break
            
        # evaluate objective function at the lower end of the solution bracket
        lowEval = objectiveFunction(egi_low, dbi, rhog, rhol, sigma, mul, p, Ugbed, Ulbed, relTol)
        
        # adjust solution bracket
        if(lowEval*midEval < 0.0):
            # lower end and midpoint evaluations have opposite signs
            # solution must be between lower end and midpoint
            # assign midpoint to upper end
            egi_high = egi
        elif(lowEval*midEval > 0.0):
            # lower end and midpoint evaluations have the same sign
            # solution must be between the midpoint and upper end
            # assign midpoint to lower end of the solution bracket
            egi_low = egi
    
    return egi
    
def lognormDistributionParameters(rhog, rhol, mul, Ugbed, Ulbed, dor, nor, nr, drec, Dc):
	"""
	returns lognormal distribution parameters

	The bubble size distribution generated at the grid is assumed 
	to be preserved through the bed and the freeboard.

	This assumption is valid when bubble coalescence is significantly inhibited
	and particles holdup in the bed is low (<25%)

	args:
		rhog: gas density (float)
		rhol: liquid density (float)
		mul: liquid viscosity (float)
		Ugbed superficial gas velocity (float)
		Ulbed: superficial liquid velocity (float)
		dor: distributor outlet orifice diameter (float)
		nor: number of outlet orifices per riser (float)
		nr: number of risers on the grid (float)
		drec: recycle line diameter (float)
		Dc: column diameter (float)
    
    return:
        mu: lognormal mean bubble size (float)
        s: lognormal standard deviation (float)
	"""
	# cross-sectional area of outlet orifice
	Aor = 0.25*pi*dor**2.0

	# column's cross-sectional area
	Ac = 0.25*pi*(Dc**2.0 - drec**2.0)

	# superificial gas velocity through the outlet orifice
	Ugor = Ugbed/(nr*nor) *(Ac/Aor)

	# superficial liquid velocity through the outlet orifice
	Ulor = Ulbed/(nr*nor) * (Ac/Aor)

	# outlet orifice based Reynolds number
	Re = (rhol*Ulor*dor)/mul

	# outlet orifice based gas Froude number
	Fr = (Ugor)/((g*dor)**0.5)

	# gas-liquid density ratio
	rhoRatio = rhog/rhol

	# mean bubble size
	muStar = -0.5305 - 0.1469*log(rhoRatio) - 0.4229*log(Re) + 0.345*log(Fr)

	# dimensionless mean bubble size
	mu = muStar + log(1000.0*dor)

	# standard deviation
	s = 0.7594 - 0.0426*log(rhoRatio) - 0.0201*log(Re) + 0.0596*log(Fr)

	return mu, s

def GLSeparator(rhog, rhol, mul, sigma, Ugbed, Ulbed, dor, nor, nr, Dc, drec, R, vsep, p, classWidth, alpha, absTol, relTol):
    """
    returns gas-liquid separation efficiency of the recycle pan
    
    args:
        rhog: gas density (float)
        rhol: liquid density (float)
        mul: liquid viscosity (float)
        sigma: gas-liquid interfacial tension (float)
        Ugbed: superficial gas velocity (float)
        Ulbed: superficial liquid velocity (float)
        dor: outlet orifice diameter (float)
        nor: number of outlet orifices per riser (float)
        nr: number of risers on the grid (float)
        Dc: column diameter (float)
        drec: recycle line diameter (float)
        R: liquid recycle ratio (float)
        vsep: separator volume (float)
        p: dimensionless operating pressure (float)
        classWidth: bubble class width (float)
        alpha: confidence interval for generating bubble size distribution (float)
        absTol: absolute tolerance (float)
        relTol: relative tolerance (float)
    
    return:
        eta: gas-liquid separation efficiency (float)
    """
    # compute lognormal distribution parameters
    mu, s = lognormDistributionParameters(rhog, rhol, mul, Ugbed, Ulbed, dor, nor, nr, drec, Dc)
    
    # compute minimum and maximum bubble chord length
    cbiMin, cbiMax = lognorm.interval(alpha, s, loc = 0.0, scale = exp(mu))
    
    # number of bubble classes
    nClasses = int(round((cbiMax - cbiMin)/classWidth) + 1)
    
    # bubble chord length
    cbi = linspace(cbiMin, cbiMax, nClasses)
    
    # probability density of individual bubble classes
    pdf = lognorm.pdf(cbi, s, loc=0.0, scale = exp(mu))
    
    # column's cross-sectional area
    Ac = 0.25*pi*(Dc**2.0 - drec**2.0)
    
    sumpdfr = 0.0
    sumpdfb = 0.0
    
    for i in range(nClasses):
        # convert chord length to volume equivalent diameter
        dbi = chordToDiameter(cbi[i], rhog, rhol, sigma, relTol)
        
        # calculate bubble class holdup
        egi = bisect(dbi, rhog, rhol, sigma, mul, p, Ugbed, Ulbed, absTol, relTol)
        
        # calculate bubble slip velocity
        usi = slipVelocity(dbi, egi, rhog, rhol, sigma, mul, p, relTol)
        
        # bubble Reynolds number
        Reb = (rhol*usi*(0.001*dbi))/mul
        
        # bubble Eotvos number
        Eo = ((rhol - rhog)*g*(0.001*dbi)**2.0)/sigma
        
        # constant 
        Bi = (Eo**0.5)/(-0.0037*Reb*tanh(dbi))
        
        # liquid residence time through the recycle pan
        kprime = vsep/(R*Ulbed*Ac)
        
        # individual gas-liquid separation efficiency
        etai = 0.29*log(kprime) + Bi
        
        # probability of entraining a bubble class in the liquid recycle
        pdfr = (1.0 - max(min(etai, 1.0), 0.0))*R*pdf[i]
        
        # accumulate bed and recycle probabilities
        sumpdfr = sumpdfr + pdfr
        sumpdfb = sumpdfb + pdf[i]
    
    # separation efficiency
    eta = 1.0 - sumpdfr/(R*sumpdfb)
    
    return eta

def recycleFlowRates(Ugbed, Ulbed, Dc, drec, R, eta):
    """
    returns recycled gas and liquid flow rates
    
    args:
        Ugbed: superficial gas velocity (float)
        Ulbed: superficial liquid velocity (float)
        Dc: column diameter (float)
        drec: recycle diameter (float)
        R: liquid recycle ratio (float)
        eta: gas-liquid separation efficiency (float)
    
    return:
        qgr: recycle gas flow rate (float)
        qlr: recycle liquid flow rate (float)
    """
    # column diameter
    Ac = 0.25*pi*(Dc**2.0 - drec**2.0)
    
    # recycled gas flow rate
    qgr = R*(1.0 - eta)*Ugbed*Ac
    
    # recycled liquid flow rate
    qlr = R*Ulbed*Ac
    
    return qgr, qlr
    
def bedFlowRates(qgt, qlvr, Dc, drec, qgr, qlr):
    """
    returns bed flow rates
    
    args:
        qgt: inlet gas flow rate (float)
        qlvr: inlet liquid flow rate (float)
        Dc: Column Diameter (float)
        drec: recycle line diameter (float)
        qgr: recycled gas flow rate (float)
        qlr: recycled liquid flow rate (float)
    
    return:
        Ugbed: bed superficial gas velocity (float)
        Ulbed: bed superficial liquid velocity (float)
    """
    # column's cross-sectional area
    Ac = 0.25*pi*(Dc**2.0 - drec**2.0)
    
    # bed superficial gas velocity
    Ugbed = (qgt + qgr)/Ac
    
    # bed superficial liquid velocity
    Ulbed = (qlvr + qlr)/Ac
    
    return Ugbed, Ulbed
    
def freeboardGasHoldup(rhog, rhol, sigma, mul, p, Ugbed, Ulbed, dor, nor, nr, Dc, drec, absTol, relTol, alpha, classWidth):
    """
    returns freeboard gas holdup or bubble column gas holdup
    
    args:
        rhog: gas density (float)
        rhol: liquid density (float)
        sigma: gas-liquid interfacial tension (float)
        mul: liquid viscosity (float)
        p: dimensionless operating pressure (float)
        Ugbed: superficial gas velocity (float)
        Ulbed: superficial liquid velocity (float)
        dor: outlet orifice diameter (float)
        nor: number of outlet orifices per riser (float)
        nr: number of risers on the grid (float)
        Dc: column diameter (float)
        drec: recycle line diameter (float)
        absTol: absolute tolerance (float)
        relTol: relative tolerance (float)
        alpha: confidence interval for generating bubble size distribution (float)
        classWidth: bubble class width (float)
    
    return:
        egfb: freeboard gas holdup (float)
        dbi: volume equivalent bubble diameters (numpy array)
        usi: bubble slip velocities (numpy array)
        pdf: probability densities (numpy array)
    """
    # compute lognormal distribution parameters
    mu, s = lognormDistributionParameters(rhog, rhol, mul, Ugbed, Ulbed, dor, nor, nr, drec, Dc)
    
    # compute minimum and maximum bubble chord length
    cbiMin, cbiMax = lognorm.interval(alpha, s, loc = 0.0, scale = exp(mu))
    
    # number of bubble classes
    nClasses = int(round((cbiMax - cbiMin)/classWidth) + 1)
    
    # chord length distribution
    cbi = linspace(cbiMin, cbiMax, nClasses)
    
    # probability density of individual bubble classes
    pdf = lognorm.pdf(cbi, s, loc = 0.0, scale = exp(mu))
    
    # bubble sizes
    dbi = zeros(nClasses)
    
    # bubble slip velocities
    usi = zeros(nClasses)
    
    egnum = 0.0
    egden = 0.0
    for i in range(nClasses):
        # convert chord length to volume equivalent diameter
        dbi[i] = chordToDiameter(cbi[i], rhog, rhol, sigma, relTol)
        
        # calculate bubble class holdup
        egi = bisect(dbi[i], rhog, rhol, sigma, mul, p, Ugbed, Ulbed, absTol, relTol)
        
        # calculate bubble slip velocity
        usi[i] = slipVelocity(dbi[i], egi, rhog, rhol, sigma, mul, p, relTol)
        
        # accumulate numerator and denominator
        egnum = egnum + egi*pdf[i]
        egden = egden + pdf[i]
    
    # overall gas holdup
    egfb = egnum/egden
    
    return egfb, dbi, usi, pdf
    