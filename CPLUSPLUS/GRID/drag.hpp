#ifndef DRAG_HPP
#define DRAG_HPP

template<typename realType> class drag
{
public:
// DRAG MODELS
	// single bubble drag model for pure liquid systems
	// see Tomiyama et al. (1998) 
	realType pureTomiyama
	(
		const realType& Re, 
		const realType& Eo
	) const;
	
	// single bubble drag model for slightly contaminated liquid systems
	// see Tomiyama et al. (1998)
	realType slightlyTomiyama
	(
		const realType& Re, 
		const realType& Eo
	) const;
	
	// single bubble drag model for fully contaminated liquid systems
	// see Tomiyama et al. (1998)
	realType contaminatedTomiyama
	(
		const realType& Re, 
		const realType& Eo
	) const;
	
	// single bubble drag model for fully contaminated liquid systems
	// see Mach et al (2020)
	realType MachDrag
	(
		const realType& Re, 
		const realType& Eo,
		const realType& beta0
	) const;

// SWARM CORRECTION MODELS
	// swarm correction model proposed in Tomiyama et al. (1995)
	realType TomiyamaSwarm
	(
		const realType& egi,
		const realType& l
	) const;
	
	// swarm correction model proposed in Mach et al. (2020)
	// for fully contaminated polydisperse swarm
	realType MachSwarm
	(
		const realType& Eo, 
		const realType& egi,
		const realType& beta1,
		const realType& beta2,
		const realType& dimlessP
	) const;

	realType NoSwarmCorrection() const;
};
		
template<typename realType> realType drag<realType>::pureTomiyama
(
	const realType& Re, 
	const realType& Eo
) const
{	// Hadamard-Rybczynky drag coefficient 
	// for a single bubble rising under
	// creeping flow conditions in 
	// a pure liquid medium
	realType Cd0 = 16.0/Re;
	
	// correction for inertial effects
	realType fRe = 1.0 + 0.15*pow(Re, 0.687);

	// bounding inertial effects
	realType fReMax = 48.0/Re;

	// correction for buoyant and surface effects
	realType fEo = 8.0/3.0*Eo/(Eo + 4.0);

	// all effects combined
	realType Cdinf = fmax(fmin(Cd0*fRe, fReMax), fEo);
	return Cdinf;
}

template<typename realType> realType drag<realType>::slightlyTomiyama
(
	const realType& Re, 
	const realType& Eo
) const
{
	// drag coefficient of a rigid sphere
	// at creeping flow conditions.
	// contamination with surface active compound
	// was reported to render bubble surface rigid.
	realType Cd0 = 24.0/Re;
	
	// correction for inertial effects
	realType fRe = 1.0 + 0.15*pow(Re, 0.687);

	// bounding inertial effects
	realType fReMax = 72.0/Re;
	
	// correction for buoyant and surface effects
	realType fEo = 8.0/3.0*Eo/(Eo + 4.0);
	
	// all effects combined
	realType Cdinf  = fmax(fmin(Cd0*fRe, fReMax), fEo);
	
	return Cdinf;
}

template<typename realType> realType drag<realType>::contaminatedTomiyama
(
	const realType& Re, 
	const realType& Eo
) const
{
	// drag coefficient of a rigid sphere
	// at creeping flow conditions.
	realType Cd0 = 24.0/Re;
	
	// correction for inertial effects
	realType fRe = 1.0 + 0.15*pow(Re, 0.687);
	
	// correction for buoyant and surface effects
	realType fEo = (8.0*Eo)/(3.0*(Eo + 4.0));
	
	// all effects combined
	realType Cdinf = fmax(Cd0*fRe, fEo);
	return Cdinf;
}

template<typename realType> realType drag<realType>::MachDrag
(
	const realType& Re, 
	const realType& Eo,
	const realType& beta0
) const
{
	// Unlike Tomiyama models, here inertial 
	// and buoyant effects are accounted for with ReEo,
	// which also represents the ratio of form drag to viscous drag. 
	// The piece-wise function seen in Tomiyama models 
	// was abandoned to better match our data 
	realType Cdinf = 24.0/Re*(1.0 + pow(Re*Eo, beta0));
	return Cdinf;
}

template<typename realType> realType drag<realType>::TomiyamaSwarm
(
	const realType& egi,
	const realType& l
)const
{
	// I believe l is lift coefficient
	// should be set to zero if the drag coefficient
	// is not dependent on radial position.
	// Otherwise, this needs to be determined beforehand
	realType fswarm = (1.0 - egi)*pow(1.0 - egi, 3.0 - l);
	return fswarm;
}

template<typename realType> realType drag<realType>::MachSwarm
(
	const realType& Eo, 
	const realType& egi,
	const realType& beta1,
	const realType& beta2,
	const realType& dimlessP
) const
{
	// This swarm correction model
	// accounts for pressure and size/shape effects
	// via a non-dimensional ratio dimlessP/Eo,
	// where dimlessP is dimensionless pressure defined
	// relative to standard pressure and the non-dimensional 
	// ratio signifies the relative importance of pressure
	// and buoyant effects.
	realType m = beta1*(1.0 - exp(beta2*dimlessP/Eo));
	realType fswarm = pow(1.0 - egi, m);
	return fswarm;
}

template<typename realType> realType drag<realType>::NoSwarmCorrection() const
{
	realType fswarm = 1.0;
	return fswarm;
}
	
#endif
