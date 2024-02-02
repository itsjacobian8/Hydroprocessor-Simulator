#ifndef BED_HPP
#define BED_HPP

#include "fluidization.hpp"

template<typename realType> class bed: public fluidization<realType>
{
private:
// private data
	// column diameter
	realType Dc;
	// recycle line diameter
	realType drec;	
	// mass of catalyst inventory
	realType mcat;

	// wall effect parameter
	realType k;
	// particle-liquid terminal velocity
	realType upinf;
	// Richardson-Zaki coefficient
	realType n;
	
// private member functions
	// computes wall effect parameter
	inline void wallEffect()
	{
		k = 1.0 - 1.33*(this->dvol/Dc);
	}
	// computes particle-liquid terminal velocity
	inline void terminalVelocity()
	{
		realType den1 = 18.0/pow(this->Arp, 2.0/3.0); 
		realType den2 = (2.335 - 1.744*this->phi)/pow(this->Arp, 1.0/6.0);
		realType Repinf = cbrt(this->Arp)/(den1 + den2);
		upinf = (Repinf*this->mul)/(this->rhol*this->dvol);
		
	}
	// computes Richardson-Zaki coefficient
	inline void RichardsonZaki()
	{
		realType A = 0.043*pow(this->Arp, 0.57)*(1.0 - 1.24*pow(this->dvol/Dc, 0.27));
		n = (4.8 + 2.4*A)/(1.0 + A);		
	}

protected:
// protected data
	// column's cross-sectional area
	realType Ac;
public:
// public data
	// particles holdup
	realType epsilonp;
	// bed height
	realType hbed;

// public member function
	// constructor
	bed
	(
		realType _rhog, 
		realType _rhol, 
		realType _mul, 
		realType _rhop, 
		realType _Vp, 
		realType _Ap, 
		realType _mcat, 
		realType _Dc, 
		realType _drec
	);
	
	// updates parameters
	virtual void updateParameters();

	// computes particles holdup
	void particlesHoldup(const realType& Ugbed, const realType& Ulbed);

	// computes bed height
	void bedHeight();
};

template<typename realType> bed<realType>::bed
(
	realType _rhog,
	realType _rhol,
	realType _mul,
	realType _rhop,
	realType _Vp,
	realType _Ap,
	realType _mcat,
	realType _Dc,
	realType _drec
)
:
fluidization<realType>::fluidization(_rhog, _rhol, _mul, _rhop, _Vp, _Ap),
mcat(_mcat),
Dc(_Dc),
drec(_drec)
{}

template<typename realType> void bed<realType>::updateParameters()
{
	fluidization<realType>::updateParameters();
	wallEffect();
	RichardsonZaki();
	terminalVelocity();	
	Ac = 0.25*this->pi*(pow(Dc, 2) - pow(drec, 2));
	
	colorizedHeading("\nCALCULATED BED PARAMETERS");
	colorizedStringOutput("wall effect parameter [-] = ");
	colorizedNumericOutput(k);
	colorizedStringOutput("Richardson-Zaki coefficient [-] = ");
	colorizedNumericOutput(n);
	colorizedStringOutput("terminal velocity [m/s] = ");
	colorizedNumericOutput(upinf);	
}

template<typename realType> void bed<realType>::particlesHoldup
(
	const realType& Ugbed, 
	const realType& Ulbed
)
{
	realType A = pow(Ulbed/(k*upinf), 1.0/n)*(1.0 + 0.22*pow(Ugbed/Ulbed, 0.92));
	epsilonp = 1.0 - A;

}

template<typename realType> void bed<realType>::bedHeight()
{
	hbed = mcat/(this->rhop*Ac*epsilonp);
}

#endif
