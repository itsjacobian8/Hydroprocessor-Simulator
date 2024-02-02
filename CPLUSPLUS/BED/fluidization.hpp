#ifndef FLUIDIZATION_HPP
#define FLUIDIZATION_HPP

template<typename realType> class fluidization
{
private:
// private data
	// particle volume
	realType Vp;	
	// particle surface area
	realType Ap;
	// minimum void fraction
	realType emfo;

// private member functions
	inline void minimumVoid()
	{
		emfo = 0.415/cbrt(phi);
	}
protected:
// protected data	
	// gas density
	realType rhog;
	// liquid density
	realType rhol;
	// liquid viscosity
	realType mul;
	// particle density
	realType rhop;
	
	// volume equivalent diameter
	realType dvol;
	// particle sphericity
	realType phi;
	// Archimedes number
	realType Arp;
	// constants
	const realType g = 9.81;
	const realType pi = M_PI;
// protected member functions
	// updates volume equivalent diameter
	inline void volumeEquivalentDiameter()
	{
		dvol = cbrt((6.0*Vp)/pi);
	}
	// updates particle sphericity
	inline void sphericity()
	{
		phi = cbrt(pi*pow(6.0*Vp, 2))/Ap;
	}
	// updates particle-liquid Archimedes number
	inline void ArchimedesNumber()
	{
		Arp = (rhol*(rhop - rhol)*g*pow(dvol, 3))/pow(mul, 2);
	}
public:
// public data
	// gas free minimum fluidization velocity
	realType Umfo;	
	// gas-liquid-solid minimum fluidization velocity
	realType Umf;

// public member functions
	// constructor
	fluidization
	(
		realType _rhog, 
		realType _rhol,
		realType _mul,
		realType _rhop, 
		realType _Vp, 
		realType _Ap
	);
	
	// update protected data
	virtual void updateParameters();	
	
	// gas free minimum fluidization
	void LSMinimumFluidization();
	
	// gas-liquid-solid minimum fluidization velocity
	void GLSMinimumFluidization(const realType& Ugbed);
};

template<typename realType> fluidization<realType>::fluidization
(
	realType _rhog,
	realType _rhol,
	realType _mul,
	realType _rhop,
	realType _Vp,
	realType _Ap
)
:
rhog(_rhog),
rhol(_rhol),
mul(_mul),
rhop(_rhop),
Vp(_Vp),
Ap(_Ap)
{}

template<typename realType> void fluidization<realType>::updateParameters()
{
	volumeEquivalentDiameter();
	sphericity();
	ArchimedesNumber();
	minimumVoid();
	
	colorizedHeading("\nCALCULATED FLUIDIZATION PARAMETERS");
	colorizedStringOutput("volume equivalent diameter [m] = ");
	colorizedNumericOutput(dvol);
	colorizedStringOutput("particle sphericity [-] = ");
	colorizedNumericOutput(phi);
	colorizedStringOutput("Archimedes number [-] = ");
	colorizedNumericOutput(Arp);
	colorizedStringOutput("minimum void fraction [-] = ");
	colorizedNumericOutput(emfo);
}

template<typename realType> void fluidization<realType>::LSMinimumFluidization()
{
	// constants
	realType C1 = (150.0*(1.0 - emfo))/(3.5*phi);
	realType C2 = pow(emfo, 3)/1.75;	
	// minimum Reynolds number
	realType Remfo = sqrt(pow(C1, 2) + C2*Arp) - C1;
	// gas free minimum fluidization velocity
	Umfo = (Remfo*mul)/(rhol*dvol);
}

template<typename realType> void fluidization<realType>::GLSMinimumFluidization
(
	const realType& Ugbed
)
{
	// initialize with gas free value
	Umf = Umfo;
	// machine precision
	realType tol = std::numeric_limits<realType>::epsilon();
	
	realType Umf_prev, factor, emf, alphamf, rhom, Argls, C1, C2, Remf, relError;
	do
	{
		// assign current value to previous value
		Umf_prev = Umf;
		// adjust minimum void fraction
		factor = 1.0 - Umf_prev/Umfo;
		emf = emfo*(1.0 - 0.34*factor + 0.22*pow(factor, 2));
		// minimum gas holdup
		alphamf = (0.16*Ugbed)/(emf*(Ugbed + Umf_prev));
		// gas-liquid mixture density
		rhom = rhog*alphamf + rhol*(1.0 - alphamf);
		// gas-liquid-solid Archimedes number
		Argls = (rhol*(rhop - rhom)*g*pow(dvol, 3))/pow(mul, 2);
		// constants
		C1 = (150.0*(1.0 - emf))/(3.5*phi);
		C2 = (pow(emf, 3)*pow(1.0 - alphamf, 3))/1.75;
		// minimum Reynolds number
		Remf = sqrt(pow(C1, 2) + C2*Argls) - C1;
		// update minimum fluidization velocity
		Umf = (Remf*mul)/(rhol*dvol);
		// update relative error
		relError = fabs((Umf - Umf_prev)/Umf);
	}while(relError > tol);
}

#endif
