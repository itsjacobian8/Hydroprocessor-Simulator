#ifndef fluidization_hpp_
#define fluidization_hpp_

#include<cmath>
#include<limits>
#include<iostream>

using namespace std;

template<typename realType> class fluidization
{
private:
// private data
	// particle volume
	realType Vp;	
	// particle surface area
	realType Ap;

//private member functions
	// computes volume equivalent diameter
	void volumeEquivalentDiameter();
	// computes particle sphericity
	void sphericity();
	// computes particle-liquid Archimedes number
	void ArchimedesNumber();
	// computes minimum void fraction
	void minimumVoid();

protected:
// protected data	
	// gas density
	realType rhog;
	// liquid density
	realType rhol;
	// particle density
	realType rhop;
	// liquid viscosity
	realType mul;

	// volume equivalent diameter
	realType dvol;
	// particle sphericity
	realType phi;
	// Archimedes number
	realType Arp;
	// minimum void fraction
	realType emfo;

	// constants
	const realType g = 9.81;
	const realType pi = M_PI;
public:
// public data
	// gas free minimum fluidization velocity
	realType Umfo;	
	// gas-liquid-solid minimum fluidization velocity
	realType Umf;

// public member functions
	// gets user inputs
	virtual void getInputs();
	// updates calculated parameters
	virtual void updateParameters();
	// gas free minimum fluidization
	void LSMinimumFluidization();
	// gas-liquid-solid minimum fluidization velocity
	void GLSMinimumFluidization(const realType& Ugbed);
};

template<typename realType> void fluidization<realType>::getInputs()
{
	cout<<"\nUSER INPUTS FOR FLUIDIZATION PARAMETERS"<<endl;
	cout<<"---------------------------------------"<<endl;
	
	cout<<"Enter gas density [kg/m^3]: ";
	cin>>rhog;
	cout<<"[RECEIVED] gas density = "<<rhog<<endl;

	cout<<"\nEnter liquid density [kg/m^3]: ";
	cin>>rhol;
	cout<<"[RECEIVED] liquid density = "<<rhol<<endl;
	
	cout<<"\nEnter liquid viscosity [Pa.s]: ";
	cin>>mul;
	cout<<"[RECEIVED] liquid viscosity = "<<mul<<endl;
	
	cout<<"\nEnter particle density [kg/m^3]: ";
	cin>>rhop;
	cout<<"[RECEIVED] particle density = "<<rhop<<endl;

	cout<<"\nEnter particle volume [m^3]: ";
	cin>>Vp;
	cout<<"[RECEIVED] particle volume = "<<Vp<<endl;

	cout<<"\nEnter particle surface area [m^2]: ";
	cin>>Ap;
	cout<<"[RECEIVED] particle surface area = "<<Ap<<endl;
}

template<typename realType> void fluidization<realType>::volumeEquivalentDiameter()
{
	dvol = cbrt((6.0*Vp)/pi);
}

template<typename realType> void fluidization<realType>::sphericity()
{
	phi = cbrt(pi*pow(6.0*Vp, 2))/Ap;
}

template<typename realType> void fluidization<realType>::ArchimedesNumber()
{
	Arp = (rhol*(rhop - rhol)*g*pow(dvol, 3))/pow(mul, 2);

}

template<typename realType> void fluidization<realType>::minimumVoid()
{
	emfo = 0.415/cbrt(phi);
}

template<typename realType> void fluidization<realType>::updateParameters()
{
	cout<<"\nCALCULATED FLUIDIZATION PARAMETERS."<<endl;
	cout<<"-----------------------------------"<<endl;
	
	volumeEquivalentDiameter();
	cout<<"volume equivalent diameter [m] = "<<dvol<<endl;

	sphericity();
	cout<<"particle sphericity [-] = "<<phi<<endl;

	ArchimedesNumber();
	cout<<"Archimedes number [-] = "<<Arp<<endl;

	minimumVoid();
	cout<<"minimum void fraction [-] = "<<emfo<<endl;
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

template<typename realType> void fluidization<realType>::GLSMinimumFluidization(const realType& Ugbed)
{
	// initialize with gas free value
	Umf = Umfo;
	
	realType Umf_prev, emf, alphamf, rhom, Argls, C1, C2, Remf, relError, relTol;
	relTol = numeric_limits<realType>::epsilon();
	do
	{
		// assign current value to previous value
		Umf_prev = Umf;

		// adjust minimum void fraction
		emf = emfo*(1.0 - 0.34*(1.0 - Umf_prev/Umfo) + 0.22*pow(1.0 - Umf_prev/Umfo, 2));
		
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
	}while(relError > relTol);
}

#endif
