#ifndef bed_hpp_
#define bed_hpp_

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
	void wallEffect();
	// computes particle-liquid terminal velocity
	void terminalVelocity();
	// computes Richardson-Zaki coefficient
	void RichardsonZaki();

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
	// gets user inputs
	virtual void getInputs();

	// updates calculated parameters
	virtual void updateParameters();

	// computes particles holdup
	void particlesHoldup(const realType& Ugbed, const realType& Ulbed);

	// computes bed height
	void bedHeight();
};

template<typename realType> void bed<realType>::getInputs()
{
	// first get fluidization parameters
	fluidization<realType>::getInputs();

	// now get remaining bed parameters
	cout<<"\nUSER INPUTS FOR REMAINING BED PARAMETERS"<<endl;
	cout<<"----------------------------------------"<<endl;
	cout<<"Enter column diameter [m]: ";
	cin>>Dc;
	cout<<"RECEIVED] column diameter = "<<Dc<<endl;

	cout<<"\nEnter recycle line diameter [m]: ";
	cin>>drec;
	cout<<"[RECEIVED] recycle line diameter = "<<drec<<endl;
	
	cout<<"\nEnter mass of catalyst inventory [kg]: ";
	cin>>mcat;
	cout<<"[RECEIVED] mass of catalyst inventory = "<<mcat<<endl;
	
	Ac = 0.25*this->pi*(pow(Dc, 2) - pow(drec, 2));
	cout<<"[CALCULATED] column's cross-sectional area = "<<Ac<<endl;
}

template<typename realType> void bed<realType>::wallEffect()
{
	k = 1.0 - 1.33*(this->dvol/Dc);
}

template<typename realType> void bed<realType>::RichardsonZaki()
{
	realType A = 0.043*pow(this->Arp, 0.57)*(1.0 - 1.24*pow(this->dvol/Dc, 0.27));
	n = (4.8 + 2.4*A)/(1.0 + A);
}

template<typename realType> void bed<realType>::terminalVelocity()
{
	realType Repinf = cbrt(this->Arp)/(18.0/pow(this->Arp, 2.0/3.0) + (2.335 - 1.744*this->phi)/pow(this->Arp, 1.0/6.0));

	upinf = (Repinf*this->mul)/(this->rhol*this->dvol);
}

template<typename realType> void bed<realType>::updateParameters()
{
	// update fluidization parameters first
	fluidization<realType>::updateParameters();

	// now update bed parameters
	cout<<"\nCALCULATED BED PARAMETERS"<<endl;
	cout<<"-------------------------"<<endl;

	wallEffect();
	cout<<"wall effect parameter = "<<k<<endl;
	
	RichardsonZaki();
	cout<<"Richardson-Zaki coefficient = "<<n<<endl;

	terminalVelocity();
	cout<<"terminal velocity = "<<upinf<<endl;
}

template<typename realType> void bed<realType>::particlesHoldup(const realType& Ugbed, const realType& Ulbed)
{
	epsilonp = 1.0 - (pow(Ulbed/(k*upinf), 1.0/n)*(1.0 + 0.22*pow(Ugbed/Ulbed, 0.92)));

}

template<typename realType> void bed<realType>::bedHeight()
{
	hbed = mcat/(this->rhop*Ac*epsilonp);
}

#endif
