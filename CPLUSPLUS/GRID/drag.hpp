#ifndef drag_hpp_
#define drag_hpp_

#include<cmath>
#include<string>
#include<iostream>

using namespace std;

template<typename realType> class drag
{
private:
// private data
	// operating pressure
	realType P = 0.1;
	
	// dimensionless operating pressure
	realType dimlessP = 1.0;
	
	// l exponent in Tomiyama's swarm
	realType lExponent = 0.0;

	// switch for swarm correction
	string swarmCorrection;
	
	// drag model functional
	string dragModel;
	
	// swarm correction model
	string swarmModel;

// private member functions
	// Drag Models
	realType MachDrag(const realType& Re, const realType& Eo) const;
	realType pureTomiyama(const realType& Re, const realType& Eo) const;
	realType slightlyTomiyama(const realType& Re, const realType& Eo) const;
	realType contaminatedTomiyama(const realType& Re, const realType& Eo) const;
	
	// Swarm Correction Models
	realType MachSwarm(const realType& Eo, const realType& egi) const;
	realType TomiyamaSwarm(const realType& egi) const;
	realType NoSwarmCorrection() const;
public:
// public memmber function
	// get user inputs
	virtual void getInputs();
	
	// overall drag coefficient
	realType DragCoefficient(const realType& Re, const realType& Eo, const realType& egi) const;
};

template<typename realType> void drag<realType>::getInputs()
{
	cout<<"\nUSER INPUTS FOR DRAG MODEL"<<endl;
	cout<<"--------------------------"<<endl;
	// drag model
	cout<<"The following drag models are available:"<<endl;
	cout<<"1. Tomiyama (1998) model for pure systems."<<endl;
	cout<<"2. Tomiyama (1998) model for slightly contaminated systems."<<endl;
	cout<<"3. Tomiyama (1998) model for fully contaminated systems."<<endl;
	cout<<"4. Mach (2020) model for contaminated polydisperse swarm."<<endl;

	cout<<"Please enter 1, 2, 3, or 4: ";
	int drag_model;
	cin>>drag_model;
	
	switch(drag_model)
	{
		case 1:
			dragModel = "pureTomiyama";
			cout<<"You've chosen Tomiyama's model for pure systems"<<endl;
			break;
		case 2:
			dragModel = "slightlyTomiyama";
			cout<<"You've chosen Tomiyama's model for slightly contaminated systems."<<endl;
			break;
		case 3:
			dragModel = "contaminatedTomiyama";
			cout<<"You've chosen Tomiyama's model for fully contaminated systems"<<endl;
			break;
		case 4:
			dragModel = "MachDrag";
			cout<<"You've chosen Mach's model for contaminated polydisperse systems."<<endl;
			break;
		default:
			cout<<"Numeric input between 1 and 4 expected."<<endl;
			throw;
	}

	// swarm correction
	cout<<"Would you like to use a swarm correction model? (yes/no): "<<endl;
	
	cin>>swarmCorrection;
	
	if(swarmCorrection == "yes" || swarmCorrection=="y")
	{
		cout<<"\nYou've chosen to use a swarm corrrection model."<<endl;
		cout<<"The following swarm correction models are available:"<<endl;
		cout<<"1. Tomiyama's swarm correction model."<<endl;
		cout<<"2. Mach(2020) swarm correction model for contaminated polydisperse systems."<<endl;
		
		cout<<"Please enter 1 or 2:";
		int swarm_model;
		cin>>swarm_model;
		switch(swarm_model)
		{
			case 1:
				swarmModel = "TomiyamaSwarm";
				cout<<"You've chosen Tomiyama's swarm correction model."<<endl;
				cout<<"Now enter the l exponent:";
				cin>>lExponent;
				cout<<"[RECEIVED] l exponent = "<<lExponent<<endl;
				break;
			case 2:
				swarmModel = "MachSwarm";
				cout<<"You've chosen Mach's swarm correction model."<<endl;
				cout<<"Now enter operating pressure [MPa]: ";
				cin>>P;
				dimlessP = P/0.1;
				cout<<"[RECEIVED] operating pressure = "<<P<<endl;
				break;
			default:
				cout<<"Numeric input (1 or 2) expected."<<endl;
				throw;
		}
	}
	else if (swarmCorrection == "no" || swarmCorrection=="n")
	{
		swarmModel = "NoSwarmCorrection";
		cout<<"\nYou've chosen to NOT use a swarm correction model."<<endl;
	}
	else
	{
		cout<<"Only acceptable inputs are 'yes' or 'no'"<<endl;
		throw;
	}
	
}

template<typename realType> realType drag<realType>::pureTomiyama(const realType& Re, const realType& Eo) const
{	
	// drag coefficient
	realType Cdinf = fmax(fmin(16.0/Re*(1.0 + 0.15*pow(Re, 0.687)), 48.0/Re), 8.0/3.0 * Eo/(Eo + 4.0));
	
	return Cdinf;
}

template<typename realType> realType drag<realType>::slightlyTomiyama(const realType& Re, const realType& Eo) const
{
	realType Cdinf = fmax(fmin(24.0/Re*(1.0 + 0.15*pow(Re, 0.687)), 72.0/Re), 8.0/3.0*Eo/(Eo + 4.0));
	return Cdinf;
}

template<typename realType> realType drag<realType>::contaminatedTomiyama(const realType& Re, const realType& Eo) const
{
	realType Cdinf = fmax(24.0/Re*(1.0 + 0.15*pow(Re, 0.687)), 8.0/3.0*Eo/(Eo + 4.0));
	return Cdinf;
}

template<typename realType> realType drag<realType>::MachDrag(const realType& Re, const realType& Eo) const
{
	realType Cdinf = 24.0/Re*(1.0 + pow(Re*Eo, 0.47));
	return Cdinf;
}

template<typename realType> realType drag<realType>::TomiyamaSwarm(const realType& egi) const
{
	realType fswarm = (1.0 - egi)*pow(1.0 - egi, 3.0 - lExponent);
	return fswarm;

}

template<typename realType> realType drag<realType>::MachSwarm(const realType& Eo, const realType& egi) const
{
	realType exponent = 0.24*(1.0 - exp(-0.13*dimlessP/Eo));
	realType fswarm = pow(1.0 - egi, exponent);
	return fswarm;
}

template<typename realType> realType drag<realType>::NoSwarmCorrection() const
{
	realType fswarm = 1.0;
	return fswarm;
}

template<typename realType> realType drag<realType>::DragCoefficient(const realType& Re, const realType& Eo, const realType& egi) const
{
	// single bubble drag coefficient
	realType Cdinf;
	if(dragModel == "pureTomiyama")
	{
		Cdinf = pureTomiyama(Re, Eo);
	}
	else if(dragModel == "slightlyTomiyama")
	{
		Cdinf = slightlyTomiyama(Re, Eo);
	}
	else if(dragModel == "contaminatedTomiyama")
	{
		Cdinf = contaminatedTomiyama(Re, Eo);
	}
	else if(dragModel == "MachDrag")
	{
		Cdinf = MachDrag(Re, Eo);
	}

	// swarm correction factor
	realType fswarm;
	if(swarmModel == "MachSwarm")
	{
		fswarm = MachSwarm(Eo, egi);
	}
	else if(swarmModel == "TomiyamaSwarm")
	{
		fswarm = TomiyamaSwarm(egi);
	}
	else if(swarmModel == "NoSwarmCorrection")
	{
		fswarm = NoSwarmCorrection();
	}

	// combined drag coefficient
	realType Cdi = Cdinf*fswarm;

	return Cdi;
}

#endif
