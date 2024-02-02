ifndef BUBBLESIZEDISTRIBUTION_HPP
#define BUBBLESIZEDISTRIBUTION_HPP

#include "drag.hpp"

#include<string>
#include<vector>
#include<boost/math/tools/config.hpp>
#include<boost/math/distributions/lognormal.hpp>

template<typename realType> class sizeDistribution : public drag<realType>
{
private:
// private data	
	// gas density
	realType rhog;
	// liquid density
	realType rhol;
	// liquid viscosity
	realType mul;
	// surface tension
	realType sigma;
	// name of drag model
	string dragModel;
	// name of swarm correction model
	string swarmModel;
	// name of bubble aspect ratio model
	string aspectRatioModel;

	// gravitational constant
	const realType g = 9.81;

// private member functions
	/*Aspect ratio model proposed in Wellek et al. (1966)
	* Model was based on single isolated bubbles 
	*/
	realType Wellek
	(
		const realType Eo
	) const;

	/*Aspect ratio model proposed in Besagni & Deen (2020)
	* -----------------------------------------------------
	* Model was based on a large dataset of single isolated bubbles
	* and several bubbles rising together under different conditions. 
	* Maximum gas holdup was less than 11%.
	* Some gas holdups were not reported (possibly for single isolated bubbles),
	* so the maximum gas holdup could be higher than 11%. 
	* Consult individual works for the accuracy of experimental methods
	* used to determine bubble shape, both in dilute and dense bubbly flows
	*/
	realType BesagniDeen
	(
		const realType& Re, 
		const realType& Eo
	) const;

	// converts chord length to volume equivalent diameter
	realType convertChordToDiameter
	(
		const realType& cbi
	) const;

	// iteratively computes bubble slip velocity
	realType slipVelocity
	(
		const realType& dbi, 
		const realType& egi
	) const;

	// objective function for iteratively solving for class holdup
	realType objectiveFunction
	(
		const realType& egi, 
		const realType& dbi,
		const realType& Ugbed,
		const realType& Ulbed
	) const;

	// custom bisection method
	realType bisect
	(
		const realType& dbi,
		const realType& Ugbed,
		const realType& Ulbed
	) const;
public:
// public data
	// vector for storing cbi, dbi, usi, egi, & pdf
	vector<vector<realType>> bsvd(2, vector<realType>(5, 0.));

// public member functions
	// constructor
	bubbleSizeDistribution
	(
		realType _rhog, 
		realType _rhol, 
		realType _mul, 
		realType _sigma,
		realType _P,
		realType _lExponent,
		string _dragModel,
		string _swarmModel,
		string _aspectRatioModel
	);

	// returns a vector
	vector<vector<realType>> sizeDistribution
	(
		const realType& Ugbed, 
		const realType& Ulbed,  
		const realType& mu, 
		const realType& s,
		const realType& classWidth
	)const;

	// modifies a dynamic array that is passed as a pointer
	void sizeDistribution
	(
		const realType& Ugbed,
		const realType& Ulbed,
		const realType& mu,
		const realType& s,
		const realType& cbiMin,
		const realType& classWidth,
		const int& nClasses,
		realType **bsvd
	)const;

	// modifies a vector that is passed by reference
	void sizeDistribution
	(
		const realType& Ugbed,
		const realType& Ulbed,
		const realType& mu,
		const realType& s,
		const realType& cbiMin,
		const realType& classWidth,
		const int& nClasses,
		vector<vector<realType>>& bsvd
	)const;

	// modifies a vector that is declared as a public member
	void sizeDistribution
	(
		const realType& Ugbed,
		const realType& Ulbed,
		const realType& mu,
		const realType& s,
		const realType& classWidth
	);
};

template<typename realType> bubbleSizeDistribution<realType>::bubbleSizeDistribution
(
	realType _rhog,
	realType _rhol,
	realType _mul,
	realType _sigma,
	realType _P,
	realType _lExponent,
	string _dragModel,
	string _swarmModel,
	string _aspectRatioModel
)
:
drag<realType>::drag(_P, _lExponent, _dragModel, _swarmModel),
rhog(_rhog),
rhol(_rhol),
mul(_mul),
sigma(_sigma),
aspectRatioModel(_aspectRatioModel)
{
}

template<typename realType> realType bubbleSizeDistribution<realType>::Wellek
(
	const realType& Eo
) const
{
	realType E = 1.0/(1.0 + 0.163*pow(Eo, 0.757));
	return E;
}

template<typename realType> realType bubbleSizeDistribution<realType>::BesagniDeen
(
	const realType& Re,
	const realType& Eo
) const
{
	realType E = pow(1.0 + 0.45*(Re*Eo), -0.08);
	return E;
}

template<typename realType> realType bubbleSizeDistribution<realType>::convertChordToDiameter
(
	const realType& cbi,
	const realType& usi,
)const
{
	// bubble Eotvos number
	realType Eo = ((rhol - rhog)*g*pow(0.0015*cbi, 2))/sigma;
	// bubble Reynolds number
	realType Re = (rhol*usi*0.0015*cbi)/mul;
	// bubble aspect ratio
	realType E; 
	if(aspectRatioModel == "Wellek") E = Wellek(Eo);
	if(aspectRatio == "BesagniDeen") E = BesagniDeen(Re,Eo);
	// initial estimate of bubble diameter
	dbi = 1.5*cbi*pow(E, -2.0/3.0);
	// tolerance
	realType tol = numeric_limits<realType>::epsilon();
	
	realType dbi_prev, relError;
	do
	{
		// assign current value to previous value
		dbi_prev = dbi;
		// update bubble Reynolds number
		Re = (rhol*usi*0.001*dbi_prev)/mul;
		// update bubble Eotvos number
		Eo = ((rhol - rhog)*g*pow(0.001*dbi, 2))/sigma;
		// update aspect ratio
		if(aspecRatioModel == "Wellek") E = Wellek(Eo);
		if(aspectRatioModel == "BesagniDeen") E = BesagniDeen(Re, Eo);
		// update bubble diameter
		dbi = 1.5*cbi*pow(E, -2.0/3.0);
		// update relative error
		relError = fabs((dbi - dbi_prev)/dbi);
	}while(relError > tol);
	return dbi;
}

template<typename realType> realType bubbleSizeDistribution<realType>::slipVelocity
(
	const realType& dbi,
	const realType& egi
) const
{
	// initialize slip velocity
	realType usi = 0.0001;
	// tolerance
	realType tol = numeric_limits<realType>::epsilon();

	realType usi_prev, Re, Eo, E, Cdi, relError;
	do
	{
		// assign current value to previous value
		usi_prev = usi;		
		// bubble Reynolds number
		Re = (rhol*usi*0.001*dbi)/mul;
		// bubble Eotvos number
		Eo = ((rhol - rhog)*g*pow(0.001*dbi, 2))/sigma;
		// bubble aspect ratio
		if(aspectRatioModel == "Wellek") E = Wellek(Eo);
		if(aspectRatioModel == "BesagniDeen") E = BesagniDeen(Re, Eo);
		// drag coefficient
		Cdi = drag<realType>::DragCoefficient(Re, Eo, egi);
		// update bubble slip velocity
		usi = sqrt(4.0/3.0*(g*0.001*dbi)/Cdi*pow(E, 2.0/3.0)*(rhol - rhog)/rhol);
		// update relative error
		relError = fabs((usi - usi_prev)/usi);
	}while(relError > tol);	
	return usi;
}

template<typename realType> realType bubbleSizeDistribution<realType>::objectiveFunction
(
	const realType& egi,
	const realType& dbi,
	const realType& Ugbed,
	const realType& Ulbed
) const
{
	// bubble slip velocity
	realType usi = slipVelocity(dbi, egi);

	// interstitial gas velocity assuming:
	// 1. negligible slip between phases
	// 2. all gas flows through the selected bubble class.
	// 
	// Assumption 1. will be corrected by forcing
	// the slip between the phases to be equal to
	// the bubble slip velocity calculated above
	//
	// Assumption 2. will be corrected by weighing the
	// contributions of individual classes to overall 
	// gas holdup and gas-liquid separation efficiency 
	realType ugas = Ugbed/fmax(egi, 1.0e-8);

	// interstitial liquid velocity 
	realType uliquid = Ulbed/fmax(1.0 - egi, 1.0e-8);
	
	// slip between the phases
	realType uslip = ugas - uliquid; 
	
	// difference between bubble slip velocity and the slip between the phases.
	// This will be minimized to correct one of the assumptions.
	// note that the two slip velocities can be interlinked via class holdup,
	// and although the effects of above assumptions on quantities
	// of interest (overall gas holdup, separation efficiency, etc)
	// are corrected, bubble slip velocities should also be scaled by 
	// their respective probability densities if they are of interest
	realType error = usi - uslip;
	
	return error;
}

template<typename realType> realType bubbleSizeDistribution<realType>::bisect
(
	const realType& dbi,
	const realType& Ugbed,
	const realType& Ulbed
) const
{
	// initial solution bracket
	// gas holdup must be between 0 and 1
	realType egi_min = 0.0;
	realType egi_max = 1.0;
	// machine precision
	realType tol = numeric_limits<realType>::epsilon();
	// absolute error
	realType absError = egi_max - egi_min;

	reaType minEval, midEval, egi;
	while(absError > tol)
	{
		// mid point
		egi = 0.5*(egi_max + egi_min);
		// midpoint evaluation of the objective function
		midEval = objectiveFunction(egi, dbi, Ugbed, Ulbed);
		// check if an exact root has been located
		if(midEval == 0.) break;
		// lower end evaluation of the objective function
		minEval = objectiveFunction(egi_min, dbi, Ugbed, Ulbed);
		// adjust solution bracket
		if(minEval*midEval < 0.) egi_max = egi;
		if(minEval*midEval > 0.) egi_min = egi;
		// update absolute error
		absError = egi_max - egi_min;
	}
	
	return egi;
}

#endif
