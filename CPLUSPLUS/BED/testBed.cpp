#include<cmath>
#include<limits>
#include<string>
#include<iostream>

#include "../misc/colorize.hpp"
#include "bed.hpp"

#include<boost/math/tools/config.hpp>
#include<boost/multiprecision/cpp_bin_float.hpp>

using namespace std;
using namespace boost::multiprecision;

typedef cpp_bin_float_quad quadPrec;

// struct for holding input parameters
template<typename realType> struct inputs{
	realType rhog, rhol, mul, rhop, Vp, Ap, mcat, Dc, drec, Ugbed, Ulbed;
};

// obtains input parameters from user
template<typename realType> void getInputs(inputs<realType>& params);

int main()
{
	// check numeric limits for defined type
	colorizedHeading("CHECKING NUMERIC LIMITS FOR SELECTED FLOAT.");
	colorizedStringOutput("minimum value [-] = ");
	colorizedNumericOutput(numeric_limits<quadPrec>::min());
	colorizedStringOutput("maximum value [-] = ");
	colorizedNumericOutput(numeric_limits<quadPrec>::max());	
	colorizedStringOutput("machine precision [-] = ");
	colorizedNumericOutput(numeric_limits<quadPrec>::epsilon());
	
	// obtain user inputs
	inputs<quadPrec> params;
	getInputs<quadPrec>(params);

	// initialize class
	bed<quadPrec> testCase
	(	params.rhog, 
		params.rhol, 
		params.mul, 
		params.rhop,
		params.Vp,
		params.Ap,
		params.mcat,
		params.Dc,
		params.drec
	);
	
	// update calculated parameters
	testCase.updateParameters();

	// compute minimum fluidization velocities
	testCase.LSMinimumFluidization();	
	testCase.GLSMinimumFluidization(params.Ugbed);	
	
	colorizedHeading("\nMINIMUM FLUIDIZATION VELOCITIES");	
	colorizedStringOutput("gas free minimum fluidization velocity [m/s] = ");
	colorizedNumericOutput(testCase.Umfo);	
	colorizedStringOutput("gas-liquid-solid minimum fluidization velocity [m/s] = ");
	colorizedNumericOutput(testCase.Umf);
	
	// compute bed properties
	testCase.particlesHoldup(params.Ugbed, params.Ulbed);
	testCase.bedHeight();

	colorizedHeading("\nCALCULATED BED PROPERTIES");
	colorizedStringOutput("particles holdup [-] = ");
	colorizedNumericOutput(testCase.epsilonp);
	colorizedStringOutput("bed height [-] = ");
	colorizedNumericOutput(testCase.hbed);
	return 0;
}

template<typename realType> void getInputs(inputs<realType>& params)
{
	colorizedHeading("\nUSER INPUTS FOR FLUIDIZATION CLASS");
	colorizedInput(params.rhog, "gas density [kg/m^3]: ");
	colorizedInput(params.rhol, "liquid density [kg/m^3]: ");
	colorizedInput(params.mul, "liquid viscosity [Pa.s]: ");
	colorizedInput(params.rhop, "particle density [kg/m^3]: ");
	colorizedInput(params.Vp, "particle volume [m^3]: ");
	colorizedInput(params.Ap, "particle surface area [m^2]: ");
	colorizedInput(params.mcat, "mass of catalyst inventory [kg]: ");
	colorizedInput(params.Dc, "column diameter [m]: ");
	colorizedInput(params.drec, "recycle line diameter [m]: ");
	colorizedInput(params.Ugbed, "superficial gas velocity [m/s]: ");
	colorizedInput(params.Ulbed, "superficial liquid velocity [m/s]: ");
}	
