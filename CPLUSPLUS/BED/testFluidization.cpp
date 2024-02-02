#include<cmath>
#include<limits>
#include<string>
#include<iostream>

#include "../misc/colorize.hpp"
#include "fluidization.hpp"

#include<boost/math/tools/config.hpp>
#include<boost/multiprecision/cpp_bin_float.hpp>

using namespace std;
using namespace boost::multiprecision;

typedef cpp_bin_float_quad quadPrec;

// struct for holding input parameters
template<typename realType> struct inputs{
	realType rhog, rhol, mul, rhop, Vp, Ap, Ugbed;
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
	fluidization<quadPrec> testCase
	(	params.rhog, 
		params.rhol, 
		params.mul, 
		params.rhop,
		params.Vp,
		params.Ap
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
	colorizedInput(params.Ugbed, "superficial gas velocity [m/s]: ");
}	
