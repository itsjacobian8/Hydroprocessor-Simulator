#include "fluidization.hpp"

#include<boost/math/tools/config.hpp>
#include<boost/multiprecision/cpp_bin_float.hpp>

using namespace boost::multiprecision;
typedef cpp_bin_float_quad quadPrec;

int main()
{
	// check numeric limits for defined type
	cout<<"CHECKING NUMERIC LIMITS FOR CHOSEN FLOAT"<<endl;
	cout<<"----------------------------------------"<<endl;
	cout<<"minimum value = "<<numeric_limits<quadPrec>::min()<<endl;
	cout<<"maximum value = "<<numeric_limits<quadPrec>::max()<<endl;
	cout<<"machine precision = "<<numeric_limits<quadPrec>::epsilon()<<endl;
	
	// initialize class
	fluidization<quadPrec> testCase;

	// obtain user inputs
	testCase.getInputs();

	// update calculated fluidization parameters
	testCase.updateParameters();
	
	cout<<"\nMINIMUM FLUIDIZATION VELOCITIES"<<endl;
	cout<<"-------------------------------"<<endl;
	
	// compute gas free minimum fluidization velocity
	testCase.LSMinimumFluidization();
	cout<<"gas free minimum fluidization velocity = "<<testCase.Umfo<<endl;
	
	// compute gas-liquid-solid minimum fluidization velocity
	quadPrec Ugbed = 0.072;
	testCase.GLSMinimumFluidization(Ugbed);
	cout<<"gas-liquid-solid minimum fluidization velocity = "<<testCase.Umf<<endl;

	return 0;
}
