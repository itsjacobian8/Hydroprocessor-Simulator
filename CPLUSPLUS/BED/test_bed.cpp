#include "bed.hpp"

#include<boost/math/tools/config.hpp>
#include<boost/multiprecision/cpp_bin_float.hpp>

using namespace boost::multiprecision;

typedef cpp_bin_float_quad quadPrec;

int main()
{
	// checking limits for selected type
	cout<<"\nCHECKING LIMITS FOR CHOSEN FLOAT"<<endl;
	cout<<"---------------------------------"<<endl;
	cout<<"minimum value = "<<numeric_limits<quadPrec>::min()<<endl;
	cout<<"maximum value = "<<numeric_limits<quadPrec>::max()<<endl;
	cout<<"machine precision = "<<numeric_limits<quadPrec>::epsilon()<<endl;

	// initializing class
	bed<quadPrec> testCase;

	// obtain user inputs
	testCase.getInputs();

	// update calculated parameters
	testCase.updateParameters();

	// calculate particles holdup
	quadPrec Ugbed = 0.072;
	quadPrec Ulbed = 0.027;
	testCase.particlesHoldup(Ugbed, Ulbed);
	cout<<"\nBED PROPERTIES"<<endl;
	cout<<"--------------"<<endl;
	cout<<"particles holdup = "<<testCase.epsilonp<<endl;
	
	// calculate bed height
	testCase.bedHeight();
	cout<<"bed height [m] = "<<testCase.hbed<<endl;
	
	return 0;
}
