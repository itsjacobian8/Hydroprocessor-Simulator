#include "drag.hpp"
#include<limits>
#include<boost/math/tools/config.hpp>
#include<boost/multiprecision/cpp_bin_float.hpp>

using namespace boost::multiprecision;

typedef cpp_bin_float_quad quadPrec;

int main()
{
	// checking limits
	cout<<"\nCHECKING LIMITS FOR CHOSEN FLOAT"<<endl;
	cout<<"--------------------------------"<<endl;
	cout<<"minimum value = "<<numeric_limits<quadPrec>::min()<<endl;
	cout<<"maximum value = "<<numeric_limits<quadPrec>::max()<<endl;
	cout<<"machine precision = "<<numeric_limits<quadPrec>::epsilon()<<endl;
	
	// initializing class
	drag<quadPrec> testCase;
	
	// obtain user inputs
	testCase.getInputs();

	// compute drag coefficient
	quadPrec Re = 100.0;
	quadPrec Eo = 1.0;
	quadPrec egi = 0.20;
	
	quadPrec Cdi = testCase.DragCoefficient(Re, Eo, egi);
	cout<<"\ndrag coefficient = "<<Cdi<<endl;

	return 0;
}
