#include<cmath>
#include<limits>
#include<iostream>
#include "drag.hpp"

#include<boost/math/tools/config.hpp>
#include<boost/multiprecision/cpp_bin_float.hpp>

using namespace std;
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
	
	// dimensionless numbers
	quadPrec Re = 100.0;
	quadPrec Eo = 1.0;
	
	// gas holdup
	quadPrec egi = 0.20;
	
	// lift coefficient in Tomiyama's swarm model
	quadPrec l = 0.0;

	// initializing drag class
	drag<quadPrec> testCase;
	
	// calculate single bubble drag coefficient
	quadPrec Cdinf = testCase.contaminatedTomiyama(Re, Eo);

	// calculate swarm correction factor
	quadPrec fswarm = testCase.TomiyamaSwarm(egi, l);

	// calculate combined drag coefficient
	quadPrec Cdi = Cdinf*fswarm;
	cout<<"\ndrag coefficient = "<<Cdi<<endl;

	return 0;
}
