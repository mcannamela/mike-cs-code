#include "gasMixture.h"
#include <iostream>
#include <vector>

using namespace std;
typedef vector<double> dVec;

int main(void)
{	
	gasMixture g;
	g = gasMixture(2);
	cout << "\n"<< g.Density(10000) << "\n";
	cout << g.Enthalpy(10000) << "\n";
	cout << g.Conductivity(10000) << "\n";
	
	return 0;
};