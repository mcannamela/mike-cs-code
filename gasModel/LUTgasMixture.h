#ifndef GAS_MIXTURE_H
#define GAS_MIXTURE_H
#include <vector>
#include <string>
#include <cmath>
#include "lookupTable.h"
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;
using namespace std;

typedef vector<double> dVec;

class gasMixture //gasMixture class models a mixture of pure gasses, returning temperature dependent properties
{

	private:
		//data 
		lookupTable 			mixProperties;	//table to hold the mixture properties
		
		//maps property name to line of table
		int						temperatureIdx;
		int						conductivityIdx; 
		int						densityIdx; 		
		int						enthalpyIdx;
		int						specificHeatIdx; 
		int						viscosityIdx;
						
	public:
		//constructors
		gasMixture()
		{
			init();
			//mixProperties = lookupTable("argon.txt",6);
		}
		gasMixture(fs::path p)
		{
		init();
		mixProperties = lookupTable(p,6);
		}
		
		void init(void)
		{
			temperatureIdx 		= 0;
			conductivityIdx		= 1;
			densityIdx 				= 2;
			enthalpyIdx				= 3;
			specificHeatIdx		= 4;
			viscosityIdx				= 5;
		}
		
		//getters, property functions take a temperature and return the property of the mixture at that temperature
		double 		Conductivity(double temperature)					
						{return mixProperties.getVal(temperatureIdx, conductivityIdx, temperature);}
		double 		Density(double temperature)							
						{return mixProperties.getVal(temperatureIdx, densityIdx, temperature);}
		double 		Enthalpy(double temperature)							
						{return mixProperties.getVal(temperatureIdx, enthalpyIdx, temperature);}
		double 		SpecificHeat(double temperature)					
						{return mixProperties.getVal(temperatureIdx, specificHeatIdx, temperature);}
		double 		Viscosity(double temperature)							
						{return mixProperties.getVal(temperatureIdx, viscosityIdx, temperature);}
		double		Temperature(double enthalpy)
						{return mixProperties.getVal(3, temperatureIdx, enthalpy);}
};


#endif