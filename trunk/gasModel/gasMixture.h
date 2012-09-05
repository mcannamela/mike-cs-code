#ifndef GAS_MIXTURE_H
#define GAS_MIXTURE_H
#include <vector>
#include <string>
#include <cmath>
#include "lookupTable.h"

using namespace std;

typedef vector<double> dVec;

class gasMixture //gasMixture class models a mixture of pure gasses, returning temperature dependent properties
{

	private:
		//data 
		lookupTable 			mixProperties;	//table to hold the mixture properties
		dVec					massFraction;	//mass fraction of component gasses, this is set second after volumeFraction
		dVec					volumeFraction;	//volume fraction of component gasses at some temperature, this is set first
		dVec					molFraction;	//mol fraction of component gasses, this is set last after massFraction
		dVec					molarMass;		//molar mass of component gasses
		vector<lookupTable>		pureGas;		//tables of properties for each component
		vector<string>			gasNames;		//name of each component
		
		//maps property name to line of table
		int						temperatureIdx;
		int						conductivityIdx; 
		int						densityIdx; 		
		int						enthalpyIdx;
		int						specificHeatIdx; 
		int						viscosityIdx;
		
		//flags
		bool					mixTempSet;
		//functions
		void				mixDensity();
		void				mixEnthalpy();
		void				mixSpecificHeat();
		void				mixConductivity();
		void				mixViscosity();								//ref. C.R. Wilke, "A Viscosity Equation for Gas Mixtures", Jou. Chemical Physics, vol 18 nr 4, Apr. 1950 
		void				computeMassFraction(double);
		void				computeMolFraction();							
		vector<dVec>		computePhi(dVec);							//ref. C.R. Wilke, "A Viscosity Equation for Gas Mixtures", Jou. Chemical Physics, vol 18 nr 4, Apr. 1950 
		dVec 				computeWilkeDenominatorSum(vector<dVec>);	//ref. C.R. Wilke, "A Viscosity Equation for Gas Mixtures", Jou. Chemical Physics, vol 18 nr 4, Apr. 1950 
				
	public:
		//constructors
		gasMixture();
		gasMixture(int);
		
		//setters
		void 		setVolumeFraction(dVec,double);//(vectorOfVolumeFractions, Temperature)
		void 		setMassFraction(double);//set the mass fraction given a temperature, with volume fraction already known
		void		setMixTemperature(dVec);//temperatures at which to evaluate mix properties
		void 		addComponent(char*,dVec,double,double);//(gasName, vectorOfVolumeFractions,Temperature), calls setMassFraction() and computeMixProperties()
		
		//getters, property functions take a temperature and return the property of the mixture at that temperature
		void		writePropertyFile(char* fName)
						{mixProperties.toFile(fName);}
		double 		Conductivity(double temperature)					//wilke method
						{return mixProperties.getVal(0,1,temperature);}
		double 		Density(double temperature)							//volume fraction average
						{return mixProperties.getVal(0,2,temperature);}
		double 		Enthalpy(double temperature)							//mass fraction average
						{return mixProperties.getVal(0,3,temperature);}
		double 		SpecificHeat(double temperature)					//mass fraction average
						{return mixProperties.getVal(0,4,temperature);}
		double 		Viscosity(double temperature)							// wilke method
						{return mixProperties.getVal(0,5,temperature);}
		double		Temperature(double enthalpy)
						{return mixProperties.getVal(3, 0, enthalpy);}
};
#endif