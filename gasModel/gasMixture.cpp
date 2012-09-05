#include "gasMixture.h"

gasMixture::gasMixture()
{
	temperatureIdx 		= 0;
	conductivityIdx		= 1;
	densityIdx 			= 2;
	enthalpyIdx			= 3;
	specificHeatIdx		= 4;
	viscosityIdx		= 5;
	
	mixTempSet = false;
};


gasMixture::gasMixture(int anInt)
{
	temperatureIdx 		= 0;
	conductivityIdx		= 1;
	densityIdx 			= 2;
	enthalpyIdx			= 3;
	specificHeatIdx		= 4;
	viscosityIdx		= 5;
	

	int			nEntries = 120; 
	double		defaulttemperature []	=	{0, 300, 500,  600,  700,  800,  900,  1000,  1100,  1200,  1300,  1400,  1500,  1600,  1700,  1800,  1900,  2000,  2100,  2200,  2300,  2400,  2500,  2600,  2700,  2800,  2900,  3000,  3100,  3200,  3300,  3400,  3500,  3600,  3700,  3800,  3900,  4000,  4300,  4500,  4600,  4700,  4800,  4900,  5000,  5100,  5200,  5300,  5400,  5500,  6000,  6300,  6500,  7000,  7500,  8000,  8200,  8500,  8800,  9000,  9100,  9200,  9300,  9500,  9800,  10000,  10200,  10300,  10500,  10600,  10700,  10800,  10900,  11000,  11100,  11200,  11300,  11400,  11500,  11800,  12000,  12300,  12500,  12700,  12800,  13000,  13200,  13300,  13500,  13700,  13800,  14000,  14200,  14300,  14500,  14800,  15000,  15200,  15300,  15500,  16000,  16100,  16200,  16300,  16400,  16500,  16800,  17000,  17300,  17500,  17800,  18000,  18300,  18500,  18800,  19000,  19500,  19700,  19800,  20000};
	dVec		defaultTemperature = arr2vec(defaulttemperature,nEntries );
	
			double		T = 300;
			double		molarMassAr = .0395; //kg/mol
			double		molarMassHe = .004;
			dVec		vFrac1(1,1);
			dVec		vFrac2(2,.5);
			char*		Ar = "argon.txt";
			char*		He = "helium.txt";

	switch(anInt)
	{
		case 0://simple gas for testing eqn solvers
			nEntries = 4;
			defaultTemperature.clear();
			for(int i=0; i<nEntries;i++)
			defaultTemperature.push_back((double) i/(nEntries-1));
			setMixTemperature(defaultTemperature);
			mixTempSet = true;
			addComponent("simpleGas.txt",vFrac1,.5,1);//ficticious gas with constant, round number properties
			break;
		case 1://pure argon
			setMixTemperature(defaultTemperature);
			mixTempSet = true;
			addComponent(Ar,vFrac1,T,molarMassAr);
			break;
		case 2://pfender 32% He case
			setMixTemperature(defaultTemperature);
			mixTempSet = true;
			addComponent(Ar,vFrac1,T,molarMassAr);
			vFrac2[1] = .32;			//<--change He volume frac here
			//cout << 1-vFrac2[1];
			vFrac2[0] = 1-vFrac2[1];
			//cout <<vFrac2[0];
			addComponent(He,vFrac2,T,molarMassHe);
			break;
		
	}
	
};

void gasMixture::mixDensity()
{
	int 	nComponents = pureGas.size();
	dVec	Temperature = mixProperties.getLine(0);
	int 	nEntries	= Temperature.size();
	dVec	theProp(nEntries,0);//hold the line of the table for the computed mixed property
	for (int i = 0; i<nComponents ; i++)//loop over the component gasses
	{
		for (int j = 0; j<nEntries ; j++)//loop over all entries
		{
			theProp[j] += volumeFraction[i]*pureGas[i].getVal(0,densityIdx,Temperature[j] ); //volume weighted average of densities;
		};
	};
	mixProperties.setLine(densityIdx,&theProp);//set the line of the mixProperties table
};

void gasMixture::mixEnthalpy()
{
	int 	nComponents = pureGas.size();
	dVec	Temperature = mixProperties.getLine(0);
	int 	nEntries	= Temperature.size();
	dVec	theProp(nEntries,0);
	for (int i = 0; i<nComponents ; i++)
	{
		for (int j = 0; j<nEntries ; j++)
		{
			theProp[j]+= massFraction[i]*pureGas[i].getVal(0,enthalpyIdx,Temperature[j]); //mass weighted average of enthalpy
		};
	};
	mixProperties.setLine(enthalpyIdx,&theProp);
};

void gasMixture::mixSpecificHeat()
{
	int 	nComponents = pureGas.size();
	dVec	Temperature = mixProperties.getLine(0);
	int 	nEntries	= Temperature.size();
	dVec	theProp(nEntries,0);
	for (int i = 0; i<nComponents ; i++)
	{
		for (int j = 0; j<nEntries ; j++)
		{
			theProp[j]+= massFraction[i]*pureGas[i].getVal(0,specificHeatIdx,Temperature[j]); //mass weighted average of specific heat
		};
	};
	mixProperties.setLine(specificHeatIdx,&theProp);
};

void gasMixture::mixConductivity()
{//use wilke's method to combine conductivities
	int 			nComponents = pureGas.size();
	dVec			Temperature = mixProperties.getLine(0);
	int 			nEntries	= Temperature.size();
	dVec			theProp(nEntries,0);
	vector<dVec> 	phi;
	dVec			denomSum;//holds the sum in the denominator of wilke's formula
	dVec			k;//hold the conductivities for a particular temperature
	
	for (int i = 0; i<nEntries ; i++)//for every entry of the table, compute a mixed value
	{
		k.clear();
		phi.clear();
		
		for (int j = 0; j<nComponents ; j++)//store component conductivities in k
			k.push_back( pureGas[j].getVal(0,conductivityIdx,Temperature[i]) );
		
		phi = computePhi(k);//compute the phi for these component k's
		denomSum = computeWilkeDenominatorSum(phi);//compute the inner product of phi with the mol fractions
		for (int j = 0; j<nComponents ; j++)
		{
				theProp[i]+= k[j]/(1+denomSum[j]/molFraction[j]);//ref wilke eqn 13
		};
	};
	mixProperties.setLine(conductivityIdx,&theProp);
};

void gasMixture::mixViscosity()
{//use wilke's method to combine conductivities
	int 			nComponents = pureGas.size();
	dVec			Temperature = mixProperties.getLine(0);
	int 			nEntries	= Temperature.size();
	dVec			theProp(nEntries,0);
	vector<dVec> 	phi;
	dVec			denomSum;//holds the sum in the denominator of wilke's formula
	dVec			mu;//hold the viscosities for a particular temperature
	
	for (int i = 0; i<nEntries ; i++)//for each temperature in the table, compute the mixed viscosity
	{
		mu.clear();
		phi.clear();
		
		for (int j = 0; j<nComponents ; j++)
			mu.push_back( pureGas[j].getVal(0,viscosityIdx,Temperature[i]) );//get property values at current temperature
		
		phi = computePhi(mu);
		denomSum = computeWilkeDenominatorSum(phi);
		for (int j = 0; j<nComponents ; j++)
		{
			theProp[i]+= mu[j]/(1+denomSum[j]/molFraction[j]);
		};
	};
	mixProperties.setLine(viscosityIdx,&theProp);
};

vector<dVec> gasMixture::computePhi(dVec theProps)
{
	//compute the matrix phi from wilke's 1950 paper
	int 					nComponents = pureGas.size();
	dVec					z(nComponents,0);
	vector<dVec>	phi;

	//initialize phi to a vector of zero vectors
	phi.clear();
	for (int i=0; i<nComponents; i++)
	phi.push_back(z);
	
	//compute each phi(i,j)
	for (int i=0; i<nComponents; i++)
	{
		for (int j=0; j<nComponents; j++)
		{
		if (i==j)
			continue;
		else
		phi[i][j] = pow( 1+pow(theProps[i]/theProps[j],.5)*pow(molarMass[j]/molarMass[i],.25) ,2)/pow(8*(1+molarMass[i]/molarMass[j]),.5);
		};
	};
	return phi;
};

dVec gasMixture::computeWilkeDenominatorSum(vector<dVec> phi)
{//compute the double sum denominator in eqn 13 of wilke's paper
	int 	nComponents = pureGas.size();
	dVec	denomSum(nComponents,0);
	for (int i=0; i<nComponents; i++)
	{
		for (int j=0; j<nComponents; j++)
		{
			denomSum[i] += molFraction[j]*phi[i][j];
		};
	};
	return denomSum;
};

void gasMixture::computeMassFraction(double T)//called whenever volumeFraction is changed
{
	int 	nComponents = pureGas.size();
	double 	total = 0;
	massFraction.clear();
	for (int i = 0; i<nComponents ; i++)//loop over all component gasses
	{
		massFraction.push_back( volumeFraction[i]*pureGas[i].getVal(0,densityIdx,T) ); //weight each density at temperature T by the volume fraction to obtain the mass fraction
		total+=massFraction[i];
	};
	for (int i = 0; i<nComponents ; i++)//normalize
			massFraction[i]=massFraction[i]/total;	
			
	computeMolFraction();
};

void gasMixture::computeMolFraction()//called whenever massFraction is changed
{
	int 	nComponents = pureGas.size();
	double 	total = 0;
	molFraction.clear();
	
	for (int i = 0; i<nComponents ; i++)//loop over all component gasses
	{
		molFraction.push_back( massFraction[i]/molarMass[i] ); //mols are the molar mass times mass
		total+=molFraction[i];
	};
	for (int i = 0; i<nComponents ; i++)//normalize
			molFraction[i]=molFraction[i]/total;	
};

//public fns
void gasMixture::setMassFraction(double T)
{
	computeMassFraction(T);
};
	
void gasMixture::setVolumeFraction(dVec vFrac, double T)//(vectorOfVolumeFractions, Temperature)
{
	volumeFraction.clear();
	
	for (unsigned int i =0 ; i<vFrac.size() ; i++)
	volumeFraction.push_back(vFrac[i]);
	
	computeMassFraction(T);
};

void gasMixture::setMixTemperature(dVec T)//temperatures at which to evaluate mix properties
{
	mixProperties.setLine(0,&T);
	mixTempSet = true;
	
	if (pureGas.size()>0)//update mix properties if we already have component gasses
	{
		mixDensity();
		mixEnthalpy();
		mixConductivity();
		mixSpecificHeat();
		mixViscosity();
	}
};

void gasMixture::addComponent(char* gasName, dVec vFrac, double T, double mMass)
{	
	string	gName(gasName);
	int 	nProps = 6;
	lookupTable gasTable;
	
	molarMass.push_back(mMass);//need to know the molar mass for later
	gasNames.push_back(gName);//keep track of what components are called
	
	try
		{gasTable.fromFile(gasName,nProps);}
	catch(char* eStr) //probably a bad filename
		{cout << eStr;};
	
	pureGas.push_back(gasTable);
	setVolumeFraction(vFrac,T);//if a gas is added then the volume fractions must change accordingly
	
	if (mixTempSet)
	{
		mixDensity();
		mixEnthalpy();
		mixConductivity();
		mixSpecificHeat();
		mixViscosity();
	}
};