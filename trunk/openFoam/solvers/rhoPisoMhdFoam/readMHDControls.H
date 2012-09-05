#include <string>
const dictionary& MHDControls = mesh.solutionDict().subDict("MHDControls");

//correctors and switches
int nBcorr(readInt(MHDControls.lookup("nCorrectors")));
int EM_scheme(readInt(MHDControls.lookup("EM_scheme")));
int J_scheme(readInt(MHDControls.lookup("J_scheme")));
bool lorenzForceOn(readBool(MHDControls.lookup("lorenzForceOn")));
bool radiationOn(readBool(MHDControls.lookup("radiationOn")));
bool fullOutputOn(readBool(MHDControls.lookup("fullOutputOn")));

//clamp values
dimensionedScalar pMax(MHDControls.lookup("pMax"));
dimensionedScalar pMin(MHDControls.lookup("pMin"));
dimensionedScalar jouleHeatingMax(MHDControls.lookup("jouleHeatingMax"));
dimensionedScalar cathodeJouleHeatingMax(MHDControls.lookup("cathodeJouleHeatingMax"));

string gasDir(MHDControls.lookup("gasDir"));
scalar gasMolarMass(readScalar(MHDControls.lookup("gasMolarMass")));

//whether patches should fix the internal value of sigma according to the wall temperature
Foam::wordList patches= MHDControls.subDict("patches").toc();

//parameters of conductivity doping for fixed attachment point model
bool dopingOn(readBool(MHDControls.subDict("sigmaDoping").lookup("dopingOn")));
bool dopingTwoSided(readBool(MHDControls.subDict("sigmaDoping").lookup("dopingTwoSided")));
scalar xDope(readScalar(MHDControls.subDict("sigmaDoping").lookup("xDope")));
scalar yDope(readScalar(MHDControls.subDict("sigmaDoping").lookup("yDope")));
scalar zDope(readScalar(MHDControls.subDict("sigmaDoping").lookup("zDope")));
scalar dopingRadius(readScalar(MHDControls.subDict("sigmaDoping").lookup("dopingRadius")));
scalar dopingLength(readScalar(MHDControls.subDict("sigmaDoping").lookup("dopingLength")));
scalar dopingAmp(readScalar(MHDControls.subDict("sigmaDoping").lookup("dopingAmp")));

//anode boundary condition parameters
bool   fixedAnodeConductivity(readBool(MHDControls.subDict("anodeBC").lookup("fixedAnodeConductivity")));
bool   electronTemperatureBoundaryConductivity(readBool(MHDControls.subDict("anodeBC").lookup("electronTemperatureBoundaryConductivity")));

scalar anodeElectronTemperature(readScalar(MHDControls.subDict("anodeBC").lookup("anodeElectronTemperature")));
scalar anodeSigmaTe(readScalar(MHDControls.subDict("anodeBC").lookup("anodeSigmaTe")));
scalar anodeMaxBLSigma(readScalar(MHDControls.subDict("anodeBC").lookup("anodeMaxBLSigma")));
scalar anodeWorkFunction(readScalar(MHDControls.subDict("anodeBC").lookup("anodeWorkFunction")));//in eV!

scalar BLThickness(readScalar(MHDControls.subDict("anodeBC").lookup("BLThickness")));
scalar BLTransitionLength(readScalar(MHDControls.subDict("anodeBC").lookup("BLTransitionLength")));
scalar anodeWallTemperature(readScalar(MHDControls.subDict("anodeBC").lookup("anodeWallTemperature")));
scalar totalCurrent(readScalar(MHDControls.subDict("anodeBC").lookup("totalCurrent")));

scalar anodeHeat = totalCurrent*(anodeWorkFunction);//+1.723e-4*anodeElectronTemperature);
dimensionedScalar twokBOvere("twokBOvere", dimensionSet(1, 2, -3, -1, 0, -1, 0), 1.723e-4);
//[1 -1 -3 0 0 0 0] - [0 -3 0 1 0 1 0]



//cathode boundary condition parameters
scalar spotRadius(readScalar(MHDControls.subDict("cathodeBC").lookup("spotRadius")));
scalar heatToCathode(readScalar(MHDControls.subDict("cathodeBC").lookup("heatToCathode")));
scalar sigmaCath(readScalar(MHDControls.subDict("cathodeBC").lookup("sigmaCath")));
vector cathodePoint(vector(MHDControls.subDict("cathodeBC").lookup("cathodePoint")));
scalar cathodeWallTemperature(readScalar(MHDControls.subDict("cathodeBC").lookup("cathodeWallTemperature")));
scalar cathShapePower = 6;
scalar cathShapeFactor = .75;
scalar cathodeHeatFlux = heatToCathode/(3.14*spotRadius*spotRadius)/cathShapeFactor;

scalar spongeFrictionCoefficient(readScalar(MHDControls.subDict("spongeLayer").lookup("spongeFrictionCoefficient")));

scalar spongeOnset(readScalar(MHDControls.subDict("spongeLayer").lookup("spongeOnset")));
scalar spongeLength(readScalar(MHDControls.subDict("spongeLayer").lookup("spongeLength")));

vector spongeVelocity(vector(MHDControls.subDict("spongeLayer").lookup("spongeVelocity")));

Info << "MHDControls are \n"
	
	<<"nBcorr"<<"  "<<nBcorr<<" "
	<<"EM_scheme"<<"  "<<EM_scheme<<" "
	<<"J_scheme"<<"  "<<J_scheme<<" "
	<<"lorenzForceOn"<<"  "<<lorenzForceOn<<" "
	<<"dopingOn"<<"  "<<dopingOn<<" "
	<<"xDope"<<"  "<<xDope<<" "
	<<"yDope"<<"  "<<yDope<<" "
	<<"zDope"<<"  "<<zDope<<" "
	<<"dopingRadius"<<"  "<<dopingRadius<<" "
	<<"dopingLength"<<"  "<<dopingLength<<" "
	<<"dopingAmp"<<"  "<<dopingAmp<<" "
	<<"anodeElectronTemperature"<<"  "<<anodeElectronTemperature<<" "
	<<"fixedAnodeConductivity"<<"  "<<fixedAnodeConductivity<<" "
	<<"electronTemperatureBoundaryConductivity"<<"  "<<electronTemperatureBoundaryConductivity<<" "
	<<"anodeSigmaTe"<<"  "<<anodeSigmaTe<<" "
	<<"anodeMaxBLSigma"<<"  "<<anodeMaxBLSigma<<" "
	<<"anodeWorkFunction"<<"  "<<anodeWorkFunction<<" "
	

	<<"BLThickness"<<"  "<<BLThickness<<" "
	<<"BLTransitionLength"<<"  "<<BLTransitionLength<<" "
	<<"anodeWallTemperature"<<"  "<<anodeWallTemperature<<" "
	<<"totalCurrent"<<"  "<<totalCurrent<<" "
	<<"spotRadius"<<"  "<<spotRadius<<" "
	<<"heatToCathode"<<"  "<<heatToCathode<<" "
	<<"cathShapePower"<<"  "<<cathShapePower<<" "
	<<"cathodeWallTemperature"<<"  "<<cathodeWallTemperature<<" "
	
	<<"sigmaCath"<<"  "<<sigmaCath<<" "
	
	<<"cathodeHeatFlux"<<"  "<<cathodeHeatFlux<<" "
	
	
	<<"\n";
	
