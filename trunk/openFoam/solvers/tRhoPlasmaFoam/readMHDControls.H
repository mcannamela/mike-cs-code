#include <string>
const dictionary& MHDControls = mesh.solutionDict().subDict("MHDControls");

dimensionedScalar pMax
	(
		MHDControls.lookup("pMax")
	);

dimensionedScalar pMin
	(
		MHDControls.lookup("pMin")
	);

scalar plasmaMolarMass(readScalar(MHDControls.lookup("plasmaMolarMass")));
scalar ambientMolarMass(readScalar(MHDControls.lookup("ambientMolarMass")));


bool radiationOn(readBool(MHDControls.lookup("radiationOn")));
bool mixingOn(readBool(MHDControls.lookup("mixingOn")));
bool turbulenceOn(readBool(MHDControls.lookup("turbulenceOn")));
string gasDir(MHDControls.lookup("gasDir"));
string ambientGasDir(MHDControls.lookup("ambientGasDir"));
string gasBlendTemperatureFile(MHDControls.lookup("gasBlendTemperatureFile"));


scalar spongeFrictionCoefficient(readScalar(MHDControls.subDict("spongeLayer").lookup("spongeFrictionCoefficient")));
scalar spongeCoolingCoefficient(readScalar(MHDControls.subDict("spongeLayer").lookup("spongeFrictionCoefficient")));

scalar spongeOnset(readScalar(MHDControls.subDict("spongeLayer").lookup("spongeOnset")));
scalar spongeLength(readScalar(MHDControls.subDict("spongeLayer").lookup("spongeLength")));

scalar spongeTemperature(readScalar(MHDControls.subDict("spongeLayer").lookup("spongeTemperature")));
vector spongeVelocity(vector(MHDControls.subDict("spongeLayer").lookup("spongeVelocity")));
