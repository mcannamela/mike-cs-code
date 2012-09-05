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
scalar smagorinskyConstant(readScalar(MHDControls.lookup("smagorinskyConstant")));

bool radiationOn(readBool(MHDControls.lookup("radiationOn")));
bool mixingOn(readBool(MHDControls.lookup("mixingOn")));
bool turbulenceOn(readBool(MHDControls.lookup("turbulenceOn")));
string gasDir(MHDControls.lookup("gasDir"));
string ambientGasDir(MHDControls.lookup("ambientGasDir"));
