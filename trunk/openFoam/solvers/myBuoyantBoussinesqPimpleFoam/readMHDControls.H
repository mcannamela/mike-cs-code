#include <string>
const dictionary& IFDict = mesh.solutionDict().subDict("MHDControls");
Foam::wordList inlets = IFDict.subDict("inlets").toc();
scalar tanTheta(readScalar(IFDict.lookup("tanTheta")));

//vector cathodePoint(vector(MHDControls.subDict("cathodeBC").lookup("cathodePoint")));


Info << "inletForce variables are \n"
	
	<<"inlets"<<"  "<<inlets<<" "
	<<"tanTheta"<<"  "<<tanTheta<<" "
	
	<<"\n";
	
