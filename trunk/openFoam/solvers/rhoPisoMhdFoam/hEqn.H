{
    //Info<< "min h is " << Hmin << endl;
    
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, H)
      + fvm::div(phi, H)
      //- fvm::laplacian(turbulence->alphaEff(), H)
      - fvm::laplacian(kappa/C, H)
      - fvm::laplacian((kappa_eddy)/C, H)
     ==
        jouleHeating
		- Qrad
		+ QAnode
		+ QCathode
		
		+fvc::laplacian(kappa, Temperature)
		+fvc::laplacian(kappa_eddy, Temperature)
		-fvc::laplacian(kappa/C,H)
		-fvc::laplacian(kappa_eddy/C,H)
	
	//	-fvc::laplacian(twokBOvere*SIGMA*Te, phee)
     // + DpDt
    );
    
    hEqn.solve();
    
    H.correctBoundaryConditions();
    
    fvScalarMatrix TeEqn
    (
        fvm::ddt(rho, Te)
      + fvm::div(phi, Te)
      //- fvm::div(phiJ*twokBOvere/C_e, Te)
      //- fvm::laplacian(turbulence->alphaEff(), H)
      - fvm::laplacian(kappa_e/C_e, Te)
      - fvm::laplacian((kappa_eddy)/C_e, Te)
     ==
        jouleHeating/C_e
        -Temperature*(Te-Temperature)/(Te_relaxationConstant)		
		//+fvc::laplacian(kappa_e/C_e, Temperature)
		//+fvc::laplacian(kappa_eddy/C_e, Temperature)
		//-fvc::laplacian(kappa_e/C_e,Te)
		//-fvc::laplacian(kappa_eddy/C_e,Te)
	
		-fvc::laplacian(twokBOvere*SIGMA*Te/C_e, phee)
    );
   // TeEqn.relax();
    TeEqn.solve();
    
	#include "updateProperties.H"
	#include "updateDoping.H"
	#include "updateElectrodeBCs.H"
	/*
	dT = Temperature-T;
	for (int i =0; i<10;i++)
	{
		h = (T+i*dT/10)*thermo.Cp();
		thermo.correct();
		//dT = Temperature-T;
	}*/
	h = Temperature*thermo.Cp();
	T.internalField() = Temperature.internalField();
	thermo.correct();
	
	
}
