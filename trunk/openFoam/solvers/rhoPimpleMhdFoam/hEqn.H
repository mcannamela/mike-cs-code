
    
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
	
		//-fvc::laplacian(twokBOvere*SIGMA*Temperature, phee)
    
     //   - (fvc::ddt(rho, K) + fvc::div(phi, K))
     // + dpdt
    );

   // hEqn.relax();
    hEqn.solve();
    
    H.correctBoundaryConditions();
    
    fvScalarMatrix TeEqn
    (
        fvm::ddt(rho, Te)
      + fvm::div(phi, Te)
      - fvm::laplacian(kappa_e/C_e, Te)
      - fvm::laplacian((kappa_eddy)/C_e, Te)
     ==
        jouleHeating/C_e
        -Temperature*(Te-Temperature)/(Te_relaxationConstant)	
        //+QAnode/C_e	
		//-fvc::laplacian(twokBOvere*SIGMA*Te/C_e, phee)
    );
   // TeEqn.relax();
    TeEqn.solve();
    
	#include "updateProperties.H"
	#include "updateDoping.H"
	#include "updateElectrodeBCs.H"
	
	PSI = 1.0/(gasConstant*Temperature);
	RHO = p*PSI;
	phiu = (fvc::interpolate(U) & mesh.Sf());
	phi = fvc::interpolate(RHO)*(fvc::interpolate(U) & mesh.Sf());
	
	h = Temperature*thermo.Cp();
	//T.internalField() = Temperature.internalField();
	
	TDiffusion = fvc::laplacian(kappa, Temperature);
	turbTDiffusion = fvc::laplacian(kappa_eddy, Temperature);
	

	thermo.correct();
	
	
   

