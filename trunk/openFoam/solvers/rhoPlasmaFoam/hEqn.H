{
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, H)
      + fvm::div(phi, H)
      //- fvm::laplacian(turbulence->alphaEff(), H)
      - fvm::laplacian((K+K_eddy)/C_p, H)
     ==
        DpDt
        + fvc::laplacian(K+K_eddy, T)
		-fvc::laplacian((K+K_eddy)/C_p,H)
		-Qrad
    );

    hEqn.relax();
    hEqn.solve();

    H.correctBoundaryConditions();
    #include "updateProperties.H"
	
	/*dT = T-T;
	for (int i =0; i<10;i++)
	{
		h = (T+i*dT/10)*thermo.Cp();
		thermo.correct();
		//dT = T-T;
	}*/
	//h = T*thermo.Cp();
	//thermo.correct();
	#include "updateThermo.H"
	Info<< "			max temperature is " << Tmax << "\n";
	DiffusionResidual = fvc::laplacian(K+K_eddy, T)-fvc::laplacian((K+K_eddy)/C_p,H);
	TDiffusion = fvc::laplacian(K+K_eddy, T);
	HDiffusion = fvc::laplacian((K+K_eddy)/C_p,H);
    
    
}
