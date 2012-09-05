{
    double blend = 0;
    double arg, num, denom;
    vector cc;
	forAll(spongeCooling.internalField(), celli)
	{
		cc = mesh.C().internalField()[celli];
			
		arg = 2.0*double(cc.component(0)-spongeOnset)/spongeLength;
		arg = std::max(std::min(arg, 20.0), -20.0);
		//Info << "arg "<<"is "		<<arg		<<"\n";
		num = (std::exp(arg)-1);
		//Info << "num is "		<<"is "		<<num		<<"\n";
		denom = (1e-10+std::exp(arg)+1);
		//Info << "denom is "		<<"is "				<<denom		<<"\n";
		
		blend = (.5+.5*num/denom);
			
		spongeCooling.internalField()[celli] = blend*spongeCoolingCoefficient*(T.internalField()[celli]-spongeTemperature);
	}
	
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, H)
      + fvm::div(phi, H)
      //- fvm::laplacian(turbulence->alphaEff(), H)
      - fvm::laplacian((K+K_eddy)/C_p, H)
     ==
        //DpDt
        fvc::laplacian(K+K_eddy, Temperature)
		-fvc::laplacian((K+K_eddy)/C_p,H)
		-Qrad
	//	-spongeCooling
    );

    hEqn.relax();
    hEqn.solve();

    #include "updateProperties.H"
    H.correctBoundaryConditions();
    
	/*
	dT = T-Temperature;	
	for (int i =0; i<10;i++)
	{
		h = (T+i*dT/10)*thermo.Cp();
		thermo.correct();
		dT = T-Temperature;
	}
	h = Temperature*thermo.Cp();
	thermo.correct();*/
	T.internalField() = Temperature.internalField();
	#include "updateThermo.H"
	Info<< "			max temperature is " << Tmax << "\n";
	DiffusionResidual = fvc::laplacian(K+K_eddy, Temperature)-fvc::laplacian((K+K_eddy)/C_p,H);
	TDiffusion = fvc::laplacian(K+K_eddy, Temperature);
	HDiffusion = fvc::laplacian((K+K_eddy)/C_p,H);
    
    
}
