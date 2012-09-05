{
	D = (MU+MU_eddy)/(.7*rho);
	
	fvScalarMatrix fEqn
    (
        fvm::ddt(rho, f)
      + fvm::div(phi, f)
      ==
      fvc::laplacian(D*rho, f)
      //- fvm::laplacian(.7*MU, f)
            
      //  fvc::laplacian((D*f+ dimensionedScalar("fSmall", dimensionSet(0,2,-1,0,0,0,0), 1e-6)),rho) 
      
    );
    fEqn.solve();
    
    
    /*fvScalarMatrix pfEqn
    (
        fvm::ddt(rhoPlasma)
        ==
      - fvc::div(phi, f)
      + fvm::laplacian(D, rho*f)
     
    );
    pfEqn.solve();
    
    fvScalarMatrix afEqn
    (
        fvm::ddt(rhoAmbient)
        ==
      - fvc::div(phi, (1-f))
      + fvm::laplacian(D, rho*(1-f))
     
    );
    afEqn.solve();
    
    f = rhoPlasma/(rhoPlasma+rhoAmbient)
    rhoPlasma = rho*f;
    rhoAmbient = rho*(1-f);
    */
    /*f = max(f, dimensionedScalar("pfMin", dimensionSet(0,0,0,0,0,0,0), 0.0));
	f = min(f, dimensionedScalar("pfMax", dimensionSet(0,0,0,0,0,0,0), 1.0));
	forAll(f.internalField(), celli){
		if (f.internalField()[celli]>1.0)
			Info << "unphysical plasma fraction detected: " << f.internalField()[celli] <<"\n";
	}*/
	
    f.correctBoundaryConditions();
}
