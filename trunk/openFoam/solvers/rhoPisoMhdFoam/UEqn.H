    
	double blend = 0;
	
	
	forAll(spongeFriction.internalField(), celli)
	{
		cc = mesh.C().internalField()[celli];
			
		arg = 2.0*double(cc.component(0)-spongeOnset)/spongeLength;
		arg = std::max(std::min(arg, 20.0), -20.0);
		//Info << "arg "<<"is "		<<arg		<<"\n";
		num = (std::exp(arg)-1);
		//Info << "num is "		<<"is "		<<num		<<"\n";
		denom = (1e-10+std::exp(arg)+1);
		//Info << "denom is "		<<"is "				<<denom		<<"\n";
		
		blend = 1-(.5+.5*num/denom);
			
		spongeFriction.internalField()[celli] = blend*spongeFrictionCoefficient*(U.internalField()[celli]-spongeVelocity*anodeDistance.internalField()[celli]*cathodeDistance.internalField()[celli]/(1*1));
	}
    
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      //- fvc::div(phiB, 2.0*DBU*B)
      //+ fvc::grad(DBU*magSqr(B))
      
      //+ turbulence->divDevRhoReff(U)
      - fvm::laplacian(MU+MU_eddy, U)
      - (lorenzForce)
      //- (100*MU*fvc::grad(fvc::div(surfaceScalarField("phiU", phi/fvc::interpolate(rho)))))
     // ==
      //-spongeFriction
    );

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }
