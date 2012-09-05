    // Solve the momentum equation
    /*forAll(inletVelocity.boundaryField(), patchi)
	{
		labelUList faceCells;

		if (mesh.boundary()[patchi].name() == "inlet1")
		{
			//Info<< "Computing anode heat flux" << endl;
			faceCells = mesh.boundary()[patchi].faceCells();
			vector inletDir(0, -.174, -.985);
			forAll(inletVelocity.boundaryField()[patchi], facei)
			{
				inletVelocity.internalField()[faceCells[facei]] = mag(U.boundaryField()[patchi][facei])*inletDir;	
			}
			
		}
	
		if (mesh.boundary()[patchi].name() == "inlet2")
		{
			//Info<< "Computing anode heat flux" << endl;
			faceCells = mesh.boundary()[patchi].faceCells();
			vector inletDir(.985, -.174, 0);
			forAll(inletVelocity.boundaryField()[patchi], facei)
			{
				inletVelocity.internalField()[faceCells[facei]] = mag(U.boundaryField()[patchi][facei])*inletDir;	
			}
			
		}
	
		if (mesh.boundary()[patchi].name() == "inlet3")
		{
			//Info<< "Computing anode heat flux" << endl;
			faceCells = mesh.boundary()[patchi].faceCells();
			vector inletDir(-.985, -.174, 0);
			forAll(inletVelocity.boundaryField()[patchi], facei)
			{
				inletVelocity.internalField()[faceCells[facei]] = mag(U.boundaryField()[patchi][facei])*inletDir;	
			}
			
		}
	
		if (mesh.boundary()[patchi].name() == "inlet4")
		{
			//Info<< "Computing anode heat flux" << endl;
			faceCells = mesh.boundary()[patchi].faceCells();
			vector inletDir(0, -.174, .985);
			forAll(inletVelocity.boundaryField()[patchi], facei)
			{
				inletVelocity.internalField()[faceCells[facei]] = mag(U.boundaryField()[patchi][facei])*inletDir;	
				
				
			}
		}
	
	}*/
    
    U.correctBoundaryConditions();
    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff(U)
      -inletForce
    );

    UEqn.relax();
    
    /*forAll(U.boundaryField(), patchi)
	{
		if (mesh.boundary()[patchi].name() == "inlet1" || 
			mesh.boundary()[patchi].name() == "inlet2" || 
			mesh.boundary()[patchi].name() == "inlet3" || 
			mesh.boundary()[patchi].name() == "inlet4" )
		{		
			UEqn.setValues(mesh.boundary()[patchi].faceCells(), inletVelocity.boundaryField()[patchi].patchInternalField());
		}
	}*/

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );
    }
    U.correctBoundaryConditions();
