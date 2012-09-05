{
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
	
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rAU));

    U = rAU*UEqn.H();
    U.correctBoundaryConditions();

    phi = (fvc::interpolate(U) & mesh.Sf())
        + fvc::ddtPhiCorr(rAU, U, phi);

    surfaceScalarField buoyancyPhi(rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());
    phi -= buoyancyPhi;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phi)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi -= p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U -= rAU*fvc::reconstruct((buoyancyPhi + p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            
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
        }
    }

    #include "continuityErrs.H"

    p = p_rgh + rhok*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rhok*gh;
    }
}
