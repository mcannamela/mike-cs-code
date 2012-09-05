/////////////////////////////////////////////////////////////
////approximate the electric field with a poisson equation///
/////////////////////////////////////////////////////////////


if (EM_scheme>0)
{
	fvScalarMatrix pheeEqn0
			(
			   fvm::laplacian( phee)   
			);



	fvScalarMatrix pheeEqn1
			(
			   fvm::laplacian(SIGMA, phee) 
			);

	fvScalarMatrix pheeEqn2
			(
			   fvm::laplacian(phee) == 
			   -(1/SIGMA)*(fvc::grad(SIGMA) & fvc::grad(phee))
			);

	fvScalarMatrix pheeEqn3
			(
			   fvm::laplacian((SIGMA), phee) 
			   == 
				//- fvc::div((SIGMA)*(U ^ B) )
				- fvc::div((SIGMA)*fvc::ddt(A))
			);
		


	switch(EM_scheme)
	{
		case 0:
			pheeEqn0.solve();
			break;

		case 1:
			pheeEqn1.solve();
			break;
		case 2:
			pheeEqn2.solve();
			break;
		default:
			pheeEqn3.solve();
			break;
	}
	
	phiJ = fvc::interpolate((SIGMA))*fvc::snGrad(phee)*mesh.magSf();
	
	
	
	phee.correctBoundaryConditions();
	SIGMA.correctBoundaryConditions();
	
		
	if (J_scheme == 0)
	{
		E = -fvc::grad(phee);
		E.correctBoundaryConditions();
		
		J = (SIGMA)*E;
		J.correctBoundaryConditions();
	}
	else if (J_scheme == 1)
	{
		E = -fvc::grad(phee) - fvc::ddt(A);
		E.correctBoundaryConditions();
		
		J = (SIGMA)*(E + (U ^ B));
		J.correctBoundaryConditions();
	}

	
	magE = mag(E);
	jouleHeating = mag(J)*(mag(J)/(SIGMA));
	jouleHeating = min(jouleHeating, jouleHeatingMax);
	
	forAll(mesh.cellZones(), czi)
	{
		//Info << "re-computing joule heating for cellZone " << czi<<endl;
		const labelList& cellLabels = mesh.cellZones()[czi];
		forAll(cellLabels, cli)
		{
			label celli = cellLabels[cli];
			jouleHeating.internalField()[celli] = mag(J.internalField()[celli])*(mag(J.internalField()[celli])/SIGMA.internalField()[celli]);
			//if ((mag(J.internalField()[celli])*mag(J.internalField()[celli])/SIGMA.internalField()[celli])>5e11)
				//Info << jouleHeating.internalField()[celli] << endl;
		}
	}
	jouleHeating = min(jouleHeating, cathodeJouleHeatingMax);
	//Info << "max jouleHeating is "<< max(jouleHeating)<<endl;
	

	
	/////////////////////////////////////////////////////////////
}
else
{
	phee*=0;
	A*=0;
	B*=0;
	J*=0;
	E*=0;
	jouleHeating*=0;
}
    forAll(JdV.internalField(), celli)
    {
		JdV.internalField()[celli]=J.internalField()[celli]*mesh.V()[celli];
	}
	
if (lorenzForceOn>0)
	lorenzForce = J ^ B;
else
	lorenzForce*= 0;
