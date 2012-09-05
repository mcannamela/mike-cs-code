if (EM_scheme>0)
{
	fvVectorMatrix AEqn
			(
			fvm::laplacian(A) 
			== 
			  - mu0*J
			  
			);
	/*fvVectorMatrix AEqn
			(
			mu0*SIGMA*fvm::ddt(A)-fvm::laplacian(A) 
			== 
			   mu0*SIGMA*fvc::grad(phee)
			  
			);*/
    AEqn.solve();
    
    snGradA = fvc::snGrad(A);
    lapA = fvc::laplacian(A);
    forAll(lapA.internalField(), celli)
    {
		lapA.internalField()[celli]*=mesh.V()[celli];
	}
    
    //Info << "A has been solved."<< endl;
    A.correctBoundaryConditions();
    //Info << "BCs on A corrected."<< endl;
    B = fvc::curl(A);
    //Info << "B computed!"<< endl;
    B.correctBoundaryConditions();
    //Info << "BCs on B corrected."<< endl;
    
    
}
else
{

	A*=0;
	B*=0;

}

if (lorenzForceOn)
lorenzForce = J ^ B;
else
lorenzForce*= 0;
