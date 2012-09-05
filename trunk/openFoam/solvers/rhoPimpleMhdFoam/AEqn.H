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
divA= fvc::div(A);
if (lorenzForceOn)
lorenzForce = J ^ B;
else
lorenzForce*= 0;
