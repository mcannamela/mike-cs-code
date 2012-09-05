rho = thermo.rho();

volScalarField rUA = 1.0/UEqn.A();
U = rUA*UEqn.H();

if (transonic)
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        )
    );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho*rUA, p)
        );

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    //Info<< "estimating phi\n";
    phi =
        fvc::interpolate(rho)*
        (
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        );
        
	//Info<< "Done, will now correct p\n";
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho*rUA, p)
          //+ fvc::laplacian(100*MU*rUA*rho, fvc::div(surfaceScalarField("phiU", phi/fvc::interpolate(rho))))
          
        );

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
        }
    }
}

//formerly included rhoEqn.H here 
{
     solve(fvm::ddt(rho) + fvc::div(phi));
}
Info << "p min, max: "<< min(p).value()
	<<" " << max(p).value() << "\n";
Info<< "rho max/min : " << max(rho).value()
    << " " << min(rho).value() << endl;
p = max(p, pMin);
p = min(p, pMax);

p.correctBoundaryConditions();


#include "compressibleContinuityErrs.H"

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();

DpDt = fvc::DDt( surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
