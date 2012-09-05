//rho = thermo.rho();
#include "updateThermo.H"
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
rhot.internalField() = rho.internalField();

U = rAU*UEqn().H();

/*if (pimple.nCorr() <= 1)
{
    UEqn.clear();
}*/
if (pimple.nCorrPISO() <= 1)
{
    UEqn.clear();
}


if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, rho, U, phi)
        )
    );

    for (int nonOrth=0; nonOrth<=pimple.nNonOrthCorr(); nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho*rAU, p)
        );

        /*pEqn.solve
        (
            mesh.solver(p.select(pimple.finalInnerIter(corr, nonOrth)))
        );*/
        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (nonOrth == pimple.nNonOrthCorr())
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    phi =
        fvc::interpolate(rho)*
        (
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, rho, U, phi)
        );

    for (int nonOrth=0; nonOrth<=pimple.nNonOrthCorr(); nonOrth++)
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho*rAU, p)
        );

        /*pEqn.solve
        (
            mesh.solver(p.select(pimple.finalInnerIter(corr, nonOrth)))
        );*/
        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
        p = max(p, pMin);
		p = min(p, pMax);

        if (nonOrth == pimple.nNonOrthCorr())
        {
            phi += pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
//#include "compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

Info << "p min, max: "<< min(p).value()
	<<" " << max(p).value() << "\n";
// Recalculate density from the relaxed pressure
//rho = thermo.rho();
p = max(p, pMin);
p = min(p, pMax);

p.correctBoundaryConditions();
pt.internalField() = p.internalField();

#include "updateThermo.H"
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
rhot.internalField() = rho.internalField();
Info<< "rho min, max : " << min(rho).value()
    << " " << max(rho).value() << endl;
   


U -= rAU*fvc::grad(p);
U.correctBoundaryConditions();

DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
