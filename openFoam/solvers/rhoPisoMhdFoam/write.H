if (runTime.outputTime())
{
	   	
	if( fullOutputOn)
	   	{
	    volVectorField gradSigma
	    (
		IOobject
		(
		    "gradSigma",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		fvc::grad((SIGMA))
	    );
	    
	    volVectorField gradPhee
	    (
		IOobject
		(
		    "gradPhee",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		fvc::grad(phee)
	    );
	    
	    volVectorField minusGradP
	    (
		IOobject
		(
		    "minusGradP",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		-fvc::grad(p)
	    );
	    
	    volVectorField UxB
	    (
		IOobject
		(
		    "UxB",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		(U ^ B)
	    );
	    
	    volVectorField ddtA
	    (
		IOobject
		(
		    "ddtA",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		fvc::ddt(A)
	    );
	    
	    	    
	    volScalarField pressureWork
	    (
		IOobject
		(
		    "DpDt",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		DpDt
	    );
	    
	    
	    
	    volScalarField hDiffusionResidual
	    (
		IOobject
		(
		    "diffusionResidual",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		 fvc::laplacian(kappa, Temperature)
		  -fvc::laplacian(kappa/(C),H)
	    );
	    
	    volScalarField hDiffusion
	    (
		IOobject
		(
		    "hDiffusion",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		
		  fvc::laplacian(kappa/(C),H)
	    );
	    
	    volScalarField TDiffusion
	    (
		IOobject
		(
		    "TDiffusion",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		
		  fvc::laplacian(kappa, Temperature)
	    );
	    volScalarField turbTDiffusion
	    (
		IOobject
		(
		    "turbTDiffusion",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		
		  fvc::laplacian(kappa_eddy, Temperature)
	    );
	    volScalarField TeTransport
	    (
		IOobject
		(
		    "TeTransport",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		
		  fvc::laplacian(twokBOvere*SIGMA*Te, phee)
	    );
	    
	    
	    volScalarField d_rhoH_dt
	    (
		IOobject
		(
		    "d_rhoH_dt",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		fvc::ddt(rho, H)
	    );
	    Info << "writing now..." << endl;
        runTime.write();
	    }
	    else
	    {
			Info << "writing now..." << endl;
        runTime.write();
		}
		
}
