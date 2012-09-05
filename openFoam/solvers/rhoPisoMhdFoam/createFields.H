    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

    linear_enthalpy_LUT_gas plasma(gasDir);//"/home/wichtelwesen/air/");
    Info<< "plasma LUT constructed\n" << endl;
    //?????????????????????????????????????????????????????
    // what determines which of these will write? 
    // if i change the value of p, will that change thermo.p_? 
    volScalarField& p = thermo.p();
    volScalarField& h = thermo.h();
    Info<< "dim of h " << h.dimensions() << nl << endl;
    const volScalarField& psi = thermo.psi();
    //?????????????????????????????????????????????????????

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );
    double massSurplus = 0;
    double initialMass = 0;
    double surplus = 0;
    /*volScalarField PSI
    (
        IOobject
        (
            "PSI",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0,0,0,0 )
    );*/
    
    //dimensionedScalar gasConstant("gasConstant", dimensionSet(0,2,-2,-1,0,0,0), 8314.0/28.0); 
    
    volScalarField H
    (
        IOobject
        (
            "H",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    double Hmin = 1e15;
    
    volScalarField Qrad
    (
        IOobject
        (
            "Qrad",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -3, 0,0,0,0 )
    );
    volScalarField QAnode
    (
        IOobject
        (
            "QAnode",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -3, 0,0,0,0 )
    );
    volScalarField QCathode
    (
        IOobject
        (
            "QCathode",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -3, 0,0,0,0 )
    );
    Qrad*=0;
    QAnode*=0;
    QCathode*=0;
    
    volScalarField Temperature
    (
        IOobject
        (
            "Temperature",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    double lastTMax = 30000;
     volScalarField Te
    (
        IOobject
        (
            "Te",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    /*volScalarField dT
    (
        IOobject
        (
            "dT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Temperature
    );*/
    
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    //double Tmax = 0;


    Info<< "\nReading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

#   include "compressibleCreatePhi.H"

  


    Info<< "Reading field A\n" << endl;
    volVectorField lapA
    (
        IOobject
        (
            "lapA",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -2, 0 ,0, -1, 0)
    );
    surfaceVectorField snGradA
		(
		    IOobject
		    (
		        "snGradA",
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    mesh,
        dimensionSet(1, 0, -2, 0, 0, -1, 0)
    );
    
    volVectorField J
		(
		    IOobject
		    (
		        "J",
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    mesh,
		    dimensionSet(0, -2, 0, 0, 0, 1, 0 )
		);

	J*=0;	
	    volVectorField JdV
		(
		    IOobject
		    (
		        "JdV",
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    J
		);
    forAll(JdV.internalField(), celli)
    {
		JdV.internalField()[celli]*=mesh.V()[celli];
	}
    
    volVectorField A
    (
        IOobject
        (
            "A",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    snGradA = fvc::snGrad(A);
    lapA = fvc::laplacian(A);
    forAll(lapA.internalField(), celli)
    {
		lapA.internalField()[celli]*=mesh.V()[celli];
	}
	
	
    
    volVectorField B
    (
        IOobject
        (
            "B",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1 ,0, -2, 0, 0, -1, 0 )
    );
    
    

		
	volVectorField lorenzForce
		(
		    IOobject
		    (
		        "lorenzForce",
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    J ^ B
		);
	volVectorField spongeFriction
		(
		    IOobject
		    (
		        "spongeFriction",
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    mesh,
		    dimensionSet(1,-2,-2,0,0,0,0)
		);
	spongeFriction*=0;
	surfaceScalarField phiJ
		(
		    IOobject
		    (
		        "phiJ",
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    fvc::interpolate(J) & mesh.Sf()
		);
		
    volScalarField phee
    (
        IOobject
        (
            "phee",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
	volVectorField E
		(
		    IOobject
		    (
		        "E",
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    mesh,
		    dimensionSet(1, 1 ,-3, 0, 0, -1, 0)
		 );	
	volScalarField magE
		(
		    IOobject
		    (
		        "magE",
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    mesh,
		    dimensionSet(1, 1 ,-3, 0, 0, -1, 0)
		 );	
	magE*=0;
	
	volScalarField jouleHeating
    (
        IOobject
        (
            "jouleHeating",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1 ,-3, 0, 0, 0, 0)
    );
    jouleHeating*=0;

	
	forAll(Temperature.boundaryField(), patchi)
	{
		fvPatchScalarField& pT = Temperature.boundaryField()[patchi];
		fvPatchScalarField& pH = H.boundaryField()[patchi];

			forAll(pT, facei)
			{
				pH[facei] = plasma.h(pT[facei]);
			}	
	}
    //update enthalpy using temperature field
	forAll(Temperature.internalField(), celli)
	{
		H.internalField()[celli] = plasma.h(Temperature.internalField()[celli]);
	}
	h = Temperature*thermo.Cp();
	thermo.correct();
	
	#include "updateProperties.H"
	
	Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Creating field DpDt\n" << endl;
    volScalarField DpDt =
        fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
	Info<< "done\n" << endl;
