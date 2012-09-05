    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();
    
    linear_enthalpy_LUT_gas plasma(gasDir);//"/home/wichtelwesen/air/");
    Info<< "plasma LUT constructed\n" << endl;

    volScalarField& p = thermo.p();
    volScalarField& h = thermo.h();
    const volScalarField& psi = thermo.psi();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );
    volScalarField RHO
    (
        IOobject
        (
            "RHO",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1,-3,0,0,0,0,0)
    );
    double massSurplus = 0;
    double initialMass = 0;
    double surplus = 0;
    
    volScalarField PSI
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
    );
    
    dimensionedScalar gasConstant("gasConstant", dimensionSet(0,2,-2,-1,0,0,0), 8314.0/gasMolarMass); 

    Info<< "Reading field U\n" << endl;
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
   
    //Te.internalField() = Temperature.internalField();
    
    volScalarField dT
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
    );
    

    //#include "compressibleCreatePhi.H"
    surfaceScalarField phi
 (
     IOobject
     (
         "phi",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
    ),
     linearInterpolate(rho*U) & mesh.Sf()
 );
    
    
    volScalarField divA
    (
        IOobject
        (
            "divA",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, 0, -2, 0, 0, -1, 0)
    );
    Info<< "Reading field A\n" << endl;
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
    divA= fvc::div(A);
    
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
	
	volVectorField tipForce
		(
		    IOobject
		    (
		        "tipForce",
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    mesh,
		    dimensionSet(1,-2,-2,0,0,0,0)
		);
	tipForce*=0;
	volVectorField outletForce
		(
		    IOobject
		    (
		        "outletForce",
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    mesh,
		    dimensionSet(1,-2,-2,0,0,0,0)
		);
	outletForce*=0;
	volVectorField nookForce
		(
		    IOobject
		    (
		        "nookForce",
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    mesh,
		    dimensionSet(1,-2,-2,0,0,0,0)
		);
	nookForce*=0;
	

		
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
    
    phiJ *=0;// fvc::interpolate((SIGMA))*fvc::snGrad(phee)*mesh.magSf();
    
    /////////////////////////////////////////////////////////////////////
    //////////////////////initialize H, T
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
	//T.internalField() = Temperature.internalField();
	thermo.correct();
	
	PSI = 1.0/(gasConstant*Temperature);
	RHO = p*PSI;
    /////////////////////////////////////////////////////////////////////
    
    //make sure all properties are valid
    #include "updateProperties.H"
    
    /*
    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    //\\\\\\\\\\\\\\\\\initialize current density\\\\\\\\\\\\\\\\\\\\\
    E = -fvc::grad(phee);
	E.correctBoundaryConditions();
	
	J = (SIGMA)*E;
	J.correctBoundaryConditions();
	
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
			jouleHeating.internalField()[celli] = mag(J.internalField()[celli])*
				(mag(J.internalField()[celli])/SIGMA.internalField()[celli]);
		}
	}
	jouleHeating = min(jouleHeating, cathodeJouleHeatingMax);
    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    */

    //clamps on density
    dimensionedScalar rhoMax(pimple.dict().lookup("rhoMax"));
    dimensionedScalar rhoMin(pimple.dict().lookup("rhoMin"));

    //initialize turbulence model
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
    
    
    //sundry fields left over, the pressure work and the kinetic energy
    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt("dpdt", fvc::ddt(p));

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));
    
    surfaceScalarField phiu("phiu", (fvc::interpolate(U) & mesh.Sf()));
	
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
	
	
