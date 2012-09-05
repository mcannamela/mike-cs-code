    Info<< "Reading thermophysical properties\n" << endl;
   
    Info<< "Reading field kappat\n" << endl;
    volScalarField kappat
    (
        IOobject
        (
            "kappat",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field T\n" << endl;
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

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

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
    
    /*volVectorField inletVelocity
    (
        IOobject
        (
            "inletVelocity",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );*/
    volVectorField inletForce
    (
        IOobject
        (
            "inletForce",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    inletForce*=0;
    
    labelUList faceCells;
    vectorField  nf;
    vector forceDir;
    
    forAll(inlets, ini)
    {
		label patchi = mesh.boundary().findPatchID(inlets[ini]);
		faceCells = mesh.boundary()[inlets[ini]].faceCells();
		const scalarField faceArea = mesh.boundary()[inlets[ini]].magSf();
		nf = mesh.boundary()[inlets[ini]].nf();
		forceDir = vector(IFDict.subDict("inlets").lookup(inlets[ini]));
		forAll(inletForce.boundaryField()[patchi], facei)
		{
			inletForce.internalField()[faceCells[facei]] = forceDir*(sqr(U.boundaryField()[patchi][facei] & nf[facei])*faceArea[facei]/(mesh.V()[faceCells[facei]]*tanTheta) );
		}
	}
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

    #include "createPhi.H"

    #include "readTransportProperties.H"

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    // Kinematic density for buoyancy force
    volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh
        ),
        1.0 - beta*(T - TRef)
    );

    // kinematic turbulent thermal thermal conductivity m2/s


    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rhok*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }

