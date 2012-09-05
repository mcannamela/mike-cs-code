    Info<< "Reading thermophysical properties\n" << endl;
	volScalarField C_p
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, -1, 0, 0, 0)
    );
    
    
    volScalarField K
    (
        IOobject
        (
            "K",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, 1, -3, -1, 0, 0, 0)
    );
    volScalarField K_eddy
    (
        IOobject
        (
            "K_eddy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, 1, -3, -1, 0, 0, 0)
    );
    K_eddy*=0;
    
    volScalarField MU
    (
        IOobject
        (
            "MU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0, 0, 0)
    );
    volScalarField MU_eddy
    (
    IOobject
        (
            "MU_eddy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0, 0, 0)
    );
    MU_eddy*=0;
    volScalarField D
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -1, 0, 0, 0, 0)
    );
    
    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();
    
    linear_enthalpy_LUT_gas plasma(gasDir);
    Info<< "plasma LUT constructed\n" << endl;
    linear_enthalpy_LUT_gas ambientGas(ambientGasDir);
    Info<< "ambient gas LUT constructed\n" << endl;
    
    linear_two_D_LUT gasBlendTemperature(gasBlendTemperatureFile);
    Info<< "gas blend temperature LUT constructed\n" << endl;

    volScalarField& pt = thermo.p();
    volScalarField& h = thermo.h();
    //const volScalarField& psit = thermo.psi();
    volScalarField rhot
    (
        IOobject
        (
            "rhot",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
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
        mesh,
        dimensionSet(1,-3,0, 0, 0, 0, 0)
    );
   
    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0,-2,2, 0, 0, 0, 0)
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
    double Htmp = 0;
    
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
    volScalarField spongeCooling
    (
        IOobject
        (
            "spongeCooling",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -3, 0,0,0,0 )
    );
    
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
    double Tmax = 0;
    double Ttmp = 0;
    
    volScalarField DiffusionResidual
    (
        IOobject
        (
            "DiffusionResidual",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -3, 0,0,0,0 )
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
        mesh,
        dimensionSet(1, -1, -3, 0,0,0,0 )
    );
    volScalarField HDiffusion
    (
        IOobject
        (
            "HDiffusion",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -3, 0,0,0,0 )
    );

    volScalarField f
    (
        IOobject
        (
            "f",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    forAll(f.internalField(), celli)
    {
		if (f.internalField()[celli]>1.0)
			f.internalField()[celli] = 1.0;
	}

    #include "compressibleCreatePhi.H"

    dimensionedScalar rhoMax( mesh.solutionDict().subDict("PIMPLE").lookup("rhoMax"));
    dimensionedScalar rhoMin( mesh.solutionDict().subDict("PIMPLE").lookup("rhoMin"));
    
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
    
    
    //update enthalpy using temperature field
	forAll(T.internalField(), celli)
	{
		T.internalField()[celli] = std::max(T.internalField()[celli], 250.0);
		H.internalField()[celli] = plasma.h(Temperature.internalField()[celli])*f.internalField()[celli]
			+(1-f.internalField()[celli])*ambientGas.h(Temperature.internalField()[celli]);
	}
    #include "updateThermo.H"
    Info<< "Creating field DpDt\n" << endl;
    volScalarField DpDt
    (
        IOobject
        (
            "DpDt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p)
    );
    
    /*volScalarField DpDt
    (
        fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p)
    );*/
    Info<<"init H\n";
    
    
	
