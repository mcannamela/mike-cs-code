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
    
    /*autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();*/
    
    linear_enthalpy_LUT_gas plasma(gasDir);
    Info<< "plasma LUT constructed\n" << endl;
    linear_enthalpy_LUT_gas ambientGas(ambientGasDir);
    Info<< "ambient gas LUT constructed\n" << endl;

    /*volScalarField& p = thermo.p();
    volScalarField& h = thermo.h();
    const volScalarField& psi = thermo.psi();
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
    );*/
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
   /* volScalarField rhoPlasma
    (
        IOobject
        (
            "rhoPlasma",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1,-3,0, 0, 0, 0, 0)
    );
    volScalarField rhoAmbient
    (
        IOobject
        (
            "rhoAmbient",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1,-3,0, 0, 0, 0, 0)
    );*/
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
    double Tmax = 0;
    
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
        T
    );*/
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

    dimensionedScalar rhoMax
    (
        mesh.solutionDict().subDict("PIMPLE").lookup("rhoMax")
    );

    dimensionedScalar rhoMin
    (
        mesh.solutionDict().subDict("PIMPLE").lookup("rhoMin")
    );
    //update enthalpy using temperature field
	forAll(T.internalField(), celli)
	{
		H.internalField()[celli] = plasma.h(T.internalField()[celli]);//*f.internalField()[celli]
			//+(1-f.internalField()[celli])*ambientGas.h(T.internalField()[celli]);
	}
    #include "updateThermo.H"
    Info<< "Creating field DpDt\n" << endl;
    volScalarField DpDt
    (
        fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p)
    );
    Info<<"init H\n";
    
    
	
