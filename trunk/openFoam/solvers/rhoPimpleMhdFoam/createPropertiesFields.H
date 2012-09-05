    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    

    dimensionedScalar mu0("mu0", dimensionSet(1, 1, -2, 0, 0, -2, 0), 1.2566e-6);

    volScalarField C
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
    dimensionedScalar C_e("C_e", dimensionSet(0, 2 , -2, -1, 0, 0, 0), 416.0);
    
    volScalarField kappa
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, 1, -3, -1, 0, 0, 0)
    );
    
    dimensionedScalar kappa_e("kappa_e", dimensionSet(1, 1, -3, -1, 0, 0, 0), 4.0);
    
    dimensionedScalar Te_relaxationConstant("Te_relaxationConstant", dimensionSet(-1, 3, 1, 1, 0, 0, 0), 5000*3e-5);
    
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
    volScalarField kappa_eddy
    (
        IOobject
        (
            "kappa_eddy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, 1, -3, -1, 0, 0, 0)
    );
    kappa_eddy*=0;
    kappa_eddy+=dimensionedScalar("dummyKappa", dimensionSet(1, 1, -3, -1, 0, 0, 0), 1.0);
    
    
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

    volScalarField SIGMA
    (
        IOobject
        (
            "SIGMA",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volScalarField SIGMA_DOPING
    (
        IOobject
        (
            "SIGMA_DOPING",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(-1, -3, 3, 0, 0, 2, 0)
    );
    volScalarField SIGMA_Te
    (
        IOobject
        (
            "SIGMA_Te",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(-1, -3, 3, 0, 0, 2, 0)
    );
    volScalarField SIGMA_BL
    (
        IOobject
        (
            "SIGMA_BL",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(-1, -3, 3, 0, 0, 2, 0)
    );
    
   
    
    volScalarField betaSigma
    (
        IOobject
        (
            "betaSigma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, 0, 0, 0, 0, 0)
    );
    
    volScalarField cathodeDistance
    (
        IOobject
        (
            "cathodeDistance",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volScalarField anodeDistance
    (
        IOobject
        (
            "anodeDistance",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volVectorField anodeDirection
    (
        IOobject
        (
            "anodeDirection",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(anodeDistance)
    );
    anodeDirection/=(1e-15+mag(anodeDirection));
    
    
    betaSigma*=0;
    SIGMA_DOPING*=0;
	SIGMA_BL*=0;
	//anodeDistance.internalField() += mesh.C().component(1);
	
	forAll(mesh.boundary(), patchi){
		forAll(mesh.boundary()[patchi], facei)
			{
			SIGMA_DOPING.boundaryField()[patchi][facei] = 0;//SIGMA.boundaryField()[patchi][facei];
			
			if (mesh.boundary()[patchi].name() == "anode" || mesh.boundary()[patchi].name() == "cathode")
				SIGMA_Te.boundaryField()[patchi][facei] = 15000;
			else
				SIGMA_Te.boundaryField()[patchi][facei] = 1e-1;
			
			betaSigma.boundaryField()[patchi][facei] = 0;
			}
	}
	
   	vector cc;
   	double num, denom, arg;
   	bool takesDope;
   	
   	Info << "prepping doping conductivity and conductivity blending function"<<"\n";
	forAll(SIGMA_DOPING.internalField(), celli)
	{
		
		cc = mesh.C().internalField()[celli];
		
		arg = 2.0*double(anodeDistance.internalField()[celli] - BLThickness)/BLTransitionLength;
		arg = std::max(std::min(arg, 20.0), -20.0);
		//Info << "arg "<<"is "		<<arg		<<"\n";
		num = (std::exp(arg)-1);
		//Info << "num is "		<<"is "		<<num		<<"\n";
		denom = (1e-10+std::exp(arg)+1);
		//Info << "denom is "		<<"is "				<<denom		<<"\n";
		
		betaSigma.internalField()[celli] = .5+.5*num/denom;
		                                      
		//Info << "blend at "		<<cc		<<"is "		<<betaSigma.internalField()[celli]		<<"\n";
		
		if (dopingTwoSided)
			takesDope = (pow(cc.component(0)-xDope, 2)+pow(cc.component(1)-yDope, 2))<pow(3*dopingRadius,2);
		else
			takesDope = (pow(cc.component(0)-xDope, 2)+pow(cc.component(1)-yDope, 2))<pow(3*dopingRadius,2);
			
		if (takesDope)
		{
			if (dopingTwoSided)
			{
				/*SIGMA_DOPING.internalField()[celli] = dopingAmp*std::exp(double(-pow(.5*(cc.component(0)-xDope)/dopingRadius,2)))*
									   std::exp(double(-pow(.5*(cc.component(1)-yDope)/dopingRadius,2)))*
									 ( std::exp(double(-pow(.5*(cc.component(2)-zDope)/dopingLength,2)))+
									   std::exp(double(-pow(.5*(cc.component(2)+zDope)/dopingLength,2))) );*/
									   SIGMA_DOPING.internalField()[celli] = dopingAmp*std::exp(double(-pow(.5*(cc.component(0)-xDope)/dopingRadius,2)))*
									   (std::exp(double(-pow(.5*(cc.component(2)-zDope)/dopingLength,2)))+std::exp(double(-pow(.5*(cc.component(2)+zDope)/dopingLength,2))))*
									   std::exp(double(-pow(.5*(cc.component(1)-yDope)/dopingRadius,2)));
            }
		    else
		    {
				SIGMA_DOPING.internalField()[celli] = dopingAmp*std::exp(double(-pow(.5*(cc.component(0)-xDope)/dopingRadius,2)))*
									   std::exp(double(-pow(.5*(cc.component(2)-zDope)/dopingLength,2)))*
									   std::exp(double(-pow(.5*(cc.component(1)-yDope)/dopingRadius,2))); 
								   }
			//Info << "doping at "		<<cc		<<"is "		<<SIGMA_DOPING.internalField()[celli]		<<"\n";
		}
		else
			{SIGMA_DOPING.internalField()[celli] = 0;}
	}
    
    Info << "properties fields initialized"<<"\n";

