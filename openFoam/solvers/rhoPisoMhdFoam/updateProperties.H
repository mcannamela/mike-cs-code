Info<< "\n.................Begin properties update...................\n";
{
scalarField& HCells = H.internalField();
	scalarField& TCells = Temperature.internalField();
	//vectorField& JCells = J.internalField();
	scalarField& QradCells = Qrad.internalField();
	scalarField& sigmaCells = SIGMA.internalField();
	scalarField& sigmaTeCells = SIGMA_Te.internalField();
	//scalarField& sdCells = SIGMA_DOPING.internalField();
	scalarField& muCells = MU.internalField();
	scalarField& KCells = kappa.internalField();
	scalarField& CCells = C.internalField();
	scalarField& pheeCells = phee.internalField();


double pheeMax = -1e37;
double pheeMin = 1e37;
//double kMax = 1e-37;
double Tmax = 0;
double Tmin = 1e37;

//vector cc;



forAll(pheeCells, celli){
	if (pheeCells[celli]>pheeMax)
		pheeMax = pheeCells[celli];
	if (pheeCells[celli]<pheeMin)
		pheeMin = pheeCells[celli];
	}
reduce(pheeMax, maxOp<double>());
reduce(pheeMin, minOp<double>());

Info <<"	phee min, max: "<< pheeMin
	<<"     "<< pheeMax<<"\n";


Hmin = 1e15;

//Info<< "Computing cell values\n" << endl;
forAll(TCells, celli)
{
	
	TCells[celli] = plasma.T(HCells[celli]);
	
	if (TCells[celli]<150)
		Info << "				low temperature detected: T = " << TCells[celli] << endl;
		
	TCells[celli] = max(TCells[celli], 250);
	HCells[celli] = plasma.h(TCells[celli]);
	
	//Te.internalField()[celli] = betaSigma.internalField()[celli]*TCells[celli]+(1-betaSigma.internalField()[celli])*anodeElectronTemperature;
    sigmaCells[celli] = plasma.sigma(HCells[celli]);
    sigmaTeCells[celli] = plasma.sigma(plasma.h(min(Te.internalField()[celli],lastTMax )));
    //sigmaCells[celli] = plasma.sigma(plasma.h(Te.internalField()[celli]));
    
	muCells[celli] = plasma.mu(HCells[celli]);
	KCells[celli] = plasma.k(HCells[celli]);
	CCells[celli] = plasma.c(HCells[celli]);
	
	if( TCells[celli]>Tmax)
		{Tmax = TCells[celli];}
	if( TCells[celli]<Tmin)
		{Tmin = TCells[celli];}
	if( HCells[celli]<Hmin)
		{Hmin = HCells[celli];}
		
	
	
	if (radiationOn)
		QradCells[celli] = radiatiativeLoss(TCells[celli]);
	else
		QradCells[celli] = 0;
	
	
}
reduce(Tmax, maxOp<double>());
reduce(Tmin, minOp<double>());
reduce(Hmin, minOp<double>());

lastTMax = Tmax;


forAll(Temperature.boundaryField(), patchi)
{
	
	fvPatchScalarField& pT = Temperature.boundaryField()[patchi];
	
	fvPatchScalarField& pQrad = Qrad.boundaryField()[patchi];
	fvPatchScalarField& pH = H.boundaryField()[patchi];
	fvPatchScalarField& psigma = SIGMA.boundaryField()[patchi];
	fvPatchScalarField& psigmaBL = SIGMA_BL.boundaryField()[patchi];
	
	fvPatchScalarField& pmu = MU.boundaryField()[patchi];
	fvPatchScalarField& pK = kappa.boundaryField()[patchi];
	fvPatchScalarField& pC = C.boundaryField()[patchi];

	
	/*if (!psigma.fixesValue())
	{
		Info<< "		patch "<< patchi << " does not fix sigma\n" << endl;
		
	}*/
	
	if (pT.fixesValue())
	{
		forAll(pT, facei)
		{
			
			if (!psigma.fixesValue())
			    {
					psigma[facei] = plasma.sigma(pH[facei]);
					//psigma[facei] = plasma.sigma(plasma.h(Te.boundaryField()[patchi][facei]));
			    }
			psigmaBL[facei] = psigma[facei];    
			pmu[facei] = plasma.mu(pH[facei]);
			pK[facei] = plasma.k(pH[facei]);
			
			
			pC[facei] = plasma.c(pH[facei]);
			pQrad[facei] = radiatiativeLoss(pT[facei]);
		}
	}
	else
	{
		forAll(pT, facei)
		{
			pT[facei] = plasma.T(pH[facei]);
			if (!psigma.fixesValue())
			    {
					psigma[facei] = plasma.sigma(pH[facei]);
					//psigma[facei] = plasma.sigma(plasma.h(Te.boundaryField()[patchi][facei]));
			    }
			psigmaBL[facei] = psigma[facei];    
			pmu[facei] = plasma.mu(pH[facei]);
			pK[facei] = plasma.k(pH[facei]);
			pC[facei] = plasma.c(pH[facei]);
			pQrad[facei] = radiatiativeLoss(pT[facei]);
		}
	}
}

//Info<< "max k = " << kMax<<"\n";
//Info<< "				properties have been updated!\n";
Temperature.correctBoundaryConditions();
//Info<< "				temperature BC corrected.\n";

forAll(Temperature.boundaryField(), patchi)
{
	fvPatchScalarField& pT = Temperature.boundaryField()[patchi];
	fvPatchScalarField& pH = H.boundaryField()[patchi];

		forAll(pT, facei)
		{
			pH[facei] = plasma.h(pT[facei]);
		}	
}

Info << "	T min, max: "<< Tmin << "  " << Tmax<< endl;
Info<< ".............................................................\n";
}

