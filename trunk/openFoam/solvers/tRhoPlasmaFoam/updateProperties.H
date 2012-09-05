//Info<< "				begin properties update\n";

scalarField& HCells = H.internalField();
scalarField& TCells = Temperature.internalField();
scalarField& QradCells = Qrad.internalField();
scalarField& muCells = MU.internalField();
scalarField& KCells = K.internalField();
scalarField& CCells = C_p.internalField();
scalarField& fCells = f.internalField();

Tmax = 0;
Hmin = 1e15;

//Info<< "Computing cell values\n" << endl;
forAll(TCells, celli)
{
	//Ttmp = TCells[celli];
	//Info << "old T was "<<Ttmp;
	//Info << " gasBlend of " << fCells[celli] <<", "<<HCells[celli] << " gives "<< gasBlendTemperature(fCells[celli], HCells[celli]) << "\n";
	
	//Info <<"new T is "<< TCells[celli]<<"\n";
	HCells[celli] = std::max( HCells[celli], 250000.0);
	if (TCells[celli]<150.0)
		Info << "				low temperature detected: T = " << TCells[celli] 
		<< ",  plasmaFraction is " << fCells[celli]
		<<"\n";
	TCells[celli] = std::max(gasBlendTemperature(fCells[celli], HCells[celli]), 250.0);

	muCells[celli] = plasma.mu(HCells[celli])*fCells[celli]+(1-fCells[celli])*ambientGas.mu(HCells[celli]);
	KCells[celli] = plasma.k(HCells[celli])*fCells[celli]+(1-fCells[celli])*ambientGas.k(HCells[celli]);
	CCells[celli] = plasma.c(HCells[celli])*fCells[celli]+(1-fCells[celli])*ambientGas.c(HCells[celli]);
	
	if( TCells[celli]>Tmax)
		{Tmax = TCells[celli];}
	if( HCells[celli]<Hmin)
		{Hmin = HCells[celli];}
	
	if (radiationOn)
		QradCells[celli] = radiatiativeLoss(TCells[celli]);
	else
		QradCells[celli] = 0;
	
	
}



forAll(T.boundaryField(), patchi)
{
	
	fvPatchScalarField& pT = Temperature.boundaryField()[patchi];
	fvPatchScalarField& pQrad = Qrad.boundaryField()[patchi];
	fvPatchScalarField& pH = H.boundaryField()[patchi];
		
	fvPatchScalarField& pmu = MU.boundaryField()[patchi];
	fvPatchScalarField& pK = K.boundaryField()[patchi];
	fvPatchScalarField& pC = C_p.boundaryField()[patchi];
	fvPatchScalarField& pf = f.boundaryField()[patchi];
	
	
	
	if (pT.fixesValue())
	{
		forAll(pT, facei)
		{
			pmu[facei] = plasma.mu(pH[facei])*pf[facei]+(1-pf[facei])*ambientGas.mu(pH[facei]);
			pK[facei] = plasma.k(pH[facei])*pf[facei]+(1-pf[facei])*ambientGas.k(pH[facei]);
			pC[facei] = plasma.c(pH[facei])*pf[facei]+(1-pf[facei])*ambientGas.c(pH[facei]);
			pQrad[facei] = radiatiativeLoss(pT[facei]);
		}
	}
	else
	{
		forAll(pT, facei)
		{
			pT[facei] = gasBlendTemperature(pf[facei], pH[facei]);
			pmu[facei] = plasma.mu(pH[facei])*pf[facei]+(1-pf[facei])*ambientGas.mu(pH[facei]);
			pK[facei] = plasma.k(pH[facei])*pf[facei]+(1-pf[facei])*ambientGas.k(pH[facei]);
			pC[facei] = plasma.c(pH[facei])*pf[facei]+(1-pf[facei])*ambientGas.c(pH[facei]);
			pQrad[facei] = radiatiativeLoss(pT[facei]);
		}
	}
}

Temperature.correctBoundaryConditions();

forAll(Temperature.boundaryField(), patchi)
{
	fvPatchScalarField& pT = Temperature.boundaryField()[patchi];
	fvPatchScalarField& pH = H.boundaryField()[patchi];
	fvPatchScalarField& pf = f.boundaryField()[patchi];
	
//	if (pT.fixesValue())
	//{
		forAll(pT, facei)
		{
			pH[facei] = plasma.h(pT[facei])*pf[facei]+(1-pf[facei])*ambientGas.h(pT[facei]);
		}
	//}
}

