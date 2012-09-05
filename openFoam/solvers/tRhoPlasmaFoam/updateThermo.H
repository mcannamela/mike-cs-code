{
scalarField& TCells = Temperature.internalField();
scalarField& psiCells = psi.internalField();
scalarField& pCells = p.internalField();
scalarField& rhoCells = rho.internalField();
scalarField& fCells = f.internalField();

//Info << "updating psi, rho in cells\n";
forAll(TCells, celli)
{
	psiCells[celli] = (1/(fCells[celli]/plasmaMolarMass+(1-fCells[celli])/ambientMolarMass))/(8314*TCells[celli]);
	rhoCells[celli] = psiCells[celli]*pCells[celli];
}



forAll(T.boundaryField(), patchi)
{
	
	fvPatchScalarField& pT = Temperature.boundaryField()[patchi];
	fvPatchScalarField& ppsi = psi.boundaryField()[patchi];
	fvPatchScalarField& prho = rho.boundaryField()[patchi];
	fvPatchScalarField& pp = p.boundaryField()[patchi];
	fvPatchScalarField& pf = f.boundaryField()[patchi];
//	Info << "updating psi, rho in boundary patch "<<patchi<<"\n";
	//if (pT.fixesValue())
	//{
		forAll(pT, facei)
			{
				ppsi[facei] = (1/(pf[facei]/plasmaMolarMass+(1-pf[facei])/ambientMolarMass))/(8314*pT[facei]);
				prho[facei] = ppsi[facei]*pp[facei];
			}
	//}
}

rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
}

