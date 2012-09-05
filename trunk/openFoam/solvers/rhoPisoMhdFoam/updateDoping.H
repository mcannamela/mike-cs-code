{
scalarField& sigmaCells = SIGMA.internalField();
scalarField& sdCells = SIGMA_DOPING.internalField();

dopingOn = readBool(MHDControls.subDict("sigmaDoping").lookup("dopingOn"));
if (dopingOn)
{

	//Info << "	max doping is " << max(SIGMA_DOPING) << "\n";
	forAll(sdCells, celli)
		{
			if (sdCells[celli]>sigmaCells[celli])
			sigmaCells[celli] = sdCells[celli];
		}
}
forAll(mesh.boundary(), patchi){
	//Info << "on patch "<<patchi<<endl;
	bool fixesInternal = true;
	forAll(patches, p)
	{
		//Info << "patch is " << patches[p] << " ID is " << mesh.boundary().findPatchID(patches[p]) << endl;
		if (int(mesh.boundary().findPatchID(patches[p]))==patchi)
		{
			// Info << "setting fixes internal to " << readBool(MHDControls.subDict("patches").lookup(patches[p])) <<endl;
			 fixesInternal = readBool(MHDControls.subDict("patches").lookup(patches[p]));
		 }
			 
	}
	
	//Info << patches[patchi] << " fixes internal sigma?  "<< fixesInternal<<endl;
	if (SIGMA.boundaryField()[patchi].fixesValue() && fixesInternal)
	{
		forAll(mesh.boundary()[patchi], facei)
		{
		sigmaCells[mesh.boundary()[patchi].faceCells()[facei]] = SIGMA.boundaryField()[patchi][facei];
		
		}
	}
}
//Info << "	doping updated"<< endl;
}
