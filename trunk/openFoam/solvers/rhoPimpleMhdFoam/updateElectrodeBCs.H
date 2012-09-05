{
Info << "\n::::::::::::::::::Begin electrode BC update::::::::::::::::::::\n";

vectorField& JCells = J.internalField();
//vectorField& adC = anodeDirection.internalField();
//volScalarField magJ(mag(J));
scalarField& qA = QAnode.internalField();
scalarField& qC = QCathode.internalField();
scalarField& s = SIGMA.internalField();
scalarField& sBL = SIGMA_BL.internalField();
scalarField& beta = betaSigma.internalField();
scalarField& jHCells = jouleHeating.internalField();




//scalarField& V = mesh.V().internalField();

//Info<< "	Computing anode BL conductivity" << endl;
forAll(sBL, celli)
{
	cc = mesh.C().internalField()[celli];
	if (electronTemperatureBoundaryConductivity)
	{
		if (cc.x()<.038)
			sBL[celli] = SIGMA_Te.internalField()[celli];	
		else
			sBL[celli] = s[celli];
			
		//if (cc.z()<0)
			//sBL[celli] *=1e-7;
		
    }		
	else
		{
		//sBL[celli] = anodeSigmaTe+plasma.sigma(H.internalField()[celli]);	
		sBL[celli] = plasma.sigma(H.internalField()[celli]);	
	}
}

//Info << "fixedAnodeConductivity " << fixedAnodeConductivity <<endl;
//Info << "!fixedAnodeConductivity " << !fixedAnodeConductivity <<endl;
if (electronTemperatureBoundaryConductivity)
{
	//Info<< "blending conductivity from BL to bulk" << endl;
	forAll(s,celli)
	{
		s[celli] = beta[celli]*s[celli]+(1-beta[celli])*sBL[celli];
	}
	
}

if (fixedAnodeConductivity)
{
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
			s[mesh.boundary()[patchi].faceCells()[facei]] = SIGMA.boundaryField()[patchi][facei];
			
			}
		}
	}
}

if (dopingOn)
{
	forAll(s,celli)
	{
		if(SIGMA_DOPING.internalField()[celli]>s[celli])
		{
			s[celli] = SIGMA_DOPING.internalField()[celli];
		}
		
	}
}



	
double blend = 0;
double arg, num, denom;
if (axialConductivityDecayOn)
{
	forAll(s, celli)
		{
			cc = mesh.C().internalField()[celli];
				
			arg = 2.0*double(cc.component(0)-.042)/.01;
			arg = std::max(std::min(arg, 20.0), -20.0);
			//Info << "arg "<<"is "		<<arg		<<"\n";
			num = (std::exp(arg)-1);
			//Info << "num is "		<<"is "		<<num		<<"\n";
			denom = (1e-10+std::exp(arg)+1);
			//Info << "denom is "		<<"is "				<<denom		<<"\n";
			
			blend = 1-(.5+.5*num/denom);
			s[celli]*= std::min(beta[celli]+blend,1.0)+1;

		}
}

double totalAnodeCurrent = 0;

double totalCathodeCurrent = 0;
vector cc;
double r;
forAll(QAnode.boundaryField(), patchi)
{
	//Info<< "looping over patches..." << patchi << endl;
	
	//fvPatchVectorField& pJ = J.boundaryField()[patchi];
	fvPatchScalarField& pQA = QAnode.boundaryField()[patchi];
	fvPatchScalarField& pQC = QCathode.boundaryField()[patchi];
	fvPatchScalarField& pS = SIGMA.boundaryField()[patchi];
	
	
	//const vectorField&  A = mesh.boundary()[patchi].Sf();
	label faceCelli =-1;

	
	if (mesh.boundary()[patchi].name() == "anode")
	{
		//Info<< "Computing anode heat flux" << endl;
		
		forAll(pQA, facei)
		{
			totalAnodeCurrent+= mag(phiJ.boundaryField()[patchi][facei]);	
		}
		
	}

    
	if (mesh.boundary()[patchi].name() == "cathode")
	{
		//Info<< "Computing cathode conductivity and heat flux" << endl;
		forAll(pQC, facei)
		{
			faceCelli = mesh.boundary()[patchi].faceCells()[facei];
			cc = mesh.C().internalField()[faceCelli];
			r = mag(cc-cathodePoint);
			//r = std::pow(std::pow(cc.component(1),2)+pow(cc.component(2),2),.5);
			//magic number is 2*k_b/e; anode work function must be in eV
			totalCathodeCurrent+= mag(phiJ.boundaryField()[patchi][facei]);
			if (cathodeRiggingOn)
			{
				if (r<=spotRadius)
					{
						//Info << "cathode point hit!\n";
						//qC[faceCelli] = -cathodeHeatFlux*(1-std::pow(r/spotRadius, cathShapePower))*mag(A[facei])/mesh.V()[faceCelli];
						//totalCathodeHeat+= mag(qC[faceCelli]*mesh.V()[faceCelli]);
						
						s[faceCelli] = sigmaCath;
						pS[facei] = sigmaCath;
						//kappa.internalField()[faceCelli] = 2;
						//kappa.boundaryField()[patchi][facei] = 2;
						Temperature.internalField()[faceCelli] = cathodeWallTemperature+1000;
						H.internalField()[faceCelli] = plasma.h(Temperature.internalField()[faceCelli]);
		
						//Te.internalField()[faceCelli] = cathodeWallTemperature*3;
										
						MU.internalField()[faceCelli] = plasma.mu(H.internalField()[faceCelli]);
						kappa.internalField()[faceCelli] = plasma.k(H.internalField()[faceCelli]);
						C.internalField()[faceCelli] = plasma.c(H.internalField()[faceCelli]);
						
					}
				else
					{
						s[faceCelli] = 1e-1;
						pS[facei] = 1e-1;
					}
			}
		}
		
	
			
	}
}
reduce(totalAnodeCurrent, sumOp<double>());
reduce(totalCathodeCurrent, sumOp<double>());

//Info << "computing anode heat loss" << endl;
double totalJouleHeating = 0;
double totalConductionLoss = 0;
double totalRadiativeLoss = 0;
double totalAnodeHeat = 0;
double totalCathodeHeat = 0;
QAnode*=0;

double cD = 0;

forAll(qA, celli)
{
	qA[celli] = -mag(JCells[celli])*(1-beta[celli]);
	totalAnodeHeat += mag(qA[celli]*mesh.V()[celli]);
	totalJouleHeating+=jHCells[celli]*mesh.V()[celli];
	totalRadiativeLoss+=Qrad[celli]*mesh.V()[celli];
	totalConductionLoss+= (TDiffusion.internalField()[celli]+turbTDiffusion.internalField()[celli])*mesh.V()[celli];
	
	cc = mesh.C().internalField()[celli];
	//rSq = sqr(mag(cc.y()))+sqr(mag(cc.z()));
	r = std::pow(double(mag(sqr(mag(cc.y()))+sqr(mag(cc.z())))),.5);
	cD = cathodeDistance.internalField()[celli];
	qC[celli] = -std::exp(-mag(std::pow(r/(spotRadius),6.0)))*std::exp(-mag(cD/3e-4));
	if (cathodeRiggingOn)
	{
		if (r<(1.3*spotRadius) && cD<1e-3)
			{s[celli] = std::max(s[celli], sigmaCath*std::exp(-mag(std::pow(r/(1.1*spotRadius),6.0)))*std::exp(-mag(cD/8e-4)));}
	}
	totalCathodeHeat+=mag(qC[celli]*mesh.V()[celli]);
}
reduce(totalJouleHeating, sumOp<double>());
reduce(totalRadiativeLoss, sumOp<double>());
reduce(totalConductionLoss, sumOp<double>());
reduce(totalAnodeHeat, sumOp<double>());
reduce(totalCathodeHeat, sumOp<double>());
//Info <<"renormalizing anode heat loss" << endl;
forAll(qA, celli)
{
	
	qA[celli] *= min(1.0,(totalAnodeCurrent/totalCurrent))*anodeHeat/(1e-6+totalAnodeHeat);
	if (Temperature.internalField()[celli]<anodeWallTemperature)
		qA[celli] = max(qA[celli], -jouleHeating.internalField()[celli]);
	
	qC[celli] *= min((totalCathodeCurrent/totalCurrent),1.0)*heatToCathode/totalCathodeHeat;
	
	if (Temperature.internalField()[celli]<cathodeWallTemperature)
		qC[celli] = max(qC[celli], -jouleHeating.internalField()[celli]);
}


totalAnodeHeat=0;
totalCathodeHeat=0;
forAll(qA, celli)
{
	totalAnodeHeat += qA[celli]*mesh.V()[celli];
	totalCathodeHeat+=qC[celli]*mesh.V()[celli];
}
//totalConductionLoss = mag(fvc::domainIntegrate(TDiffusion)+fvc::domainIntegrate(turbTDiffusion)).value();
reduce(totalAnodeHeat, sumOp<double>());
reduce(totalCathodeHeat, sumOp<double>());
Info << "	heat sunk to anode is "<< totalAnodeHeat << " and should be "<< -anodeHeat*min(1.0,(totalAnodeCurrent/totalCurrent))<<endl;
Info << "	heat sunk to cathode is "<< totalCathodeHeat << " and should be "<< -heatToCathode*min((totalCathodeCurrent/totalCurrent),1.0)<<endl;

Info << "	patch currents:"<< "  ";
double iTot = 0;
forAll(phiJ.boundaryField(), patchi)
{
	iTot=0;
	forAll(phiJ.boundaryField()[patchi], facei)
	{
		iTot+=phiJ.boundaryField()[patchi][facei];
	}
	reduce(iTot, sumOp<double>());
	Info << mesh.boundary()[patchi].name() << " "<< iTot << ", ";
}
Info <<endl;
Info << "	total joule heating is "<< totalJouleHeating << endl;
Info << "	total radiative loss is "<< totalRadiativeLoss << endl;
Info << "	total conductive loss is "<< totalConductionLoss << endl;
Info << "	net power is approximately "<< totalJouleHeating -totalRadiativeLoss +totalAnodeHeat +totalCathodeHeat+totalConductionLoss<< endl;
Info << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n\n";
}

