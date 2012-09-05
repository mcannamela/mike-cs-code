// Solve the Momentum equation
/*vector spongeDir(1,0,0);
double blend = 0;
forAll(spongeFriction.internalField(), celli)
{
	cc = mesh.C().internalField()[celli];
		
	arg = 2.0*double(cc.component(0)-.01)/.006;
	arg = std::max(std::min(arg, 20.0), -20.0);
	//Info << "arg "<<"is "		<<arg		<<"\n";
	num = (std::exp(arg)-1);
	//Info << "num is "		<<"is "		<<num		<<"\n";
	denom = (1e-10+std::exp(arg)+1);
	//Info << "denom is "		<<"is "				<<denom		<<"\n";
	
	blend = 1-(.5+.5*num/denom);
		
	spongeFriction.internalField()[celli] = blend*1e6*(U.internalField()[celli]-spongeDir*3.5*anodeDistance.internalField()[celli]*cathodeDistance.internalField()[celli]/(.0015*.0015));
}*/
double blend = 0;
forAll(spongeFriction.internalField(), celli)
	{
		cc = mesh.C().internalField()[celli];
			
		arg = 2.0*double(cc.component(0)-spongeOnset)/spongeLength;
		arg = std::max(std::min(arg, 20.0), -20.0);
		//Info << "arg "<<"is "		<<arg		<<"\n";
		num = (std::exp(arg)-1);
		//Info << "num is "		<<"is "		<<num		<<"\n";
		denom = (1e-10+std::exp(arg)+1);
		//Info << "denom is "		<<"is "				<<denom		<<"\n";
		
		blend = 1-(.5+.5*num/denom);
		
		if (cc.component(0)< (spongeOnset+3*spongeLength))
				spongeFriction.internalField()[celli] = blend*spongeFrictionCoefficient*(U.internalField()[celli]-spongeVelocity*anodeDistance.internalField()[celli]*cathodeDistance.internalField()[celli]/(.0015*.0015));
		else
			spongeFriction.internalField()[celli] = vector(0,0,0);
		
	}

//vector forceCenter(.0265,0,0);
//vector forceDir(-.50,0,1);
//double forceMag = 5e6;
//double forceRadius = .00135;
double r = 0;

if (runTime.time().value()<endTime)
{
	Info << "tip force is on until "<< endTime << "\n";
	forAll(tipForce.internalField(), celli)
		{
			
			cc = mesh.C().internalField()[celli];
			r = mag(cc-forceCenter);
			
			//tipForce.internalField()[celli] = forceDir*forceMag*( std::exp(-mag(std::pow(r/(forceRadius),4.0)))+std::exp(-mag(std::pow(r2/(forceRadius),4.0))) );
			tipForce.internalField()[celli] = forceDir*forceMag*( std::exp(-mag(std::pow(r/(forceRadius),4.0))) )*mesh.V()[celli]/(2.33236152e-11);
		}
	}
else
	{tipForce*=0;}
	
if (runTime.time().value()<nookEndTime)
{
	Info << "nook force is on until "<< endTime << "\n";
	forAll(nookForce.internalField(), celli)
		{
			
			cc = mesh.C().internalField()[celli];
			r = mag(cc-nookForceCenter);
			
			//tipForce.internalField()[celli] = forceDir*forceMag*( std::exp(-mag(std::pow(r/(forceRadius),4.0)))+std::exp(-mag(std::pow(r2/(forceRadius),4.0))) );
			nookForce.internalField()[celli] = nookForceDir*nookForceMag*( std::exp(-mag(std::pow(r/(nookForceRadius),4.0))) )*mesh.V()[celli]/(2.33236152e-11);
		}
	}
else
	{nookForce*=0;}
	
vector oForceDir(-1,0,0);
Info << "mDotActual-mDotSet  = "<< 	(mDotActual-mDotSet) << "\n";
forAll(outletForce.internalField(), celli)
	{
		cc = mesh.C().internalField()[celli];
		//Info << "scaling factor " << mesh.V()[celli]*rho.internalField()[celli]/(1.6*2.33236152e-11) <<"\n";
		arg = 2.0*double(cc.component(0)-.041)/.003;
		arg = std::max(std::min(arg, 20.0), -20.0);
		//Info << "arg "<<"is "		<<arg		<<"\n";
		num = (std::exp(arg)-1);
		//Info << "num is "		<<"is "		<<num		<<"\n";
		denom = (1e-10+std::exp(arg)+1);
		//Info << "denom is "		<<"is "				<<denom		<<"\n";
		
		blend = (.5+.5*num/denom);
		if (mag(mDotActual-mDotSet)>.0001 && U.internalField()[celli].component(0)>5.0)
			outletForce[celli]= (mDotActual-mDotSet)*blend*pGain*oForceDir*rho.internalField()[celli]*mesh.V()[celli]/(1.6*2.33236152e-11);
		else
			outletForce[celli] = vector(0,0,0);
		
		
	}
	
tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  //+ turbulence->divDevRhoReff(U)
  - fvm::laplacian(MU+MU_eddy, U)
  - (lorenzForce)
  ==
  -spongeFriction
  +tipForce
  +outletForce
  +nookForce
);

//UEqn().relax();

volScalarField rAU(1.0/UEqn().A());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));
    K = 0.5*magSqr(U);
}
