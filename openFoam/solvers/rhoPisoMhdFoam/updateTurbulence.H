{
	MU_eddy.internalField() = turbulence->muEff()-turbulence->mu();//(pow(smagorinskyConstant*pow(mesh.V(), 1.0/3.0),2)*sqrt( (symm(fvc::grad(U)) && symm(fvc::grad(U)) ) ));
	kappa_eddy = C*(turbulence->alphaEff()-turbulence->alpha())+dimensionedScalar("minKappa", dimensionSet(1, 1, -3, -1, 0, 0, 0), .01);//MU_eddy/.7;
	
}
