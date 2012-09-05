{
	MU_eddy.internalField() = turbulence->muEff()-turbulence->mu();//(pow(smagorinskyConstant*pow(mesh.V(), 1.0/3.0),2)*sqrt( (symm(fvc::grad(U)) && symm(fvc::grad(U)) ) ));
	K_eddy = C_p*(turbulence->alphaEff()-turbulence->alpha());//MU_eddy/.7;
	
}
