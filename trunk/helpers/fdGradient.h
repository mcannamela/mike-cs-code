//#include "gasMixture.h"
#include <vector>
#include <iostream>
typedef vector<double> 	dVec;
typedef vector<dVec> 	dMat;

//need to define your own parameterContainer struct to work with your function f, like 

/*struct parameterContainer
{
int 	parOne
double 	parTwo
};*/
struct parameterContainer
{
	
		dVec 		boundaryConstant;
		dVec 		Xcore;
		double		rTotal;
		double 		rCore;
		double 		energyIn;
		double 		momentumIn;
		gasMixture  g;
};

dMat simpleInv(dMat A) //dumb way to find the inverse of a 2x2 matrix
{
	double detA;
	dVec zVec(2,0);
	dMat invA;
	invA.resize(2,zVec);
	
	detA= A[0][0]*A[1][1]- A[0][1]*A[1][0];
	
	invA[0][0]= A[1][1]/detA;
	invA[0][1]= -A[0][1]/detA;
	invA[1][0]= -A[1][0]/detA;
	invA[1][1]= A[0][0]/detA;

	return invA;
};

//finite difference gradient of a vector valued function
dVec fdGradient(dVec x, dVec dx, dVec (*f)(dVec, parameterContainer&), parameterContainer& p)
{
	dVec 	F = f(x,p);
	int		nOutputs = F.size();
	if (nOutputs != 1)
		cout << "error, fdGradient only works on scalar functions" << endl;
	int		nInputs = x.size();
	dVec	h(0,nInputs);
	
	dMat 	Xplus(nInputs,x);//ith row is x plus a step of dx(i) in the ith direction  
	dMat 	Fplus(nInputs,F);//ith row is f evaluated at the ith element in Xplus
	dMat 	X2plus(nInputs,x);
	dMat 	F2plus(nInputs,F);
	
	dVec	G(nInputs,0);
	
	for(int i = 0; i<nInputs; i++)//compute values of f at x+dx
	{
		Xplus[i][i]	+= dx[i];
		X2plus[i][i]+= 2*dx[i];
		Fplus[i] 	=  f(Xplus[i],p);
		F2plus[i]	=  f(X2plus[i],p);
		
	};
	
	for (int i=0; i<nOutputs; i++)
	{	
		for(int j=0; j<nInputs; j++)//2nd order finite difference approximation
			{
				G[j] = -(F2plus[j][i] - 4*Fplus[j][i] +3*F[i])/(2*dx[j]);
			};
	};
	
	return G;
};

dMat fdHessian(dVec x, dVec dx, dVec (*f)(dVec, parameterContainer&), parameterContainer& p, dVec* G)
{
	dVec 	F = fdGradient(x,dx,f,p);
	*G        = F;
	int		nOutputs = F.size();
	int		nInputs = x.size();
	dVec	h(0,nInputs);
	
	dMat 	Xplus(nInputs,x);//ith row is x plus a step of dx(i) in the ith direction  
	dMat 	Fplus(nInputs,F);//ith row is f evaluated at the ith element in Xplus
	dMat 	X2plus(nInputs,x);
	dMat 	F2plus(nInputs,F);
	
	dMat	H(nOutputs,x);
	
	for(int i = 0; i<nInputs; i++)//compute values of f at x+dx
	{
		Xplus[i][i]	+= dx[i];
		X2plus[i][i]+= 2*dx[i];
		Fplus[i] 	=  fdGradient(Xplus[i],dx,f,p);
		F2plus[i]	=  fdGradient(X2plus[i],dx,f,p);
		
	};
	
	for (int i=0; i<nOutputs; i++)
	{	
		for(int j=0; j<nInputs; j++)//2nd order finite difference approximation
			{
				H[i][j] = -(F2plus[j][i] - 4*Fplus[j][i] +3*F[i])/(2*dx[j]);
			};
	};
	return H;
};