#include "cylMeshCell.h"

cylMeshCell::cylMeshCell();
{
	ind endIdx	 = 0
	int faceIdx  = 1;
	int innerIdx = 2 ;
	int outerIdx = 3 ;
	r = 0;
	z = 0;
	t = 0;
	dr = 0;
	dz = 0;
	pi = 3.141592;
	dt = 2*pi;
};

void cylMeshCell::computeVolume()
{
	V = (r*dTheta)*dr*dz;
};

void cylMeshCell::computeArea()
{
	A[faceIdx] 	= dt*( pow(r+dr,2)-pow(r,2) )/2;
	A[innerIdx] = dz*dt*r;
	A[outerIdx] = dz*dt*(r+dr);
	A[endIdx]	= dz*dr;
};

void cylMeshCell::set(char c, double val);//set one of r,z,dr,dz
{
	switch (c)
	{
		case 'r': r = val;
					break;
		case 'z': z = val;
					break;
		case 't': t = val;
					break;
	}
	computeArea();
	computeVolume();
};
void cylMeshCell::setD(char c, double val);//set one of r,z,dr,dz
{
	switch (c)
	{
		case 'r': dr = val;
					break;
		case 'z': dz = val;
					break;
		case 't': dt = val;
					break;
	}
	computeArea();
	computeVolume();
};
void cylMeshCell::set(char c,dVec vals);//set (r,dr) or (z,dz) or (t,dt)
{
	switch (c)
	{
		case 'r': 	r = vals[0];
					dr = val[1];
					break;
		case 'z': 	z = vals[0];
					dz = val[1];
					break;
		case 't': 	t = vals[0];
					dt = val[1];
					break;
	}
	computeArea();
	computeVolume();
};
void cylMeshCell::set(int idx, double val);//set x[idx]
{
	x[idx] = val;
};
void cylMeshCell::set(dVec vals);// set x
{
	for (int i=0; i<x.size(); i++)
	x[i] = vals[i];
	computeArea();
	computeVolume();
}
}
void cylMeshCell::set(dVec rr,dVec zz);//set (r,dr),(z,dz)
{
	r=rr[0];
	dr=rr[1];
	z=zz[0];
	dz=zz[1];
	computeArea();
	computeVolume();
};
