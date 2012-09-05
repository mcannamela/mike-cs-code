#ifndef MESHCELL_H
#define MESHCELL_H
#include <cmath>
#include <vector>
#include <iostream>
//here is a change
using namespace std;
typedef vector<double>	dVec;
typedef vector<dVec>	dMat;

class meshCell//
{
private:
	double 	r; //radial coord and dimension, (r, dr)
	double	z; //axial coord, (z,dz)
	double	dr;
	double 	dz;
	double	dTheta;
	dVec 	X; //dependent variables (velocity, energy, quantity of interest)
	double	V; //element volume
	dVec	A; //element areas
	
	int outerIdx;
	int innerIdx;
	int faceIdx;
	
	void computeArea();
	void computeVolume();

public:
	//constructors
	void meshCell();
		
	//setters
	void set(char, double);//set one of r,z,dr,dz
	void set(char,dVec);//set (r,dr) or (z,dz)
	void set(int, double);//set x[idx]
	void set(dVec);// set x
	void set(dVec,dVec);//set (r,dr),(z,dz)
	
	//getters
	double Volume()
		{return V;};
	double AreaOuter()
		{return A[outerIdx];};
	double AreaInner()
		{return A[innerIdx];};
	double AreaFace()
		{return A[faceIdx];};
	double R()
		{return r[1]};
	double dR()
		{return r[2]};
	double Z()
		{return z[1]};
	double dZ()
		{return z[2]};
};

#endif
