#ifndef CYLMESHCELL_H
#define CYLMESHCELL_H
#include <cmath>
#include <vector>
#include <iostream>
//here is a change
using namespace std;
typedef vector<double>	dVec;
typedef vector<dVec>	dMat;

class CylMeshCell//
{
private:
	double 	r; //radial coord and dimension, (r, dr)
	double	z; //axial coord, (z,dz)
	double 	t; //circumferential coordinate
	double	dr;
	double 	dz;
	double	dt;
	double	pi;
	dVec 	X; //dependent variables (velocity, energy, quantity of interest)
	double	V; //element volume
	dVec	A; //element areas
	
	int outerIdx;
	int innerIdx;
	int faceIdx;
	int endIdx;
	
	void computeArea();
	void computeVolume();

public:
	//constructors
	cylMeshCell();
		
	//setters
	void set(char, double);//set one of r,z,dr,dz
	void setD(char, double);//set one of r,z,dr,dz
	void set(char,dVec);//set (r,dr) or (z,dz) or (t,dt)
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
		{return r};
	double dR()
		{return dr};
	double Z()
		{return z};
	double dZ()
		{return dz};
	double T()
		{return t};
	double dT()
		{return dt};	
};

#endif
