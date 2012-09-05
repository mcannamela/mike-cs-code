#ifndef MESH_H
#define MESH_H

#include "fileHelpers.h"
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
//here is a change
using namespace std;
typedef vector<double>	dVec;
typedef vector<dVec>	dMat;

template <class T>
class mesh//axisymmetric 
{
private:
	vector<vector<T>>	cells;	
	
public:
	//constructors
	void axiSymMeshLL();
	
	//mesh expansion
	void addNodeNext(axiSymMeshLL*);
	void addNodeUp(axiSymMeshLL*);
	
	//mesh manipulations
	double sumUp(dVec);
	
	//mesh navigation
	void Next();	//next axial el
	void Prev();	//previous axial el
	void Out();		//radial next el outward
	void In();		//radial previous element inward
	void First();	//first axial element at this radius
	void Last();	//last axial element at this radius
	void Outer();	//outermost radial element at this axial location
	void Inner();	//innermost radial element at this location

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
