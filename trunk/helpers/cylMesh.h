#ifndef CYLMESH_H
#define CYLMESH_H

#include "fileHelpers.h"
#include "cylMeshCell.h"
#include <vector>
#include <fstream>
#include <iostream>
//here is a change
using namespace std;
typedef vector<double>	dVec;
typedef vector<dVec>	dMat;

class cylMesh//axisymmetric 
{

private:
	//holds the cells of the mesh
	vector<vector<cylMeshCell>>		cells;	// (z,r)
	cylMeshCell 					nullCell; //empty cell for expanding the mesh
	vector<cylMeshCell>				nullCellVec;
	
	double dZdefault;
	double dRdefault;
	//current cell index
	int rIdx;
	int zIdx;
	
public:
	//constructors
	cylMesh();
	cylMesh(double, double);
	
	//mesh expansion
	void addCell(char)
	void addCell(char, double)
		
	//mesh navigation
	void Next(void);	//next axial el
		{zIdx++;};
	void Prev(void);	//previous axial el
		{zIdx--;};
	void Out(void);		//radial next el outward
		{rIdx++;};
	void In(void);		//radial previous element inward
		{rIdx--;};
	
	//setters
	void setIdx(int z,int r)
		{zIdx = z;
		rIdx = r;};
	void set(int idx, double val);//set x[idx]
		{cells[zIdx][rIdx].set(idx,val)};
	void set(dVec vals);// set x
		{cells[zIdx][rIdx].set(vals)};
	
	//getters
	dVec Idx(void)
		{dVec idx(2,0);
		idx[0] = zIdx;
		idx[1] = rIdx
		return idx};
		return idx};
	
	double R(void)
		{return cells[zIdx][rIdx].R()};
	double dR(void)
		{return cells[zIdx][rIdx].dR()};
	double Z(void)
		{return cells[zIdx][rIdx].Z()};
	double dZ(void)
		{return cells[zIdx][rIdx].dZ()};
	
	double R(int,int)
		{return cells[zIdx][rIdx].R()};
	double dR(int,int)
		{return cells[zIdx][rIdx].dR()};
	double Z(int,int)
		{return cells[zIdx][rIdx].Z()};
	double dZ(int,int)
		{return cells[zIdx][rIdx].dZ()};
	
	double Volume()
		{return cells[zIdx][rIdx].V();};
	double AreaOuter()
		{return cells[zIdx][rIdx].AreaOuter();};
	double AreaInner()
		{return cells[zIdx][rIdx].AreaInner();};
	double AreaFace()
		{return cells[zIdx][rIdx].AreaFace();};
};

#endif
