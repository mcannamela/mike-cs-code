//implements one-D table lookup and interpolation with binary search (lines of table must be ordered), 
#ifndef LOOKUPTABLE_H
#define LOOKUPTABLE_H

#include "fileHelpers.h"
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;
using namespace std;
typedef vector<double>	dVec;
typedef vector<dVec>	dMat;

class lookupTable
{
private:
	dMat 	t;  //holds the table values
	int 	binaryBracketSearch(int,double);//(lineOfTable, valueToBracket)
	
public:
	lookupTable();
	lookupTable(char*,int);
	lookupTable(fs::path p, int nLines);
	
	void			fromFile(char*,int);//read the table from  a file
	void			toFile(char*);//write the table to a file
	void			setLine(int, dVec*);//(lineToSet,vectorToSet)
	dVec			getLine(int);// return a line of the table
	int				nLines(void)
					{return t.size();};
	void			printLine(int);//use to display the table contents in the console
	double			getVal(int,int,double);//use to retrieve a value from the table using linear interpolation
};

#endif