#include "lookupTable.h"
#include <iostream>
#include <vector>
using namespace std;
//typedef vector<double> dVec;

int main(void)
{
	
	int nLines = 6;

	char* goodName = "testTable.txt";
	char* badName = "notAFileName.txt";
	lookupTable bTable;

	try
	{
		bTable.fromFile(badName,nLines);
	}
	catch(char* exceptionStr)
	{
		cout << exceptionStr;
		bTable.fromFile(goodName,nLines);
	}
	lookupTable aTable = bTable;
	for (int i=0;i<aTable.nLines();i++)
	aTable.printLine(i);
	return 0;
};