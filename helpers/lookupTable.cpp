#include "lookupTable.h"

//make an empty table object
lookupTable::lookupTable()//put a zero in the first line of the table to start
{
	dVec aVector;
	aVector.push_back(0);
	t.resize(1,aVector);
};


//make a table from the first nLines lines in the file fName
lookupTable::lookupTable(fs::path p, int nLines)
{
	dMat		X;//to hold data from file
	ifstream 	dataFile;
	
	dataFile.open(p.file_string().c_str());
		
	if (dataFile.fail())
		throw "ERROR, file could not be opened \n";

	X = file2mat(dataFile, nLines);//read nLines lines from dataFile and put them into a vector of vectors X
	
	for(unsigned int i=0; i<X.size(); i++)
		setLine(i,&X[i]);
};
lookupTable::lookupTable(char* fName, int nLines)
{
	dMat		X;//to hold data from file
	ifstream 	dataFile;
		
	dataFile.open(fName);
	//cout << "open succeeded?  "   << !dataFile.fail() << endl;
	
	if (dataFile.fail())
		throw "ERROR, file could not be opened \n";

	X = file2mat(dataFile, nLines);//read nLines lines from dataFile and put them into a vector of vectors X
	
	for(unsigned int i=0; i<X.size(); i++)
		setLine(i,&X[i]);
};

//read the table from a file
void lookupTable::fromFile(char* fName,int nLines)
{
	lookupTable a(fName,nLines);
	t=a.t;
};

//write the table to a file
void lookupTable::toFile(char* fName)
{
	ofstream 	dataFile;
		
	dataFile.open(fName);
		
	if (dataFile.fail())
		throw "ERROR, file could not be opened \n";

	for (unsigned int i=0; i<t.size() ; i++)
	{
	for(unsigned int j=0; j<t[0].size(); j++)
		dataFile << t[i][j] << "   ";
	
	dataFile << "\n";
	}
}

//puts the vector of values tableValues into the line tableLine of the table
void lookupTable::setLine(int tableLine, dVec* tableValues)//(lineToSet,vectorToSet)
	{	int		tSize;
		dVec	aVector(1,0); //need a 0 vector to expand the table
		tSize = t.size();
	
		if ( (tSize-1) < tableLine ) //if there are not enough lines in the table...
		{
			for (int i=1; i <= (tableLine-tSize+1) ;  i++)//...then add more lines
			{
				tSize =  t.size();
				t.push_back(aVector);
			};
		};
		t[tableLine].resize((*tableValues).size(),0);//expand the line of interest out to the size of the input vector
		t[tableLine]=*tableValues;
	};
	
dVec lookupTable::getLine(int tableLine)
	{
		return t[tableLine];
	};
	
//prints a line of the table so you can see its contents
void lookupTable::printLine(int idx)
{
	cout << "line " << idx << " is: ";
	for (unsigned int i = 0; i<t[idx].size(); i++)
	cout << t[idx][i] << " ";
	cout << endl;
};

//returns a value from the table after interpolating linearly
double lookupTable::getVal(int fromLine, int toLine, double fromVal)//(inputTableLine, outputTableLine, inputValue), find the value in toLine corresponding to fromVal in fromLine
	{
		int idx; //table entry index
		double val; //value corresponding to fromVal
		
		idx = binaryBracketSearch(fromLine, fromVal);//find the relevant table entries
		
		if (idx == t[fromLine].size()-1) //if we have the very last element, extrapolate using slope of last two elements
			//return t[toLine][t[fromLine].size()-1]; //use last value
			idx=idx-1;//extrapolate
		
		val = t[toLine][idx]+(t[toLine][idx+1]-t[toLine][idx])*( fromVal-t[fromLine][idx] )/( t[fromLine][idx+1]-t[fromLine][idx] );//linear interpolation
		return val;
	};

//finds values in the table 
int lookupTable::binaryBracketSearch(int tableLine,double val)//(lineOfTable, valueToBracket)
		{
			int  last   = t[tableLine].size()-1; //upper bound on search
			int  first  = 0; //lower bound on search
			int  middle = (last+first)/2; //element in the middle of the search space
			bool found  = false; //whether value has been found yet
			int  idx   = -1;
			
			//error checking for values out of table range
			if (t[tableLine][t[tableLine].size()-1]<=val)
			{
				//cout << "warning, value requested from table greater than largest table value" << endl;
				return (t[tableLine].size()-1);//return the last valid index
			}
			else if (t[tableLine][0]>=val)
			{
				//cout << "warning, value requested from table less than smallest table value" << endl;
				return 0;//return the first index
			};
			
			while (!found && (first<last))
			{
			    if ( (last+first)%2==0 )
				middle = (last+first)/2;
				else 
				middle = (last+first+1)/2;
								
				if (  (val >= t[tableLine][middle] ) &&  (val< t[tableLine][middle+1])  )//if the value is between middle and middle+1
				{
					idx   =   middle; 
					found =   true;							//search is over
				}				
				else if (val>t[tableLine][middle])															//if the value is in the upper half
				{
					first = middle;
				}						//forget about the first half
				else if ((last-first)==1)
				{
					idx = first;
					found = true;
				}
				else																												//value is in the lower half
				{
					last = middle;	//so forget about the upper half
				}
			};	//end while

			return idx;
		};

