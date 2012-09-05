#include "fileHelpers.h"
dMat file2mat(ifstream& fName, int nVecs)
{	
	dMat X;
	dVec x;
	char			thisChar = ' ';
	string			thisNum;

	
	if (!fName.good()) //file failed to open
	cout << "file no good" << endl;
	else //file opened successfully
	{
		while (!fName.eof())									//as long as lines are left in the file and X.size() < nVecs
		{	
			x.clear();
			thisChar= ' ';
			while (true)										// so long as we are not at end of line or eof()
			{				
				if ( thisChar == '\n' || fName.eof() )	//break if line or file ends
					break;
				else if ( thisChar == ' ' || thisChar == ',' ) //if the character is white space or a comma, get the next character
				{
					fName.get(thisChar);//grab next char from file
					continue;
				}
				else									//if we have part of a number, push a number into x
				{
					thisNum.clear(); //fresh slate
					while (true)//															//as long as we keep getting numbers, append characters to thisChar
					{	
						thisNum.append(1,thisChar);//append thisChar to a string
						fName.get(thisChar);//and get the next character
						if (thisChar==' ' || thisChar==',' ||thisChar=='\n' || fName.eof()) //if we hit a delimiter, or line or file ends, push the number into x
						{
							x.push_back( atof(thisNum.c_str() ) );
							break;
						}						
					}
					
					if (thisChar=='\n' || fName.eof())	//if line or file ended, will break on next while iteration 
						continue;
					else
						fName.get(thisChar); //otherwise resume looking for numbers
				}//end else have a number
				
			}//line or file has ended
			
			X.push_back(x); //accumulate vectors in vector of vectors
				if (X.size()==nVecs)//return if we have specified number of vectors
					return X;
		}//end while X.size()< nVecs and not eof()
	}//else file opened properly
	return X;//have hit the end of the file before X.size()==nVecs, so return all lines from file
};

dMat file2vec(char* fName)
{
	dMat x;
	return x;
};

dVec arr2vec(double arr [], int size)
{
	dVec x;
	for (int i=0; i<size ; i++)
		x.push_back(arr[i]);
		
	return x;
};

void printVec(dVec a)
{
	for (unsigned int i=0; i<a.size(); i++)
	cout << a[i] << "  ";
	
	cout <<endl;
};