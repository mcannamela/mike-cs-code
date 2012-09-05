#include "fileHelpers.h"

int main(void)
{
ifstream aFile;
dMat 	 X;
int		nVecs= 3;
aFile.open("testTable.txt");
X = file2mat(aFile,nVecs);
cout << "X ";
for(int i=0; i<X.size(); i++)
printVec(X[i]);

return 0;
};