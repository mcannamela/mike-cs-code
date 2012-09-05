#ifndef FILEHELPERS
#define FILEHELPERS

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

using namespace std;
typedef vector<double> 	dVec;
typedef vector<dVec>	dMat;

//void vec2file(char* fName, dVec x)
//{
//};
//void vec2file(char* fName, dMat x)
//{
//};
//void vec2file(char* fName, dVec x, char modeFlag)
//{
//};
//void vec2file(char* fName, dMat x, char modeFlag)
//{
//};
dMat file2mat(ifstream& , int );
dMat file2vec(char* );
dVec arr2vec(double [] , int );
void printVec(dVec );

#endif
//int main(void)
//{
//ifstream aFile;
//dVec 	 x;
//aFile.open("testTable.txt");
//x = file2vec(aFile);
//cout << "x ";
//printVec(x);
//return 0;
//};