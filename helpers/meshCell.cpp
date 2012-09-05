#include "meshCell.h"

meshCell::meshCell();
{
	int faceIdx  = 0;
	int innerIdx = 1 ;
	int outerIdx = 2 ;
	r.resize(2,0);
	z.resize(2,0);
};