#include "cylMesh.h"

//constructors
	void cylMesh::cylMesh()
	{
	nullCell  		= cylMeshCell();
	nullCellVec.push_back(nullCell); 
	zIdx=0;
	rIdx=0;
	addCell('z');
	dZdefault = 1;
	dRdefault = 1;
	};
	
	void cylMesh::cylMesh(double dZ, double dR)
	{
	nullCell  		= cylMeshCell();
	nullCellVec.push_back(nullCell); 
	zIdx=0;
	rIdx=0;
	addCell('z');
	
	dZdefault = dZ;
	dRdefault = dR;
	};
	
	//mesh expansion
	void cylMesh::addCell(char c)//add in a direction using default  or previous step size
		{
			switch(c)
			{
				case 'r': 	int 		endIdx  = cells[zIdx].size()-1;
							cylMeshCell endCell = cells[zIdx][endIdx];//formerly the last radial cell in the mesh
																					
							cells[zIdx].push_back(endCell); //new radial cell is just like the previous one...
							cells[zIdx][endIdx+1].set('r', endCell.R() + endCell.dR() ); //...plus its  radial step size
							break 
							
				case 'z':	double  zVal   = 0;//will be replaced by the new zVal unless this is the first cell
							double  dzVal  = dZdefault;//will use previous z step unless this is the first cell
							int 	endIdx = cells.size();
							cells.push_back(nullCellVec);//make room for new cell with an empty cell
							if(cells.size()!=1) //as long as this is not the first cell
							{
								dzVal = cells[endIdx][0].dZ(); //use the previous axial step size
								zVal = cells[endIdx][0].Z() + cells[endIdx][0].dZ(); //and the previous axial value plus the step size
							}
							cells[endIdx+1][0].set('z',zVal);
							cells[endIdx+1][0].setD('z',dzVal);
							cells[endIdx+1][0].set('r',0);
							cells[endIdx+1][0].set('r',dRdefault);
			}
		};
	void cylMesh::addCell(char c, double stepSize)// add a cell in a direction using a specified step size
		{
			switch(c)
			{
				case 'r': 	int 		endIdx  = cells[zIdx].size()-1;
							cylMeshCell endCell = cells[zIdx][endIdx];//formerly the last radial cell in the mesh
																					
							cells[zIdx].push_back(endCell); //new radial cell is just like the previous one...
							cells[zIdx][endIdx+1].set('r', endCell.R() + endCell.dR() ); //...plus its  radial step size
							cells[zIdx][endIdx+1].setD('r', stepSize ); //set the radial step size of the new cell
							break
							
				case 'z':	double  zVal   = 0;//will be replaced by the new zVal unless this is the first cell
							double	drVal  = dRdefault;//set to default in case this is the first cell
							int 	endIdx = cells.size();
							cells.push_back(nullCellVec);//make room for new cell with an empty cell
							if(cells.size()!=1) //as long as this is not the first cell
							{
								zVal = cells[endIdx][0].Z() + cells[endIdx][0].dZ(); //and the previous axial value plus the step size
								drVal = cells[endIdx][0].dR();//use the neighboring radial step size
							}
							cells[endIdx+1][0].set('z',zVal);
							cells[endIdx+1][0].setD('z',stepSize);
							cells[endIdx+1][0].set('r',0);
							cells[endIdx+1][0].set('r',drVal);
			}
		};