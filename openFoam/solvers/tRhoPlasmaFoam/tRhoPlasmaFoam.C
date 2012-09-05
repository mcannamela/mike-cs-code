/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    rhoPimpleFoam

Description
    Transient solver for laminar or turbulent flow of compressible fluids
    for HVAC and similar applications.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"
#include "bound.H"
#include "pimpleControl.H"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <exception>
#include <cmath>
#include <unistd.h>
class linear_one_D_LUT {
	private:
		std::vector<double> X;
		
		double tMax;
		double tMin;

	
	public:
		linear_one_D_LUT(){
			tMin = 0;
			tMax = 1;
		}
		
		linear_one_D_LUT(std::string fname);
		
		
		double operator()(double t);
		
};
class linear_enthalpy_LUT_gas {
	public:
		linear_one_D_LUT k, T, c, sigma, mu, h;//,rho;
	
		linear_enthalpy_LUT_gas(std::string gasDir);
		
};


linear_one_D_LUT::linear_one_D_LUT(std::string fname){
			std::ifstream f;
			std::string numStr;
			numStr.clear();
			char ch = 'a';
		  	
			//open the file fname
			f.open(fname.c_str(), std::ifstream::in);
			Info << fname.c_str() << endl;
			Info << f.fail() << endl;
			Info << f.bad() << endl;

			//read in the arrays there, tMin, tMax on the first line, X on the second
			
			//read tMin and tMax
			for(int i=0; i<2; i++){
				//Info<< "reading tMin, tMax"<<endl;
				ch = f.get();
				//Info<< "got a char"<<endl;
				Info << ch;
				while (ch!=' ' && ch!= '\n'){
					numStr+= ch;
					ch = f.get();
					//Info << ch;
				} 
				if(i==0)
					tMin = atof(numStr.c_str());
				if(i==1)
					tMax = atof(numStr.c_str());
				numStr.clear();
				
			}
			
			//burn to next line
			while(ch!= '\n')
				ch = f.get();
				
			
			//read X
			ch = f.get();
			while (true){
					
					//if ch is a space or newline add an element to X
					if(ch==' ' || ch=='\n' || f.eof()){
						X.push_back(atof(numStr.c_str()));
						numStr.clear();
						
						//escape the while loop when we hit a newline or finish the file
						if(ch=='\n' || f.eof())
							break;
					}
					
					//append to numStr if we haven't completed an entry
					if(ch!=' ')
						numStr+= ch;
					ch = f.get();				
			}
			
			f.close();

}
		
		
double linear_one_D_LUT::operator()(double t){
			double ii = (X.size()-1)*(t-tMin)/(tMax-tMin);
			int i = floor(ii);
			double f = ii-i;
			
			try{
				return X[i]+f*(X[i+1]-X[i]);
			}
			catch(std::exception& e){
				if(i<0)
					return X[0];
			    else if(i>(X.size()-1))
					return X[X.size()-1];
				else 
					Info << "t: "<< t<< "   i: "<< i <<"\nbad index into LUT!" << endl;
					return 0x7ff0000000000000;
					
			}
			
		}
class linear_two_D_LUT {
	private:
		std::vector< std::vector<double> >  X;
		std::ifstream f;
		std::string fname;
		
		double t0Min;
		double t0Max;
		double t1Min;
		double t1Max;
		
		
	
		double ii,jj, f0,f1;
		double x00,x01,x10,x11;
		double x0,x1;
		int i_;
		int j_;		
		int k_;
	public:
		linear_two_D_LUT(){
			t0Min = 0;
			t0Max = 1;
			t1Min = 0;
			t1Max = 1;	
		}
		linear_two_D_LUT(std::string fname);
		linear_two_D_LUT(const linear_two_D_LUT& other);
       ~linear_two_D_LUT(void);
		linear_two_D_LUT& operator=(const linear_two_D_LUT& other);
		
		double operator()(double t0, double t1);
		
		double idx0(double t){
			return std::max( std::min( (size0()-1)*(t-t0Min)/(t0Max-t0Min), double(size0())-1.000001), 0.0);
		}
		double idx1(double t){
			return std::max( std::min( (size1()-1)*(t-t1Min)/(t1Max-t1Min), double(size1())-1.000001), 0.0);
		}
		
		int size0(void){
			return X[0].size();
		}
		int size1(void){
			return X.size();
		}
		
		
		double getNum(void);
		void printMe(void);
};

linear_two_D_LUT::~linear_two_D_LUT(void){
	if (f.is_open()){	
		f.close();
}
}
linear_two_D_LUT::linear_two_D_LUT(const linear_two_D_LUT& other){

	t0Min = other.t0Min;
	t0Max = other.t0Max;
	t1Min = other.t1Min;
	t1Max = other.t1Max;
	X = other.X;
	fname = other.fname;

}
linear_two_D_LUT& linear_two_D_LUT::operator=(const linear_two_D_LUT& other){
	t0Min = other.t0Min;
	t0Max = other.t0Max;
	t1Min = other.t1Min;
	t1Max = other.t1Max;
	X = other.X;
	fname = other.fname;
	return *this;
}

linear_two_D_LUT::linear_two_D_LUT(std::string fname){
	int n0, n1;
	std::vector<double> xxtmp;
	double x;
	
	//open the file fname
	f.open(fname.c_str(), std::ifstream::in);
    std::cout << "reading LUT from file "<<fname.c_str() << std::endl;
	std::cout << "fail is" << f.fail() << std::endl;
	std::cout << "bad is" << f.bad() << std::endl;
	
	if (f.fail()==1  || f.bad()==1)
		throw;

	//read in the arrays there, tMin, tMax on the first line, X on the second
	
	//read tMin and tMax
	t0Min = getNum();
	t0Max = getNum();
	n0    = getNum();
	t1Min = getNum();
	t1Max = getNum();
	n1    = getNum();

	

	for (int j=0;j<n1;j++){
		xxtmp.clear();
		for (int i=0; i<n0; i++){
			x = getNum();
			xxtmp.push_back(x);
		}
		X.push_back(xxtmp);
	}
	f.close();
}

void linear_two_D_LUT::printMe(void){
	for (int j=0;j<size1();j++){	
		for (int i=0; i<size0(); i++){
			std::cout << X[j][i] << "  ";
		}
		std::cout << std::endl;
	}

}
double linear_two_D_LUT::getNum(void){
	std::string numStr;
	numStr.clear();
	char ch = 'a';
	double num;
	
	ch = f.get();
	////std::cout << "got a ch " << ch << "\n";
	if (ch==' ' || ch== '\n'){
		while (true){
			ch = f.get();
			////std::cout << "ch is " << ch << "\n";
			if (ch!=' ' && ch!= '\n'){
				break;
			}
		}
	}
	////std::cout << "burned, will now get the num" << "\n";
	while (ch!=' ' && ch!= '\n'){
		numStr+= ch;
		ch = f.get();
	} 
	num = atof(numStr.c_str());
	////std::cout<< "got it, the num is " << num << "\n";
	numStr.clear();
	return num;
}

double linear_two_D_LUT::operator()(double t0, double t1){

			ii = idx0(t0);
			jj = idx1(t1);
			
			i_ = floor(ii);
			j_ = floor(jj);

			f0 = ii-i_;
			f1 = jj-j_;
			
			x00 = X[j_][i_];
			x01 = X[j_+1][i_];
			x10 = X[j_][i_+1];
			x11 = X[j_+1][i_+1];
			
			
			x0 = x00+f1*(x01-x00);
			x1 = x10+f1*(x11-x10);
			

			return x0+f0*(x1-x0);
			
		}		


linear_enthalpy_LUT_gas::linear_enthalpy_LUT_gas(std::string gasDir){
			
		std::string kName(gasDir), cName(gasDir), TName(gasDir), sName(gasDir), mName(gasDir), hName(gasDir);//,rName(gasDir);
		
		Info <<"constructing linear_enthalpy_LUT_gas from directory : " << gasDir <<endl;
					
		kName+="thermalConductivity.LUT";
		sName+="electricalConductivity.LUT";
		TName+="temperature.LUT";
		mName+="viscosity.LUT";
		//rName+="density.LUT";
		cName+="specificHeat.LUT";
		hName+="enthalpy.LUT";
		
		Info <<" filenames set, constructing LUTs " << endl;
		
		linear_one_D_LUT k_(kName);
		linear_one_D_LUT c_(cName);
		linear_one_D_LUT T_(TName);
		linear_one_D_LUT sigma_(sName);
		linear_one_D_LUT mu_(mName);
		//linear_one_D_LUT rho_(rName);
		linear_one_D_LUT enthalpy_(hName);
		
		Info << "complete! " << endl;
		
		
		k = k_;
		T = T_;
		c = c_;
		mu = mu_;
		//rho = rho_;
		sigma = sigma_;
		h = enthalpy_;
			
}
double radiatiativeLoss(double T)
{
		if (T<9500)
			return 0;
		else
			return 5600*(T-9500)+181*pow((T-9500),2);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readMHDControls.H"
    #include "createFields.H"
    Info<< "initializing properties...\n" << endl;
    #include "updateProperties.H"
    Info<< "...complete. will update thermo...\n" << endl;
    #include "updateThermo.H"
    
  //  rhoPlasma = rho*f;
    //rhoAmbient = rho*(1-f);
    
    Info<< "done\n" << endl;
    #include "initContinuityErrs.H"

    pimpleControl pimple(mesh);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "rhoEqn.H"
        #include "updateTurbulence.H"
       /* forAll(f.internalField(), celli)
		{
			if (f.internalField()[celli]>1.0)
				f.internalField()[celli] = 1.0;
		}*/

        // --- Pressure-velocity PIMPLE corrector loop
        //for (pimple.start(); pimple.loop(); pimple++)
        while (pimple.loop())
        {
            /*if (pimple.nOuterCorr() != 1)
            {
                p.storePrevIter();
                rho.storePrevIter();
            }*/

            #include "UEqn.H"
            #include "hEqn.H"
            // --- PISO loop
            //for (int corr=0; corr<pimple.nCorr(); corr++)
            while (pimple.correct())
            {
				
                #include "pEqn.H"
            }
            if(mixingOn){
				#include "plasmaFractionEqn.H"
		    }
		    Info << "updating turbulence...\n";
		    if(turbulenceOn){
				turbulence->correct();
				#include "updateTurbulence.H"
		    }
            
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
