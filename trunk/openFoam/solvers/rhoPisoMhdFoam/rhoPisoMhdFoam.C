/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    rhoPisoFoam

Description
    Transient PISO solver for compressible, laminar or turbulent flow.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
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
			Info << "fail is "<< f.fail() 
			<< ", bad is "<<f.bad() << endl;

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
	//from bernardi, colombo, and ghedini
	//comparison of different techniques for the fluent based 
	//treatment of the electromagnetic field in inductively coupled plasma torches
		if (T<9500)
			return 0;
		else
			return 5600*(T-9500)+181*pow((T-9500),2);
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    
    #include "readMHDControls.H"
    #include "createPropertiesFields.H"
    #include "createFields.H"
    
    Info<< "initializing properties...\n" << endl;
    #include "updateProperties.H"
    #include "updateDoping.H"
	#include "updateElectrodeBCs.H"
    Info<< "done\n" << endl;
    
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    
    
    


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" 
    << "maxDeltaT is " << maxDeltaT
    << endl;
    
    turbulence->correct();
    #include "updateTurbulence.H"
    #include "pheeEqn.H"
    
    forAll(rho.internalField(), celli)
		{initialMass+= rho.internalField()[celli]*mesh.V()[celli];}
	reduce(initialMass, sumOp<double>());

    while (runTime.run())
    {
        Info <<"\n********************************************************************\n";
		Info <<"********************************************************************\n";
        #include "readTimeControls.H"
        #include "readPISOControls.H"
        #include "readMHDControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        Info <<"********************************************************************\n\n";
        turbulence->correct();
        

        for (int bigcorr = 1; bigcorr<=nOuterCorr; bigcorr++)
        {
			Info<< "\n }}}}}}}}}}}}}}}}}}}OUTER CORRECTION "<< bigcorr << " BEGIN!{{{{{{{{{{{{{{\n\n";
			{
				// Info<< "\n ----------SOLVING CONTINUITY ------------\n\n";
				 solve(fvm::ddt(rho) + fvc::div(phi));
			}
			
			//Predict U
			Info<< "\n >>>>>>>>>>> PREDICTING U <<<<<<<<<<<<<\n";
			#include "UEqn.H"
			Info<< "------------------------------------------------------------------------------------\n";

			// energy equation and pressure correction loop
			pMax = MHDControls.lookup("pMax");
			pMin = MHDControls.lookup("pMin");
			radiationOn = readBool(MHDControls.lookup("radiationOn"));
			for (int corr=1; corr<=nCorr; corr++)
			{
				
				Info<< "\n		:::::::::energy and pressure correction loop iteration "<<corr<< " ::::::::::\n";
				#include "hEqn.H"
			   // Info<< "Done. Now solving for p\n";
				#include "pEqn.H"
			 //   Info<< "Done.\n";
			  //  #include "pheeEqn.H"
			}
			
			Info << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
			Info << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
			Info << "patch mass flux:"<< "  ";
			double mDotTotal = 0;
			surplus = 0;
			forAll(phi.boundaryField(), patchi)
			{
				mDotTotal=0;
				forAll(phi.boundaryField()[patchi], facei)
				{
					mDotTotal+=phi.boundaryField()[patchi][facei];
				}
				reduce(mDotTotal, sumOp<double>());
				Info << mesh.boundary()[patchi].name() << " "<< mDotTotal << ", ";
				if (mesh.boundary()[patchi].name()=="inlet")
					surplus+=mag(mDotTotal)*runTime.deltaT().value();
				if (mesh.boundary()[patchi].name()=="outlet")
					surplus-=mag(mDotTotal)*runTime.deltaT().value();
			}
			Info <<endl;
			double totalMass = 0;
			forAll(rho.internalField(), celli)
				{totalMass+=rho.internalField()[celli]*mesh.V()[celli];}
			reduce(totalMass, sumOp<double>());
			Info << "total mass in volume is "<< totalMass <<"\n";
			Info << "initial mass at start of run was "<< initialMass<<"\n";
			Info << "so the surplus is "<< totalMass-initialMass << " and should be "<< massSurplus+surplus <<"\n";
			Info << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
			
			nBcorr = readInt(MHDControls.lookup("nCorrectors"));
			EM_scheme = readInt(MHDControls.lookup("EM_scheme"));
			J_scheme = readInt(MHDControls.lookup("J_scheme"));
			jouleHeatingMax = MHDControls.lookup("jouleHeatingMax");
			cathodeJouleHeatingMax = MHDControls.lookup("cathodeJouleHeatingMax");
			lorenzForceOn = readBool(MHDControls.lookup("lorenzForceOn"));
			for (int Bcorr=1; Bcorr<=nBcorr; Bcorr++)
			{
				Info<< "\n///////////////////////// EM loop iteration "<<Bcorr<< " /////////////////////////\n";
				
				#include "pheeEqn.H"
				#include "AEqn.H"
				Info<< "/////////////////////////////////////////////////////////////////////////\n\n";
			}
			
			
			//solve turbulence model and update turbulent transport properties
			turbulence->correct();

			//??????????????????????????????????????????????????????????????????????????????????????????????????????????????????
			//??????update density by the thermochemical model. why do we bother computing the density with rhoEqn.H, then??????
			rho = thermo.rho();
			//??????????????????????????????????????????????????????????????????????????????????????????????????????????????????
			Info<< "------------------------------------------------------------------------------------\n";
            Info<< "------------------------------------------------------------------------------------\n";
     	}
     	massSurplus+=surplus;
        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
