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
		linear_one_D_LUT k, T, c, sigma, mu, h;
	
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
