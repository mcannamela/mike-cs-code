#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <exception>
#include <cmath>
using namespace std;
typedef vector<double> dvec;

class linear_one_D_LUT {
	public:
		linear_one_D_LUT(){
			tMin = 0;
			tMax = 1;
			
		}
		linear_one_D_LUT(string fname){
			ifstream f;
			string numStr;
			numStr.clear();
			char ch = 'a';
			
			//open the file fname
			f.open(fname.c_str());
			cout << fname.c_str()<< endl;
			//read in the arrays there, tMin, tMax on the first line, X on the second
			
			//read tMin and tMax
			for(int i=0; i<2; i++){
				//cout<< "reading tMin, tMax"<<endl;
				ch = f.get();
				//cout << ch;
				while (ch!=' ' && ch!= '\n'){
					numStr+= ch;
					ch = f.get();
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
		
		
		double operator()(double t){
			uint i = round((X.size()-1)*(t-tMin)/(tMax-tMin));
			
			try{
				return X[i];
			}
			catch(exception& e){
				if(i<0)
					return X[0];
			    else if(i>(X.size()-1))
					return X[X.size()-1];
				else 
					cout << "t: "<< t<< "   i: "<< i <<"\nbad index into LUT!" << endl;
					return 0x7ff0000000000000;
					
			}
			
		}
		
	
	private:
		dvec X;
		
		double tMax;
		double tMin;
		

		
		
};
class linear_enthalpy_LUT_gas{
	public:
		linear_one_D_LUT k, T, c, sigma, mu, rho;
		
		
		linear_enthalpy_LUT_gas(string gasDir){
			
			string kName(gasDir), cName(gasDir), TName(gasDir), sName(gasDir), mName(gasDir), rName(gasDir);
						
			kName+="thermalConductivity.LUT";
			sName+="electricalConductivity.LUT";
			TName+="temperature.LUT";
			mName+="viscosity.LUT";
			rName+="density.LUT";
			cName+="specificHeat.LUT";
			
			linear_one_D_LUT k_(kName);
			linear_one_D_LUT c_(cName);
			linear_one_D_LUT T_(TName);
			linear_one_D_LUT sigma_(sName);
			linear_one_D_LUT mu_(mName);
			linear_one_D_LUT rho_(rName);
			
			
			k = k_;
			T = T_;
			c = c_;
			mu = mu_;
			rho = rho_;
			sigma = sigma_;
			
		}
		
};


int main(void){
	//double h = 1.e8;
	double h = 3.1e5;
	
	cout << "linear LUT test:" << endl;
	linear_one_D_LUT lut("testLinearOneDTable.lut");
	cout << lut(0) << "  " << lut(1.2e-3) << "  " <<lut(1) << "  " << lut(2.1)<<endl;
	cout << "end linear LUT test" <<endl;
	
	cout <<"gas test gas"<< endl;
	linear_enthalpy_LUT_gas gas("./air/");
	cout << "enthalpy is :" << h<< endl;
	cout << "T = "<<gas.T(h) << endl;
	cout << "thermal cond = "<<gas.k(h) << endl;
	cout << "electrical cond = "<<gas.sigma(h) << endl;
	cout << "rho = "<<gas.rho(h) << endl;
	cout << "specificHeat = "<< gas.c(h) << endl;
	cout << "viscosity = "<< gas.mu(h) << endl;
	
	
	//cout << "hello world \n";
	return 0;
};
	
