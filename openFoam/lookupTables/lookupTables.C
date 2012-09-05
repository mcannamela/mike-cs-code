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

class linear_three_D_LUT {
	private:
		std::vector<std::vector<std::vector<double> > > X;
		std::ifstream f;
		
		double t0Min;
		double t0Max;
		double t1Min;
		double t1Max;
		double t2Min;
		double t2Max;
	public:
		linear_three_D_LUT(){
			t0Min = 0;
			t0Max = 1;
			t1Min = 0;
			t1Max = 1;
			t2Min = 0;
			t2Max = 1;
		}
		linear_three_D_LUT(std::string fname);
		
		double operator()(double t0, double t1, double t2);
		
		double idx0(double t){
			return (size0()-1)*(t-t0Min)/(t0Max-t0Min);
		}
		double idx1(double t){
			return (size1()-1)*(t-t1Min)/(t1Max-t1Min);
		}
		double idx2(double t){
			return (size2()-1)*(t-t2Min)/(t2Max-t2Min);
		}
		
		int size0(void){
			return X[0][0].size();
		}
		int size1(void){
			return X[0].size();
		}
		int size2(void){
			return X.size();
		}
		
		double getNum(void);
		void printMe(void);
};
void linear_three_D_LUT::printMe(void){
	for (int k=0;k<size2();k++){
		std::cout << std::endl;
		for (int j=0;j<size1();j++){	
			for (int i=0; i<size0(); i++){
				std::cout << X[k][j][i] << "  ";
			}
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;
}
double linear_three_D_LUT::getNum(void){
	std::string numStr;
	numStr.clear();
	char ch = 'a';
	double num;
	
	ch = f.get();
	//std::cout << "got a ch " << ch << "\n";
	if (ch==' ' || ch== '\n'){
		while (true){
			ch = f.get();
			//std::cout << "ch is " << ch << "\n";
			if (ch!=' ' && ch!= '\n'){
				break;
			}
		}
	}
	//std::cout << "burned, will now get the num" << "\n";
	while (ch!=' ' && ch!= '\n'){
		numStr+= ch;
		ch = f.get();
	} 
	num = atof(numStr.c_str());
	//std::cout<< "got it, the num is " << num << "\n";
	numStr.clear();
	return num;
}
linear_three_D_LUT::linear_three_D_LUT(std::string fname){
	int n0, n1, n2;
	std::vector<std::vector<double> > Xtmp;
	std::vector<double> xxtmp;
	double x;
	
	//open the file fname
	f.open(fname.c_str(), std::ifstream::in);
	std::cout << fname.c_str() << std::endl;
	std::cout << "fail is" << f.fail() << std::endl;
	std::cout << "bad is" << f.bad() << std::endl;

	//read in the arrays there, tMin, tMax on the first line, X on the second
	
	//read tMin and tMax
	t0Min = getNum();
	t0Max = getNum();
	n0    = getNum();
	t1Min = getNum();
	t1Max = getNum();
	n1    = getNum();
	t2Min = getNum();
	t2Max = getNum();
	n2    = getNum();
	
	for (int k=0;k<n2;k++){
		Xtmp.clear();
		for (int j=0;j<n1;j++){
			xxtmp.clear();
			for (int i=0; i<n0; i++){
				x = getNum();
				xxtmp.push_back(x);
			}
			Xtmp.push_back(xxtmp);
		}
		X.push_back(Xtmp);
	}
	
	f.close();
}
double linear_three_D_LUT::operator()(double t0, double t1, double t2){
			double ii,jj,kk, f0,f1,f2;
			double x000,x010,x100,x110,x001,x011,x101,x111;
			double x00,x01,x10,x11;
			double x0,x1;
			int i,j,k;
			
			ii = idx0(t0);
			jj = idx1(t1);
			kk = idx2(t2);
			
			i = floor(ii);
			j = floor(jj);
			k = floor(kk);
			f0 = ii-i;
			f1 = jj-j;
			f2 = kk-k;
			
			x000 = X[k][j][i];
			x010 = X[k][j+1][i];
			x100 = X[k][j][i+1];
			x110 = X[k][j+1][i+1];
			
			x001 = X[k+1][j][i];
			x011 = X[k+1][j+1][i];
			x101 = X[k+1][j][i+1];
			x111 = X[k+1][j+1][i+1];
			
			x00 = x000+f2*(x001-x000);
			x01 = x010+f2*(x011-x010);
			x10 = x100+f2*(x101-x100);
			x11 = x110+f2*(x111-x110);
			
			x0 = x00+f1*(x01-x00);
			x1 = x10+f1*(x11-x10);
			
			std::cout<< "f's: " <<f0<<" " << f1 << " "<<f2 <<"\n";
			std::cout<< "ii " <<ii<<" jj " << jj << " kk "<<kk<<"\n";
			std::cout<< " x000 " << x000 
			         << " x010 " << x010
			         << " x100 " << x100
			         << " x110 " << x110 << "\n"
			         << " x001 " << x001 
			         << " x011 " << x011
			         << " x101 " << x101
			         << " x111 " << x111 << "\n"
			         << " x00 "  << x00  
			         << " x01 "  << x01
			         << " x10 "  << x10
			         << " x11 "  << x11  << "\n"
			         << " x0 "   << x0
			         << " x1 "   << x1
			         << "\n";
			         
			
			return x0+f0*(x1-x0);
		}

linear_one_D_LUT::linear_one_D_LUT(std::string fname){
			std::ifstream f;
			std::string numStr;
			numStr.clear();
			char ch = 'a';
		  	
			//open the file fname
			f.open(fname.c_str(), std::ifstream::in);
			std::cout << fname.c_str() << std::endl;
			std::cout << f.fail() << std::endl;
			std::cout << f.bad() << std::endl;

			//read in the arrays there, tMin, tMax on the first line, X on the second
			
			//read tMin and tMax
			for(int i=0; i<2; i++){
				//std::cout<< "reading tMin, tMax"<<std::endl;
				ch = f.get();
				//std::cout<< "got a char"<<std::endl;
				std::cout << ch;
				while (ch!=' ' && ch!= '\n'){
					numStr+= ch;
					ch = f.get();
					//std::cout << ch;
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
					//std::cout << "t: "<< t<< "   i: "<< i <<"\nbad index into LUT!" << std::endl;
					return 0x7ff0000000000000;
					
			}
			
		}
		
		
int main(void){
	/*std::string LUTFile("test.lut");
	linear_three_D_LUT L3(LUTFile);
	L3.printMe();*/
	
	/*std::string LUTFile("test2.lut");
	linear_three_D_LUT L3(LUTFile);
	L3.printMe();
	std::cout << L3(.5, 1.5, 2.5) << "\n";*/
	
	std::string LUTFile("test3.lut");
	linear_three_D_LUT L3(LUTFile);
	L3.printMe();
	std::cout << L3(1.5, 2.5, 3.5) << "\n";
	std::cout << "size of X: " << L3.size0() << " "<<L3.size1()<< " "<<L3.size2()<<"\n";
	return 0;
}
