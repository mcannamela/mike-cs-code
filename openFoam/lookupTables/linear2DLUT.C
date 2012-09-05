#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <exception>
#include <cmath>
#include <unistd.h>

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

int main (void){
	std::cout<<"constructing LUT\n";
	linear_two_D_LUT LUT("test_2d_1.LUT");
	std::cout<<"done. printing.\n";
	LUT.printMe();
	/*linear_two_D_LUT LUT2;
	LUT2 = LUT;
	LUT2.printMe();*/
	double x0,x1;
	x0 = .25;
	x1 = 1;
	std::cout << "LUT("<<x0<<","<<x1<<") "<< LUT(x0,x1) << "\n";
	
}


