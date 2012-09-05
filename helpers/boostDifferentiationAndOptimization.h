#ifndef UBLAS_FD_GRADIENT_H
#define UBLAS_FD_GRADIENT_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <iostream>
#include <cmath>

#ifndef UVEC_TYPES 
#define UVEC_TYPES
namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> uVec;
 typedef ublas::matrix<double> uMat;
#endif
template<class T>
ublas::matrix<T> fdGradient (const ublas::vector<T>& x,const ublas::vector<T>& dx, boost::function<ublas::vector<T> (ublas::vector<T>)> f) 
{
	using namespace ublas;
	vector<T> F;
	F = f(x);
	int nIn = x.size();
	int nOut = F.size();
	matrix<T> fPlus(nIn,nOut);
	matrix<T> f2Plus(nIn,nOut);
	matrix<T> G(nOut,nIn);
	
	
	for(int i=0;i<nIn;i++)
	{
		unit_vector<T> e(nIn,i); 
		row(fPlus, i)= f(x + element_prod(e,dx) );     // assign the 2x3 element submatrix of A
		row(f2Plus, i)= f(x + 2*element_prod(e,dx) );
		std::cout << "i is   " << i << std::endl;
		std::cout << "e  " << e << std::endl;
		std::cout << "x+e.dx  " << x + element_prod(e,dx)  << std::endl;
		std::cout << "f(x+e.dx)  " << row(fPlus, i) << std::endl <<std::endl;
	}
	for(int i=0;i<nIn;i++)
	{
		for(int j=0;j<nOut;j++)
		{
			G(j,i) = -(f2Plus(i,j) - 4*fPlus(i,j) + 3*F(j))/(2*dx(i));
		}
	}
	return G;
}

//non templated version of fdGradient, written as a workaround to binding problems with fdGradient
//TODO: make binding work with templates
uMat fdGradientD (const uVec& x,const uVec& dx, boost::function<uMat (uVec)> f) 
{
//	std::cout<<"x is " << x <<std::endl;
//	std::cout<<"dx is " << dx <<std::endl;
	using namespace ublas;
	uVec F;//the function value at x
	F = row(f(x),0);
//	std::cout << "f(x) " << F << std::endl;
	int nIn = x.size();//length of vector into function
	int nOut = F.size();//length of vector out of function
	bool blah = false;//flag for debugging, can ignore if not debugging
	matrix<double> fPlus(nIn,nOut);//function value at x+deltaX
	matrix<double> f2Plus(nIn,nOut);//function value at x+2deltaX
	matrix<double> G(nOut,nIn);//holds the gradient (jacobian of f)
	
	
	for(int i=0;i<nIn;i++) //compute function at stencil points x, x+dx, x+2dx for each x(i)
	{
		unit_vector<double> e(nIn,i); 
		row(fPlus, i)= row( f(x + element_prod(e,dx)),0  );     
		row(f2Plus, i)= row(f(x + 2*element_prod(e,dx)),0 );
		//std::cout << "i is   " << i << std::endl;
		//std::cout << "e  " << e << std::endl;
	//	std::cout << "x+e.dx  " << x + element_prod(e,dx)  << std::endl;
	//	std::cout << "f(x+e.dx)  " << row(fPlus, i) << std::endl; 
	}
	
	//compute members of the jacobian using finite differences
	for(int i=0;i<nIn;i++)
	{
		for(int j=0;j<nOut;j++)
		{
			G(j,i) = -(f2Plus(i,j) - 4*fPlus(i,j) + 3*F(j))/(2*dx(i));
			/*if (nOut == 2)
			{
			std::cout << "G(" << j << "," << i <<") = " << G(j,i) << std::endl;	
			}*/
		}
	}
	
//	std::cout << "G(x)  " << G << std::endl << std::endl;
	return G;
}

/* Matrix inversion routine.
    Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
 template<class T>
bool InvertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
 	using namespace boost::numeric::ublas;
 	typedef permutation_matrix<std::size_t> pmatrix;
 	// create a working copy of the input
 	matrix<T> A(input);
 	// create a permutation matrix for the LU-factorization
 	pmatrix pm(A.size1());

 	// perform LU-factorization
 	int res = lu_factorize(A,pm);
        if( res != 0 ) return false;

 	// create identity matrix of "inverse"
 	inverse.assign(ublas::identity_matrix<T>(A.size1()));

 	// backsubstitute to get the inverse
 	lu_substitute(A, pm, inverse);

 	return true;
 }

 template<class T>
ublas::matrix<T> dumbInverse(const ublas::matrix<T>& A)
{
		ublas::matrix<T> Ainv(2,2);
		T detA;
		//std::cout << A << "\n";
		detA= A(0,0)*A(1,1)-A(0,1)*A(1,0);
		Ainv(0,0)= A(1,1);
		Ainv(0,1)= -A(0,1);
		Ainv(1,0)= -A(1,0);
		Ainv(1,1)= A(0,0);
		Ainv=Ainv/detA;
		return Ainv;
}


uVec rowMat2Vec(const uMat& A)
{
	uVec v(A.size2());
	v= row(A,1);
	return v;
}

uMat paraboloid(uVec x, double p )
{
	uMat normX(1,1);
	normX(0,0)=0;	
	for(int i=0; i<x.size();i++)
		{
			normX(0,0)+=pow(x(i),p);
		}
	return normX;
}

//perform one newton iteration on f, return the new x in a vector, put the function value into the reference F 
uVec scalarNewtonIterate(const uVec& x,const uVec& dx, boost::function<uMat (uVec)> f, uMat& F)
{
	bool verbose = false;//false suppresses output statements, turn off for speed
	bool verifyScalar = false;//check to make sure f is scalar valued (1x1 matrix), turn off for speed
	uVec xNew(x.size());//new guess goes here for returning
	uVec deltaX; //subtract from x to get xNew
	uVec J(x.size());//to hold the Jacobian 
	uMat H(x.size(),x.size());//to hold the Hessian
	uMat invH(x.size(),x.size());// to hold the inverse of the Hessian
	boost::function<uMat (uVec)>Jacobian; //wrapper for first derivative (a vector)
	boost::function<uMat (uVec)>Hessian; //wrapper for second derivative (a scalar)
	
	//function compositing, not the most efficient considering the hessian is symmetric
	Jacobian = boost::bind<uMat>(fdGradientD, _1,dx,f);
	Hessian  = boost::bind<uMat>(fdGradientD,_1,dx,Jacobian);
	
	if (verifyScalar)
	{
		if ( !( (f(x).size1() == 1) && (f(x).size2() == 1) ) )
		{
			std::cout << "error in ScalarNewtonIterate: higher order tensor detected, \n		this function only optimizes scalar functions. \n	f(x) must be 1x1 \n" ;
		};
	};

	
	J		=	row(Jacobian(x),0);//compute Jacobian
	H		=	Hessian(x);//compute Hessian
	
	//need inverse of Hessian now
	if (  (H.size1()==2) && (H.size2()==2)   )//use a dumb but reliable expression for the inverse if H is 2x2, else use LU decomposition 
	invH = dumbInverse(H);// 1/det(H) sort of thing
	else
	{
		InvertMatrix(H,invH);// LU decomposition inversion
	//	std::cout << "\n" << "inverse of Hessian is " << "\n" << "		" << invH <<"\n" << "		";
		//std::cout << "\n" << "H*invH = " << "\n" << "		" << prod(H,invH) <<"\n" << "		";
	}
	//get the update to x
	deltaX	=	prod(invH,J);//want to find the stationary point of the Jacobian
	
	//display some things if we care to know them
	if (verbose)
	{
		std::cout << "the function at x = " << x << " is " <<  row(f(x),0) << "\n";
		std::cout << "\n" << "Jacobian is " << J << "\n";
		std::cout << "\n" << "Hessian is " << "\n" << "		" << H <<"\n" << "		";
		std::cout << "\n" << "next step in x is " << -deltaX << "\n";
	}
	
	//returning values
	xNew = x - deltaX;//update x
	F = f(xNew);//put the value of the function into F
	return xNew;
} 


#endif
