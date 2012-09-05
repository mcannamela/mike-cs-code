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
namespace ublas = boost::numeric::ublas;

template<class T>
ublas::matrix<T> fdGradient (const ublas::vector<T>& x,const ublas::vector<T>& dx, boost::function<ublas::vector<T> (ublas::vector<T>)> f) 
{
	using namespace ublas;
	vector<T> F;
	F = f(x);
	int nIn = x.size();
	int nOut = F.size();
	matrix<T> xPlus(nIn,nIn);
	matrix<T> x2Plus(nIn,nIn);
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

 
ublas::matrix<T> dumbInverse(const ublas::matrix<T>& A)
{
		ublas::matrix<T> Ainv(2,2);
		T detA;
		
		detA= A(0,0)*A(1,1)-A(0,1)*A(1,0);
		Ainv(0,0)= A(1,1);
		Ainv(0,1)= -A(0,1);
		Ainv(1,0)= -A(1,0);
		Ainv(1,1)= A(0,0);
		Ainv=Ainv/detA;
		return Ainv;
}

#endif
// typedef ublas::vector<double> uVec;
// typedef ublas::matrix<double> uMat;

// uVec paraboloid(uVec x, int p)
// {
	// uVec y(1);
	// y=0*y;
	// std::cout << "argument passed to paraboloid is " << x << std::endl << std::endl;
    // std::cout << "initial value of y is " << y<< std::endl << std::endl;
	// for(int i=0; i<x.size(); i++)
		// y(0)+=pow(x(i),p);

	// return y;
// };

// int main(void)
 // {
	// using namespace std;
	// uVec x(2);
	// uVec dx(2);
	// uMat G(1,2);
	// ublas::scalar_vector<double> c(2,.1);
	// x(0)=2;
	// x(1)=2;
	// dx = c;
	// boost::function<uVec (uVec)>f;
	// f = boost::bind(paraboloid, _1, 2);

	// G = fdGradient(x,dx,f);
	// cout << G << endl;
	// return 0;
 // };
 
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/matrix_proxy.hpp>
//#include <boost/numeric/ublas/io.hpp>
//


// class button
// {
// public:

    // boost::function<void()> onClick;
// };

// class player
// {
// public:

    // void play();
    // void stop();
// };

// button playButton, stopButton;
// player thePlayer;

// void connect()
// {
    // playButton.onClick = boost::bind(&player::play, &thePlayer);
    // stopButton.onClick = boost::bind(&player::stop, &thePlayer);
// }
