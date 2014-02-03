//
//  Chebyshev_Interpolation_1D.cpp
//  
//
//  Created by Sivaram Ambikasaran on 1/2/14.
//
//

#include <cmath>
#include "Chebyshev_Interpolation_1D.hpp"

/********************************************************************************/
//      FUNCTION:               kernel1D                                        //
//                                                                              //
//      PURPOSE OF EXISTENCE:   Computes a kernel in 1D.                        //
//                                                                              //
//      PARAMETERS:                                                             //
//      x1      -       'x' location of the points in the first cluster.        //
//      n1      -       Number of points in the first cluster.                  //
//      x2      -       'x' location of the points in the second cluster.       //
//      n2      -       Number of points in the second cluster.                 //
//      K       -       Matrix with 'n1' rows and 'n2' columns, where K(n2*i+j) //
//                      is the interaction between the 'i'th point in the first //
//                      cluster and 'j'th point in the second cluster.          //
//                                                                              //
/********************************************************************************/
void kernel1D(double* x1, unsigned n1, double* x2, unsigned n2, double*& K) {
        double Rsquare;
        K       =       new double [n1*n2];
        unsigned index;
        for (unsigned i=0; i<n1; ++i) {
                index   =       i*n2;
                for (unsigned j=0; j<n2; ++j) {
                        Rsquare         =       (x1[i]-x2[j])*(x1[i]-x2[j]);
                        K[index+j]      =       1.0/(Rsquare);
                }
        }
}

/********************************************************************************/
//      FUNCTION:               Chebyshev_polynomials                           //
//                                                                              //
//	PURPOSE OF EXISTENCE:   Evaluates the first rank Chebyshev polynomials  //
//                              at the locations x in [-1,1].                   //
//                                                                              //
//      PARAMETERS:                                                             //
//      rank    -       Number of terms in the approximation.                   //
//      x       -       Location of points in the interval [-1,1].              //
//      n       -       Number of points.                                       //
//      T       -       Matrix with 'n' rows and 'rank' columns, where          //
//                      T(rank*i+j) stores the 'j'th Chebyshev polynomial       //
//                      evaluated at the 'i'th location.                        //
//                                                                              //
/********************************************************************************/
void Chebyshev_polynomials(unsigned rank, double* x, unsigned n, double*& T){
        T       =       new double [n*rank];
        unsigned index;
        if (rank>=1) {
                for (unsigned k=0; k<n; ++k) {
                        T[k*rank]       =       1.0;
                }
        }
	if(rank>=2){
                for (unsigned k=0; k<n; ++k) {
                        T[k*rank+1]     =       x[k];
                }
	}
	if(rank>2){
                for (unsigned k=0; k<n; ++k) {
                        index           =       k*rank;
                        for(int j=2;j<rank;++j){
                                T[index+j] =       2*x[k]*T[index+j-1]-T[index+j-2];
                        }
                }
	}
}

/********************************************************************************/
//      FUNCTION:               get_standard_Chebyshev_nodes                    //
//                                                                              //
//      PURPOSE OF EXISTENCE:   Obtains the standard Chebyshev nodes in [-1,1]. //
//                                                                              //
//      PARAMETERS:                                                             //
//      rank                    -       Number of Chebyshev nodes.              //
//      standardchebnodes       -       The desired set of Chebyshev nodes.     //
//                                                                              //
/********************************************************************************/
void get_standard_Chebyshev_nodes(unsigned rank, double*& standardchebnodes){
        standardchebnodes       =       new double[rank];
        const double PI         =       3.14159265358979323846264338327950;
	for(unsigned j=0; j<rank; ++j){
		standardchebnodes[j]	=	-cos(PI*(2.0*j+1.0)/2.0/rank);
	}
}

/********************************************************************************************************/
//      FUNCTION:               get_standard_Chebyshev_polynomials_evaluated_at_Chebyshev_nodes         //
//                                                                                                      //
//      PURPOSE OF EXISTENCE:   Evaluates the first rank Chebyshev polynomials at the first rank        //
//                              Chebyshev nodes.                                                        //
//                                                                                                      //
//      PARAMETERS:                                                                                     //
//      rank                    -       The number of polynomials and the number of Chebyshev nodes.    //
//      T                       -       Matrix, where T(i,j) stores the 'j'th Chebyshev polynomial      //
//                                      evaluated at the 'i'th location.                                //
//                                                                                                      //
/********************************************************************************************************/
void get_standard_Chebyshev_polynomials_evaluated_at_Chebyshev_nodes(unsigned rank, double*& T){
	double* standardchebnodes;
	get_standard_Chebyshev_nodes(rank, standardchebnodes);
	Chebyshev_polynomials(rank, standardchebnodes, rank, T);
}

/********************************************************************************/
//      FUNCTION:               get_Chebyshev_L2L_Operator                      //
//                                                                              //
//      PURPOSE OF EXISTENCE:   Obtains the interpolation or L2L operator, which//
//                              transfers information from the Chebyshev nodes  //
//                              to the points x.                                //
//                                                                              //
//      PARAMETERS:                                                             //
//      x                       -       Locations of points in interval [-1,1]. //
//      n                       -       Number of points in the interval.       //
//      x_Cheb_Nodes            -       Standard Cheb Nodes in interval [-1,1]. //
//      rank                    -       Number of Chebyshev nodes.              //
//      S                       -       Interpolation or L2L operator.          //
//                                                                              //
/********************************************************************************/
void get_Chebyshev_L2L_Operator(double* x, unsigned n, double* x_Cheb_Nodes, unsigned rank, double*& L2L) {
        double* Tx;
        double* Tcheb;

	Chebyshev_polynomials(rank, x, n, Tx);
	Chebyshev_polynomials(rank, x_Cheb_Nodes, rank, Tcheb);
        
        unsigned index1, index2;

        L2L             =       new double[n*rank];
        double scale    =       1.0/rank;
        for (unsigned i=0; i<n; ++i) {
                index1  =       i*rank;
                for (unsigned j=0; j<rank; ++j) {
                        L2L[index1+j]   =       -1.0;
                        index2          =       j*rank;
                        for (unsigned k=0; k<rank; ++k) {
                                L2L[index1+j]   =       L2L[index1+j]+2.0*Tx[index1+k]*Tcheb[index2+k];
                        }
                        L2L[index1+j]   =       scale*L2L[index1+j];
                }
        }
}


/********************************************************************************/
//      FUNCTION:               scale_Points                                    //
//                                                                              //
//      PURPOSE OF EXISTENCE:   Obtains a set of points in an interval and      //
//                              scales them.                                    //
//                                                                              //
//      PARAMETERS:                                                             //
//      INPUT:                                                                  //
//      x                       -       Intial set of points.                   //
//      center                  -       Intial center of the set of points.     //
//      radius                  -       Intial radius of the set of points.     //
//      N                       -       Number of points in the interval.       //
//      center_New              -       New center of the set of points.        //
//      radius_New              -       New radius of the set of points.        //
//                                                                              //
//      OUTPUT:                                                                 //
//      x_New                   -       New set of points.                      //
//                                                                              //
/********************************************************************************/
void scale_Points(double center, double radius, double* x, unsigned N, double center_New, double radius_New, double*& x_New) {
        x_New   =       new double[N];
        for (unsigned k=0; k<N; ++k) {
                x_New[k]        =       center_New + radius_New*(x[k]-center)/radius;
        }
}