//
//  Chebyshev_Interpolation_2D.cpp
//  
//
//  Created by Sivaram Ambikasaran on 1/10/14.
//
//

#include <cmath>
#include "Chebyshev_Interpolation_1D.hpp"
#include "Chebyshev_Interpolation_2D.hpp"

/********************************************************************************/
//      FUNCTION:               function2D                                      //
//                                                                              //
//      PURPOSE OF EXISTENCE:   Computes a function in 2D.                      //
//                                                                              //
//      PARAMETERS:                                                             //
//      x       -       'x' location of the points in the cluster.              //
//      y       -       'y' location of the points in the cluster.              //
//      n       -       Number of points in the cluster.                        //
//      f       -       Vector with the function values evalauted at the input  //
//                      'x' and 'y' locations.                                  //
//                                                                              //
/********************************************************************************/
void function2D(double*& x, double*& y, unsigned n, double*& f) {
        f       =       new double[n];
        for (unsigned j=0; j<n; ++j) {
                f[j]    =       exp(-x[j]*y[j]-2.0*y[j]);
        }
}


/********************************************************************************/
//      FUNCTION:               kernel2D                                        //
//                                                                              //
//      PURPOSE OF EXISTENCE:   Computes the kernel in 2D.                      //
//                                                                              //
//      PARAMETERS:                                                             //
//      x1      -       'x' location of the points in the first cluster.        //
//      y1      -       'y' location of the points in the first cluster.        //
//      n1      -       Number of points in the first cluster.                  //
//      x2      -       'x' location of the points in the second cluster.       //
//      y2      -       'y' location of the points in the second cluster.       //
//      n2      -       Number of points in the second cluster.                 //
//      K       -       Matrix, where K(i,j) is the interaction between the     //
//                      'i'th point in the first cluster and 'j'th point in     //
//                      the second cluster.                                     //
//                                                                              //
/********************************************************************************/
void kernel2D(double*& x1, double*& y1, unsigned n1, double*& x2, double*& y2, unsigned n2, double*& K) {
        double Rsquare;
        unsigned index;
        K       =       new double[n1*n2];
        for (unsigned j=0; j<n1; ++j) {
                index   =       j*n2;
                for (unsigned k=0; k<n2; ++k) {
                        Rsquare         =       (x1[j]-x2[k])*(x1[j]-x2[k])+(y1[j]-y2[k])*(y1[j]-y2[k]);
                        K[index+k]      =       0.5*log(Rsquare);
                }
        }
}

/********************************************************************************/
//      FUNCTION:               get_Scaled_Chebyshev_Nodes                      //
//                                                                              //
//      PURPOSE OF EXISTENCE:   Obtains the scaled Chebyshev nodes in 2D ordered//
//                              as shown below for rank = 3.                    //
//                                                                              //
//               _________________                                              //
//              |     |     |     |                                             //
//              |  7  |  8  |  9  |                                             //
//              |_____|_____|_____|                                             //
//              |     |     |     |                                             //
//              |  4  |  5  |  6  |                                             //
//              |_____|_____|_____|                                             //
//              |     |     |     |                                             //
//              |  1  |  2  |  3  |                                             //
//              |_____|_____|_____|                                             //
//                                                                              //
//      PARAMETERS:                                                             //
//                                                                              //
//      x_Center        -       The 'x' coordinate of the center of cluster.    //
//      y_Center        -       The 'y' coordinate of the center of cluster.    //
//      x_Radius        -       Radius of the cluster along the X direction.    //
//      y_Radius        -       Radius of the cluster along the Y direction.    //
//      rank            -       Number of Chebyshev nodes along one direction.  //
//      Cheb_Node       -       Chebyshev nodes in [-1,1].                      //
//      x_Cheb_Node     -       Scaled 'x' locations of the 'rank*rank'         //
//                              Chebyshev nodes in the square.                  //
//      y_Cheb_Node     -       Scaled 'y' locations of the 'rank*rank'         //
//                              Chebyshev nodes in the square.                  //
/********************************************************************************/
void get_Scaled_Chebyshev_Nodes(double x_Center, double x_Radius, double y_Center, double y_Radius, unsigned rank, double*& Cheb_Node, double*& x_Cheb_Node, double*& y_Cheb_Node) {

        unsigned RANK   =       rank*rank;

        x_Cheb_Node     =       new double[RANK];
        y_Cheb_Node     =       new double[RANK];

        unsigned jx, jy;

        for (unsigned j=0; j<RANK; ++j) {
                jx              =       j%rank;
                jy              =       j/rank;

                x_Cheb_Node[j]  =       x_Center        +       x_Radius*Cheb_Node[jx];
                y_Cheb_Node[j]  =       y_Center        +       y_Radius*Cheb_Node[jy];
        }
}

/********************************************************************************/
//      FUNCTION:               get_Chebyshev_L2L_Operator                      //
//                                                                              //
//      PURPOSE OF EXISTENCE:   Obtains the Chebyshev L2L Operator over the     //
//                              square [-1,1]^2.                                //
//                                                                              //
//      PARAMETERS:                                                             //
//      x               -       'x' location of points in the square [-1,1]^2.  //
//      y               -       'y' location of points in the square [-1,1]^2.  //
//      n               -       Total number of points.                         //
//      Cheb_Node       -       Location of Chebyshev nodes in [-1,1].          //
//      rank            -       Number of Chebyshev nodes in [-1,1].            //
//      L2L             -       Interpolation or L2L operator, which transfers  //
//                              information from parent to child.               //
/********************************************************************************/
void get_Chebyshev_L2L_Operator(double*& x, double*& y, unsigned n, double*& Cheb_Node, unsigned rank, double*& L2L) {

        double* L2Lx;
        get_Chebyshev_L2L_Operator(x, n, Cheb_Node, rank, L2Lx);

        double* L2Ly;
        get_Chebyshev_L2L_Operator(y, n, Cheb_Node, rank, L2Ly);

        unsigned RANK   =       rank*rank;

        unsigned jx, jy;

        L2L     =       new double[n*RANK];
        unsigned index1, index2;

        for (unsigned i=0; i<n; ++i) {
                index1  =       i*RANK;
                index2  =       i*rank;
                for (unsigned j=0; j<RANK; ++j) {
                        jx              =       j%rank;
                        jy              =       j/rank;

                        L2L[index1+j]   =       L2Lx[index2+jx]*L2Ly[index2+jy];
                }
        }
}