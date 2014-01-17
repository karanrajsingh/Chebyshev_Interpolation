//
//  Chebyshev_Interpolation_2D.hpp
//  
//
//  Created by Sivaram Ambikasaran on 1/10/14.
//
//

#ifndef __CHEBYSHEV_INTERPOLATION_2D__
#define __CHEBYSHEV_INTERPOLATION_2D__

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
void function2D(double*& x, double*& y, unsigned n, double*& f);


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
void kernel2D(double*& x1, double*& y1, unsigned n1, double*& x2, double*& y2, unsigned n2, double*& K);

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
//      INPUTS:                                                                 //
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
void get_Scaled_Chebyshev_Nodes(double x_Center, double x_Radius, double y_Center, double y_Radius, unsigned rank, double*& Cheb_Node, double*& x_Cheb_Node, double*& y_Cheb_Node);

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
void get_Chebyshev_L2L_Operator(double*& x, double*& y, unsigned n, double*& Cheb_Node, unsigned rank, double*& L2L);

#endif /* defined(__CHEBYSHEV_INTERPOLATION_2D__) */