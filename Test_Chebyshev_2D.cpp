//
//  Test_Chebyshev_2D.cpp
//  
//
//  Created by Sivaram Ambikasaran on 1/17/14.
//
//

#include <iostream>
#include <cstdlib>
#include "Chebyshev_Interpolation_2D.hpp"
#include "Chebyshev_Interpolation_1D.hpp"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

void get_Points_In_Standard_Square(unsigned N, double*& x, double*& y) {
        x               =       new double[N];
        y               =       new double[N];
        double RAND     =       RAND_MAX;
        for (unsigned k=0; k<N; ++k) {
                x[k]    =       2*double(rand())/RAND-1;
                y[k]    =       2*double(rand())/RAND-1;
        }
}

int main() {
        srand(time(NULL));

        //      Obtain the points in first cluster.
        unsigned n1     =       10000;
        double* x1_Standard_Location;
        double* y1_Standard_Location;
        get_Points_In_Standard_Square(n1, x1_Standard_Location, y1_Standard_Location);

        double xcenter1 =       -1;
        double xradius1 =       0.5;

        double ycenter1 =       -1;
        double yradius1 =       0.5;
        
        double* x1;
        double* y1;
        scale_Points(0, 1, x1_Standard_Location, n1, xcenter1, xradius1, x1);
        scale_Points(0, 1, y1_Standard_Location, n1, ycenter1, yradius1, y1);
        
        //      Obtain the points in second cluster.
        unsigned n2     =       10000;
        double* x2_Standard_Location;
        double* y2_Standard_Location;
        get_Points_In_Standard_Square(n2, x2_Standard_Location, y2_Standard_Location);
        
        double xcenter2 =       1;
        double xradius2 =       0.5;
        
        double ycenter2 =       1;
        double yradius2 =       0.5;

        double* x2;
        double* y2;
        scale_Points(0, 1, x2_Standard_Location, n2, xcenter2, xradius2, x2);
        scale_Points(0, 1, y2_Standard_Location, n2, ycenter2, yradius2, y2);


        //      Obtain the exact interaction between the two clusters.
        double* Kexact;
        kernel2D(x1, y1, n1, x2, y2, n2, Kexact);

        //      Obtain the standard Chebyshev nodes in the standard interval.
        double* Cheb_Nodes;
        unsigned rank   =       7;
        get_standard_Chebyshev_nodes(rank, Cheb_Nodes);

        //      Obtain the standard Chebyshev nodes in the standard square.
        double* x_Cheb_Nodes;
        double* y_Cheb_Nodes;
        get_Scaled_Chebyshev_Nodes(0, 1, 0, 1, rank, Cheb_Nodes, x_Cheb_Nodes, y_Cheb_Nodes);
        
        //      Obtain the scaled Chebyshev nodes for the first cluster.
        unsigned RANK   =       rank*rank;
        double* x1_Cheb_Nodes;
        double* y1_Cheb_Nodes;
        scale_Points(0, 1, x_Cheb_Nodes, RANK, xcenter1, xradius1, x1_Cheb_Nodes);
        scale_Points(0, 1, y_Cheb_Nodes, RANK, ycenter1, yradius1, y1_Cheb_Nodes);

        //      Obtain the scaled Chebyshev nodes for the second cluster.
        double* x2_Cheb_Nodes;
        double* y2_Cheb_Nodes;
        scale_Points(0, 1, x_Cheb_Nodes, RANK, xcenter2, xradius2, x2_Cheb_Nodes);
        scale_Points(0, 1, y_Cheb_Nodes, RANK, ycenter2, yradius2, y2_Cheb_Nodes);

        //      Obtain L2L for first cluster
        double* L2L1;
        get_Chebyshev_L2L_Operator(x1_Standard_Location, y1_Standard_Location, n1, Cheb_Nodes, rank, L2L1);

        //      Obtain L2L for second cluster
        double* L2L2;
        get_Chebyshev_L2L_Operator(x2_Standard_Location, y2_Standard_Location, n2, Cheb_Nodes, rank, L2L2);

        //      Obtain M2L
        double* M2L;
        kernel2D(x1_Cheb_Nodes, y1_Cheb_Nodes, RANK, x2_Cheb_Nodes, y2_Cheb_Nodes, RANK, M2L);

        Map<Matrix<double,Dynamic,Dynamic,RowMajor> >  Kexact_E(Kexact, n1, n2);
        Map<Matrix<double,Dynamic,Dynamic,RowMajor> >  L2L1_E(L2L1, n1, RANK);
        Map<Matrix<double,Dynamic,Dynamic,RowMajor> >  L2L2_E(L2L2, n2, RANK);
        Map<Matrix<double,Dynamic,Dynamic,RowMajor> >  M2L_E(M2L, RANK, RANK);

        cout << endl << "Number of points in the first cluster centered at (" << xcenter1 << ", " << ycenter1 << ") with side of length " << 2*xradius1 << " is: " << n1 << endl;
        cout << endl << "Number of points in the second cluster centered at (" << xcenter2 <<  ", " << ycenter2 << ") with side of length " << 2*xradius2 << " is: " << n2 << endl;
        cout << endl << "Rank of interaction considered is: " << RANK << endl;
        cout << endl << "Maximum error in the low-rank interaction between the two cluster is: " << (Kexact_E-L2L1_E*M2L_E*L2L2_E.transpose()).cwiseAbs().maxCoeff() << endl;
}