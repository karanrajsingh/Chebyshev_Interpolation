//
//  Test_Chebyshev_1D.cpp
//  
//
//  Created by Sivaram Ambikasaran on 1/17/14.
//
//

#include <iostream>
#include <cstdlib>
#include "Chebyshev_Interpolation_1D.hpp"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

void get_Points(double center, double radius, unsigned N, double*& x) {
        x               =       new double[N];
        double RAND     =       RAND_MAX;
        for (unsigned k=0; k<N; ++k) {
                x[k]    =       center + radius*(2*double(rand())/RAND-1);
        }
}

int main() {
        srand(time(NULL));

        //      Obtain the points in both the clusters.
        unsigned n1     =       5000;
        unsigned n2     =       5000;

        double* x1;
        double* x2;
        double center1  =       -1;
        double center2  =       +1;
        double radius1  =       0.5;
        double radius2  =       0.5;

        get_Points(center1, radius1, n1, x1);
        get_Points(center2, radius2, n2, x2);

        //      Obtain the exact interaction between the two clusters.
        double* Kexact;
        kernel1D(x1, n1, x2, n2, Kexact);

        //      Obtain the standard Chebyshev nodes.
        double* Cheb_Nodes;
        unsigned rank   =       20;
        get_standard_Chebyshev_nodes(rank, Cheb_Nodes);

        //      Obtain the scaled Chebyshev nodes for the first cluster.
        double* x1_Cheb_Nodes;
        scale_Points(0, 1, Cheb_Nodes, rank, center1, radius1, x1_Cheb_Nodes);

        //      Obtain the scaled Chebyshev nodes for the second cluster.
        double* x2_Cheb_Nodes;
        scale_Points(0, 1, Cheb_Nodes, rank, center2, radius2, x2_Cheb_Nodes);

        //      Obtain the standard location for points in the first cluster.
        double* x1_Standard_Location;
        scale_Points(center1, radius1, x1, n1, 0, 1, x1_Standard_Location);

        //      Obtain the standard location for points in the second cluster.
        double* x2_Standard_Location;
        scale_Points(center2, radius2, x2, n2, 0, 1, x2_Standard_Location);

        //      Obtain L2L for first cluster
        double* L2L1;
        get_Chebyshev_L2L_Operator(x1_Standard_Location, n1, Cheb_Nodes, rank, L2L1);

        //      Obtain L2L for second cluster
        double* L2L2;
        get_Chebyshev_L2L_Operator(x2_Standard_Location, n2, Cheb_Nodes, rank, L2L2);

        //      Obtain M2L
        double* M2L;
        kernel1D(x1_Cheb_Nodes, rank, x2_Cheb_Nodes, rank, M2L);

        Map<Matrix<double,Dynamic,Dynamic,RowMajor> >  Kexact_E(Kexact, n1, n2);
        Map<Matrix<double,Dynamic,Dynamic,RowMajor> >  L2L1_E(L2L1, n1, rank);
        Map<Matrix<double,Dynamic,Dynamic,RowMajor> >  L2L2_E(L2L2, n2, rank);
        Map<Matrix<double,Dynamic,Dynamic,RowMajor> >  M2L_E(M2L, rank, rank);

        cout << endl << "Number of points in the first cluster centered at " << center1 << " of length " << 2*radius1 << " is: " << n1 << endl;
        cout << endl << "Number of points in the second cluster centered at " << center2 << " of length " << 2*radius2 << " is: " << n2 << endl;
        cout << endl << "Rank of interaction considered is: " << rank << endl;
        cout << endl << "Maximum error in the low-rank interaction between the two cluster is: " << (Kexact_E-L2L1_E*M2L_E*L2L2_E.transpose()).cwiseAbs().maxCoeff() << endl;
}