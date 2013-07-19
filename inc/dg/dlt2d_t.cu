#include <iostream>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <cusp/ell_matrix.h>
#include <cusp/dia_matrix.h>

#include "blas.h"
#include "laplace.cuh"
#include "array.cuh"
#include "dlt.h"
#include "arrvec1d.cuh"
#include "evaluation.cuh"
#include "operator.cuh"
#include "operator_matrix.cuh"
#include "tensor.cuh"
#include "typedefs.cuh"

using namespace std;
using namespace dg;

const unsigned n = 3; 
const unsigned Nx = 2e0;
const unsigned Ny = 2e0;

double function( double x, double y ) { return sin(x)*sin(y);}

int main()
{
    cout << "# of Legendre coefficients: " << n<<endl;
    cout << "# of grid cells:            " << Nx*Ny<<endl;
    Grid<double,n> grid( 0, 2.*M_PI, 0, 2.*M_PI, Nx, Ny);
    ArrVec2d<double, n, HVec> hv ( evaluate( function, grid), grid.Nx());
    DVec dv( hv.data());
    DVec dv2( dv), dv3( dv);
    //cout << "Evaluated\n";
    //cout <<hv<<endl;

    Operator<double, n> forward( DLT<n>::forward);
    Operator<double, n*n> forward2d = tensor( forward, forward);
    dg::blas2::symv( 1., forward2d, dv,0., dv);
    dg::blas2::symv( forward2d, dv2, dv2);

    //HMatrix hm = tensor<double, n>( tensor(Ny, forward), tensor(Nx, forward));
    HMatrix hm = tensor( Nx*Ny, forward2d);
    DMatrix dm(hm);
    dg::blas2::symv( dm, dv3, dv3);

    ArrVec2d<double,n, HVec> hv2 ( expand( function, grid), grid.Nx()); 

    //test for equality...
    hv.data() = dv;
    cout << "Multiplied 1st version\n";
    cout << hv <<endl;
    cout << "Multiplied 2nd version\n";
    hv.data() = dv2; 
    cout << hv<<endl;
    cout << "Multiplied 3rd version\n";
    hv.data() = dv3; 
    cout << hv<<endl;
    cout << "Correct is: "<<endl;
    cout << hv2 <<endl;
    
    return 0;
}
