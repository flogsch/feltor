#include <iostream>
#include <cmath>
#include <cusp/print.h>
#include "dg/algorithm.h"
#include "fem.h"
#include "fem_weights.h"

double function( double x, double y){return sin(x)*cos(y);}
double function( double x, double y, double z){return sin(x)*cos(y)*cos(z);}

using Vector = dg::DVec;
using Matrix = cusp::coo_matrix<int,double,cusp::device_memory>;
using MassMatrix = dg::KroneckerTriDiagonal2d<Vector>;
using InvMassMatrix = dg::InverseKroneckerTriDiagonal2d<Vector>;

int main ()
{
    unsigned n = 3, Nx = 18, Ny = 24, mx = 3;
    double eps = 1e-10;
    //std::cout << "# Type in n Nx Ny mx eps!\n";
    std::cout << "# on grid " << n <<" x "<<Nx<<" x "<<Ny<<"\n";
    std::cout << "# eps and Multiply " << eps <<" " << mx<<"\n";
    dg::CartesianGrid2d gDIR( 0, 2.*M_PI, M_PI/2., 5*M_PI/2., n, Nx, Ny, dg::DIR,
            dg::DIR);
    dg::FemRefinement fem_ref( mx);
    dg::CartesianRefinedGrid2d gDIR_f( fem_ref, fem_ref, gDIR.x0(), gDIR.x1(),
            gDIR.y0(), gDIR.y1(), n, Nx,Ny, dg::DIR, dg::DIR);
    const Vector func = dg::evaluate( function, gDIR);
    const Vector v2d = dg::create::fem_inv_weights( gDIR);
    const Vector w2d = dg::create::fem_weights( gDIR);
    const Vector wf2d = dg::create::volume( gDIR_f);
    double integral = dg::blas2::dot( func, w2d, func);
    std::cout << "error of integral is "
              <<(integral-M_PI*M_PI)/M_PI/M_PI<<std::endl;
    dg::HVec Xf = dg::pullback( dg::cooX2d, gDIR_f);
    dg::HVec Yf = dg::pullback( dg::cooY2d, gDIR_f);
    Matrix inter = dg::create::interpolation( Xf, Yf, gDIR, dg::NEU, dg::NEU, "linear");
    Matrix interT = dg::transpose( inter);
    Matrix Wf = dg::create::diagonal( (dg::HVec)wf2d), project;
    Matrix Vf = dg::create::diagonal( (dg::HVec)v2d), tmp;
    cusp::multiply( interT, Wf, tmp);
    cusp::multiply( Vf, tmp, project);
    project.sort_by_row_and_column();
    //std::cout << "Project matrix\n";
    //cusp::print( project);


    Vector func_f( gDIR_f.size());
    dg::blas2::symv( inter, func, func_f);
    integral = dg::blas2::dot( func_f, wf2d, func_f);
    std::cout << "error of refined integral is "
              <<(integral-M_PI*M_PI)/M_PI/M_PI<<std::endl;
    Vector barfunc(func);
    dg::blas2::symv( project, func_f, barfunc);
    // test now should contain Sf
    Vector test( barfunc);
    dg::PCG<Vector> cg( test, 1000);
    // PCG tests fem-mass
    MassMatrix fem_mass = dg::create::fem_mass( gDIR);
    //std::cout << "S matrix\n";
    //cusp::print( fem_mass);
    unsigned number = cg.solve( fem_mass, test, barfunc, 1., w2d, eps);
    dg::blas1::axpby( 1., func, -1., test);
    double norm = sqrt(dg::blas2::dot( w2d, test) );
    double func_norm = sqrt(dg::blas2::dot( w2d, func) );
    std::cout <<"PCG Distance to true solution: "<<norm/func_norm<<"\n";
    std::cout << "using "<<number<<" iterations\n";
    InvMassMatrix inv_fem_mass = dg::create::inv_fem_mass( gDIR);
    dg::blas2::symv( inv_fem_mass, barfunc, test);
    dg::blas1::axpby( 1., func, -1., test);
    norm = sqrt(dg::blas2::dot( w2d, test) );
    std::cout <<"Thomas Distance to true solution: "<<norm/func_norm<<"\n";


    std::cout << "TEST L2C projection\n";
    Matrix interC = dg::create::interpolation( Xf, Yf, gDIR, dg::NEU, dg::NEU, "nearest");
    interT = dg::transpose( interC);
    cusp::multiply( interT, Wf, tmp);
    cusp::multiply( Vf, tmp, project);
    project.sort_by_row_and_column();
    dg::blas2::symv( project, func_f, barfunc);
    fem_mass = dg::create::fem_linear2const( gDIR);

    number = cg.solve( fem_mass, test, barfunc, 1., w2d, eps);
    dg::blas1::axpby( 1., func, -1., test);
    norm = sqrt(dg::blas2::dot( w2d, test) );
    std::cout <<"PCG Distance to true solution: "<<norm/func_norm<<"\n";
    std::cout << "using "<<number<<" iterations\n";
    inv_fem_mass = dg::create::inv_fem_linear2const( gDIR);
    dg::blas2::symv( inv_fem_mass, barfunc, test);
    dg::blas1::axpby( 1., func, -1., test);
    norm = sqrt(dg::blas2::dot( w2d, test) );
    std::cout <<"Thomas Distance to true solution: "<<norm/func_norm<<"\n";
    std::cout << "TEST 3d grid\n";
    dg::CartesianGrid3d g3d( 0, 2.*M_PI, M_PI/2., 5*M_PI/2., 0,1, n, Nx, Ny,10,
            dg::DIR, dg::DIR, dg::DIR);
    dg::DVec x = dg::evaluate( function,  g3d), y(x), z(x);
    const dg::DVec w3d = dg::create::fem_weights( g3d);
    auto split_x = dg::split(x, g3d);
    auto split_y = dg::split(y, g3d);
    auto split_z = dg::split(z, g3d);
    MassMatrix fem_mass2d = dg::create::fem_linear2const2d( g3d);
    dg::blas2::symv( fem_mass2d, x, y);
    for( unsigned i=0; i<g3d.Nz(); i++)
        dg::blas2::symv( fem_mass, split_x[i], split_z[i]);
    dg::blas1::axpby( 1., y, -1., z);
    double err = dg::blas2::dot( z, w3d, z);
    std::cout << "1. Error in 3d is "<<err<<"\n";
    InvMassMatrix inv_fem_mass2d = dg::create::inv_fem_linear2const2d( g3d);
    dg::blas2::symv( inv_fem_mass2d, y, z);
    for( unsigned i=0; i<g3d.Nz(); i++)
        dg::blas2::symv( inv_fem_mass, split_y[i], split_x[i]);
    dg::blas1::axpby( 1., z, -1., x);
    err = dg::blas2::dot( x, w3d, x);
    std::cout << "2. Error in 3d is "<<err<<"\n";
    x = dg::evaluate( function,  g3d);
    dg::blas1::axpby( 1., x, -1., z);
    err = dg::blas2::dot( z, w3d, z);
    std::cout << "Abs Error in 3d is "<<err<<"\n";


    return 0;
}
