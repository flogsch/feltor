#include <iostream>

#include "blas.h"

#include "helmholtz.h"
#include "backend/typedefs.h"
#include "multistep.h"

#include "pcg.h"
#include "dg/algorithm.h"

const double eps = 1e-3;
const double alpha = 1.;
double lhs( double x, double y){ return 0.25*sin(x)*sin(y);}
double rhs( double x, double y){ return  exp(0.25*sin(x)*sin(y))*(1.+ 0.5*sin(x)*sin(y)) ;}
const double R_0 = 1000;
const double lx = 2.*M_PI;
const double ly = 2.*M_PI;
const double lz = 2.*M_PI;
double fct(double x, double y, double z){ return sin(x-R_0)*sin(y);}
double laplace_fct( double x, double y, double z) {
    return -1./x*cos(x-R_0)*sin(y) + 2.*sin(x-R_0)*sin(y);}
double helmholtz_fct( double x, double y, double z) {
    return fct(x,y,z) - alpha*laplace_fct(x,y,z);}

dg::bc bcx = dg::DIR;
double initial( double x, double y, double z) {return sin(0);}
int main()
{

    unsigned n=2, Nx=99, Ny=99;//, Nz=1;
    //std::cout << "Type n, Nx Ny and Nz\n";
    //std::cin >> n>> Nx >> Ny >> Nz;
    dg::Grid2d grid( 0, 2.*M_PI, 0, 2.*M_PI, n, Nx, Ny, dg::DIR, dg::PER);
    const dg::DVec w2d = dg::create::weights( grid);
    const dg::DVec rho = dg::evaluate( rhs, grid);
    const dg::DVec sol = dg::evaluate( lhs, grid);
    dg::DVec x(rho.size(), 0.);

    dg::HelmholtzLN<dg::CartesianGrid2d, dg::DMatrix, dg::DVec > helmln( alpha, {grid});

    std::cout << "FIRST METHOD:\n";
    dg::PCG< dg::DVec > pcg(x, x.size());
    unsigned number = 1;//pcg.solve( helmln, x, rho, 1., w2d, eps);

    std::cout << "SECOND METHOD:\n";
    dg::ChebyshevIteration< dg::DVec > cheby(x);
    dg::DVec x_(rho.size(), 0.);
    //cheby.solve( helmln, x_, rho, 0., 1., 50, true);

    std::cout << "THIRD METHOD:\n";
    dg::LGMRES<dg::DVec> lgmres(x, 30, 3, 10);
    dg::DVec x__(rho.size(), 0.);
    //lgmres.solve(helmln, x__, rho, 1., w2d, eps, 1.);

    std::cout << "FOURTH METHOD:\n";
    dg::AndersonAcceleration<dg::DVec> anderson(x,0);
    dg::DVec x___(rho.size(), 10.);
    anderson.solve(helmln, x___, rho,  w2d, eps, eps, 100, 1e-3, 10, false);


    for (double i: x___ ){
        std::cout << i << " ";
    }
    //Evaluation
    dg::blas1::axpby( 1., sol, -1., x);
    dg::blas1::axpby( 1., sol, -1., x_);
    dg::blas1::axpby( 1., sol, -1., x__);

    std::cout << "number of iterations:  "<<number<<std::endl;
    std::cout << "ALL METHODS SHOULD DO THE SAME!\n";
    dg::exblas::udouble res;
    res.d = sqrt( dg::blas2::dot( w2d, x));
    std::cout << "error1 " << res.d<<"\t"<<res.i<<std::endl;
    res.d = sqrt( dg::blas2::dot( w2d, x_));
    std::cout << "error2 " << res.d<<"\t"<<res.i<<std::endl;
    res.d = sqrt( dg::blas2::dot( w2d, x__));
    std::cout << "error3 " << res.d<<"\t"<<res.i<<std::endl;

    /*
    std::cout << "Test 3d cylincdrical norm:\n";
    dg::CylindricalGrid3d g3d( R_0, R_0+lx, 0, ly, 0,lz, n, Nx, Ny,Nz, bcx, dg::PER, dg::PER);
    dg::DVec fct_ = dg::evaluate(fct, g3d );
    dg::DVec laplace_fct_ = dg::evaluate( laplace_fct, g3d);
    dg::DVec helmholtz_fct_ = dg::evaluate( helmholtz_fct, g3d);
    dg::DVec temp_(fct_);
    dg::Elliptic< dg::CylindricalGrid3d, dg::DMatrix, dg::DVec > laplaceM( g3d);
    dg::Helmholtz< dg::CylindricalGrid3d, dg::DMatrix, dg::DVec > helmholtz( alpha, {g3d});
    dg::blas2::symv( laplaceM, fct_, temp_);
    dg::blas1::axpby( 1., laplace_fct_, -1., temp_);
    dg::DVec w3d =  laplaceM.weights();
    std::cout << "error Laplace " << sqrt( dg::blas2::dot( w3d, temp_))<<" (Note the supraconvergence!)"<<std::endl;
    dg::blas2::symv( helmholtz, fct_, temp_);
    dg::blas1::axpby( 1., helmholtz_fct_, -1, temp_);
    std::cout << "error " << sqrt( dg::blas2::dot( w3d, temp_))<<" (Note the supraconvergence!)"<<std::endl;
    */
  




    return 0;
}
