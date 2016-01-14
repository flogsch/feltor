#include <iostream>


#include "blas.h"
#include "dx.h"
#include "sparseblockmat.cuh"
#include "evaluation.cuh"
#include "typedefs.cuh"
#include "weights.cuh"


double function( double x) { return sin(x);}
double derivative( double x) { return cos(x);}
double zero( double x) { return 0;}

double functionX( double x) { 
    if( x < 0) return sin(x);
    else if( 0 <= x && x < 2*M_PI) return cos(x);
    else return sin(x - 2*M_PI);
}
double derivativeX( double x) { 
    if( x < 0) return cos(x);
    else if( 0 <= x && x < 2*M_PI) return -sin(x);
    else return cos(x - 2*M_PI);
}
double zeroX( double x) { return 0;}

typedef dg::HVec Vector;
typedef dg::EllSparseBlockMat Matrix;

int main ()
{
    unsigned n, N;
    std::cout << "Type in n an Nx!\n";
    std::cin >> n>> N;
    std::cout << "# of Legendre nodes " << n <<"\n";
    std::cout << "# of cells          " << N <<"\n";
    dg::Grid1d<double> gPER( 0.1, 2*M_PI+0.1, n, N, dg::PER);
    dg::Grid1d<double> gDIR( 0, M_PI, n, N, dg::DIR);
    dg::Grid1d<double> gNEU( M_PI/2., 3*M_PI/2., n, N, dg::NEU);
    dg::Grid1d<double> gDIR_NEU( 0, M_PI/2., n, N, dg::DIR_NEU);
    dg::Grid1d<double> gNEU_DIR( M_PI/2., M_PI, n, N, dg::NEU_DIR);
    dg::Grid1d<double> g[] = {gPER, gDIR, gNEU, gDIR_NEU,gNEU_DIR};

    std::cout << "TEST NORMAL TOPOLOGY: YOU SHOULD SEE CONVERGENCE FOR ALL OUTPUTS!!!\n";
    for( unsigned i=0; i<5; i++)
    {
        Matrix hs = dg::create::dx( g[i], dg::centered);
        Matrix hf = dg::create::dx( g[i], dg::forward);
        Matrix hb = dg::create::dx( g[i], dg::backward);
        Matrix js = dg::create::jump( g[i].n(), g[i].N(), g[i].h(), g[i].bcx());
        const Vector func = dg::evaluate( function, g[i]);
        Vector error = func;
        const Vector w1d = dg::create::weights( g[i]);
        const Vector deri = dg::evaluate( derivative, g[i]);
        const Vector null = dg::evaluate( zero, g[i]);

        dg::blas2::symv( hs, func, error);
        dg::blas1::axpby( 1., deri, -1., error);
        std::cout << "Distance to true solution (symmetric): "<<sqrt(dg::blas2::dot( w1d, error) )<<"\n";
        dg::blas2::symv( hf, func, error);
        dg::blas1::axpby( 1., deri, -1., error);
        std::cout << "Distance to true solution (forward  ): "<<sqrt(dg::blas2::dot( w1d, error) )<<"\n";
        dg::blas2::symv( hb, func, error);
        dg::blas1::axpby( 1., deri, -1., error);
        std::cout << "Distance to true solution (backward ): "<<sqrt(dg::blas2::dot( w1d, error) )<<"\n";
        dg::blas2::symv( js, func, error);
        dg::blas1::axpby( 1., null , -1., error);
        std::cout << "Distance to true solution (jump     ): "<<sqrt(dg::blas2::dot( w1d, error) )<<"\n\n";
    }
    std::cout << "TEST X-POINT TOPOLOGY: YOU SHOULD SEE CONVERGENCE FOR ALL OUTPUTS!!!\n";
    double fx;
    dg::GridX1d gXDIR( -M_PI, 2*M_PI+M_PI, 1./4., n, N, dg::DIR);
    dg::GridX1d gXDIR0( 0, 2*M_PI, 0., n, N, dg::DIR);
    dg::GridX1d g2[] = {gXDIR, gXDIR0};
    for( unsigned i=0; i<2; i++)
    {
        Matrix hs = dg::create::dx( g2[i], dg::centered);
        Matrix hf = dg::create::dx( g2[i], dg::forward);
        Matrix hb = dg::create::dx( g2[i], dg::backward);
        Matrix js = dg::create::jump( g2[i]);
        const Vector func = dg::evaluate( functionX, g2[i]);
        Vector error = func;
        const Vector w1d = dg::create::weights( g2[i]);
        const Vector deri = dg::evaluate( derivativeX, g2[i]);
        const Vector null = dg::evaluate( zeroX, g2[i]);

        dg::blas2::symv( hs, func, error);
        dg::blas1::axpby( 1., deri, -1., error);
        std::cout << "Distance to true solution (symmetric): "<<sqrt(dg::blas2::dot( w1d, error) )<<"\n";
        dg::blas2::symv( hf, func, error);
        dg::blas1::axpby( 1., deri, -1., error);
        std::cout << "Distance to true solution (forward  ): "<<sqrt(dg::blas2::dot( w1d, error) )<<"\n";
        dg::blas2::symv( hb, func, error);
        dg::blas1::axpby( 1., deri, -1., error);
        std::cout << "Distance to true solution (backward ): "<<sqrt(dg::blas2::dot( w1d, error) )<<"\n";
        dg::blas2::symv( js, func, error);
        dg::blas1::axpby( 1., null , -1., error);
        std::cout << "Distance to true solution (jump     ): "<<sqrt(dg::blas2::dot( w1d, error) )<<"\n\n";
    }
    //for periodic bc | dirichlet bc
    //n = 1 -> p = 2      2
    //n = 2 -> p = 1      1
    //n = 3 -> p = 3      3
    //n = 4 -> p = 3      3
    //n = 5 -> p = 5      5


    
    return 0;
}
