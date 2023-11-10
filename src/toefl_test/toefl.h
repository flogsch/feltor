#pragma once
#include <exception>

#include "dg/algorithm.h"
#include "parameters.h"

namespace toefl
{

template< class Geometry,  class Matrix, class Container >
struct Explicit
{
    Explicit( const Geometry& g, const Parameters& p );

    const Container& phi() const { return m_phi;}
    const Container& var() const { return m_ype;}
    const Container& uE2() const { return m_uE2;}

    dg::Elliptic<Geometry, Matrix, Container>& laplacianM( ) {
        return m_laplaceM;
    }

    unsigned ncalls() const{ return m_ncalls;}

    void operator()( double t, const Container & y,
            Container & yp);

    
  private:
    //use chi and omega as helpers to compute square velocity in uE2
    Container m_chi, m_omega, m_uE2;
    const Container m_binv; //magnetic field

    Container m_phi, m_dxphi, m_dyphi, m_ype;
    Container m_gamma_n;

    //matrices and solvers
    dg::Elliptic<Geometry, Matrix, Container> m_laplaceM;

    std::array<Matrix,2> m_centered;
    dg::Advection< Geometry, Matrix, Container> m_adv;
    dg::ArakawaX< Geometry, Matrix, Container> m_arakawa;
    dg::PCG<Container> m_pcg;
    dg::Extrapolation<Container> m_extra;
    dg::Helmholtz<Geometry, Matrix, Container> m_helmholtz;

    Parameters m_p;

    unsigned m_ncalls = 0;

};

template< class Geometry, class M, class Container>
Explicit< Geometry, M, Container>::Explicit( const Geometry& grid, const Parameters& p ):
    m_chi( evaluate( dg::zero, grid)), m_omega(m_chi), m_uE2(m_chi),
    m_binv( evaluate( dg::LinearX( p.kappa, 1.-p.kappa*p.posX*p.lx), grid)),
    m_phi( m_chi), m_dxphi(m_phi), m_dyphi( m_phi), m_ype(m_phi),
    m_gamma_n(m_chi),
    m_laplaceM( grid,  p.diff_dir),
    m_pcg( m_phi, grid.size()),
    m_extra( 2, m_phi),
    m_helmholtz( -1., {grid, dg::centered}),
    m_adv( grid), m_arakawa(grid),
    m_p(p)   
{
    m_centered = {dg::create::dx( grid, m_p.bcx),
                  dg::create::dy( grid, m_p.bcy)};
}


template< class G, class M, class Container>
void Explicit<G, M, Container>::operator()( double t,
        const Container & y, Container & yp)
{
    m_ncalls ++ ;
    
   
    //need to compute m_phi from y here!!!
    m_extra.extrapolate( t, m_phi);
    m_pcg.solve(m_helmholtz, m_phi, y,
                m_helmholtz.precond(), m_helmholtz.weights(), m_p.eps_gamma[0]);
    m_extra.update( t, m_phi);

    dg::blas1::axpby( 1., m_phi, -1., y, m_chi); //chi = lap \phi
    m_arakawa( m_phi, m_chi, yp);
    //compute derivatives
    dg::blas2::gemv( m_arakawa.dx(), m_phi, m_dxphi);
    dg::blas2::gemv( m_arakawa.dy(), m_phi, m_dyphi);

    //gradient terms
    dg::blas1::axpby( -1, m_dyphi, 1., yp);
    
    dg::blas2::gemv( -m_p.nu, m_laplaceM, y, 1., yp);
    
    
    
    //y[0] = N_e - 1
    //y[1] = N_i - 1 || y[1] = Omega


    ///////////////////////////////////////////////////////////////////////
    


    return;
}

}//namespace dg
