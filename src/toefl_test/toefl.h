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

    const Container& phi(unsigned i ) const { return m_phi[i];}
    const Container& var(unsigned i ) const { return m_ype[i];}
    const Container& uE2() const { return m_uE2;}

    dg::Elliptic<Geometry, Matrix, Container>& laplacianM( ) {
        return m_laplaceM;
    }

    dg::Helmholtz<Geometry, Matrix, Container >&  gamma_inv() {
        return m_multi_gamma1[0];
    }
    unsigned ncalls() const{ return m_ncalls;}

    void operator()( double t, const std::array<Container,2>& y,
            std::array<Container,2>& yp);

    
  private:
    //use chi and omega as helpers to compute square velocity in uE2
    Container m_chi, m_omega, m_uE2;
    const Container m_binv; //magnetic field

    std::array<Container,2> m_phi, m_dxphi, m_dyphi, m_ype;
    std::array<Container,2> m_lapy, m_v;
    Container m_gamma_n;

    //matrices and solvers
    dg::Elliptic<Geometry, Matrix, Container> m_laplaceM;
    std::vector<dg::Elliptic<Geometry, Matrix, Container> > m_multi_pol;
    std::vector<dg::Helmholtz<Geometry,  Matrix, Container> > m_multi_gamma1;
    std::array<Matrix,2> m_centered;
    dg::Advection< Geometry, Matrix, Container> m_adv;
    dg::ArakawaX< Geometry, Matrix, Container> m_arakawa;
    dg::PCG<Container> m_pcg;
    dg::Extrapolation<Container> m_extra;
    dg::Helmholtz<Geometry, Matrix, Container> m_helmholtz;

    dg::MultigridCG2d<Geometry, Matrix, Container> m_multigrid;
    dg::Extrapolation<Container> m_old_phi, m_old_psi, m_old_gammaN;
    std::vector<Container> m_multi_chi;

    Parameters m_p;

    unsigned m_ncalls = 0;

};

template< class Geometry, class M, class Container>
Explicit< Geometry, M, Container>::Explicit( const Geometry& grid, const Parameters& p ):
    m_chi( evaluate( dg::zero, grid)), m_omega(m_chi), m_uE2(m_chi),
    m_binv( evaluate( dg::LinearX( p.kappa, 1.-p.kappa*p.posX*p.lx), grid)),
    m_phi( {m_chi, m_chi}), m_dxphi(m_phi), m_dyphi( m_phi), m_ype(m_phi),
    m_lapy(m_phi), m_v(m_phi),
    m_gamma_n(m_chi),
    m_laplaceM( grid,  p.diff_dir),
    m_pcg( m_phi[0], grid.size()),
    m_extra( 2, m_phi[0]),
    m_helmholtz( -1., {grid, dg::centered}),
    m_adv( grid), m_arakawa(grid),
    m_multigrid( grid, p.num_stages),
    m_old_phi( 2, m_chi), m_old_psi( 2, m_chi), m_old_gammaN( 2, m_chi),
    m_p(p)
    
{
    m_multi_chi= m_multigrid.project( m_chi);
    for( unsigned u=0; u<p.num_stages; u++)
    {
        m_multi_pol.push_back({ m_multigrid.grid(u),  p.pol_dir, 1.});
        m_multi_gamma1.push_back({-0.5*p.tau, { m_multigrid.grid(u), p.pol_dir}});
    }
    m_centered = {dg::create::dx( grid, m_p.bcx),
                  dg::create::dy( grid, m_p.bcy)};
}


template< class G, class M, class Container>
void Explicit<G, M, Container>::operator()( double t,
        const std::array<Container,2>& y, std::array<Container,2>& yp)
{
    m_ncalls ++ ;
    
   
    //need to compute m_phi from y here!!!
    m_extra.extrapolate( t, m_phi[0]);
    m_pcg.solve(m_helmholtz, m_phi[0], y[0],
                m_helmholtz.precond(), m_helmholtz.weights(), m_p.eps_gamma[0]);
    m_extra.update( t, m_phi[0]);

    dg::blas1::axpby( 1., m_phi[0], -1., y[0], m_chi); //chi = lap \phi
    m_arakawa( m_phi[0], m_chi, yp[0]);
    //compute derivatives
    dg::blas2::gemv( m_arakawa.dx(), m_phi[0], m_dxphi[0]);
    dg::blas2::gemv( m_arakawa.dy(), m_phi[0], m_dyphi[0]);

    //gradient terms
    dg::blas1::axpby( -1, m_dyphi[0], 1., yp[0]);
    
    dg::blas2::gemv( -m_p.nu, m_laplaceM, y[0], 1., yp[0]);
    
    
    
    //y[0] = N_e - 1
    //y[1] = N_i - 1 || y[1] = Omega


    ///////////////////////////////////////////////////////////////////////
    


    return;
}

}//namespace dg
