#pragma once
#include <exception>

#include "dg/algorithm.h"
#include "parameters.h"

namespace mima
{

template< class Geometry,  class Matrix, class Container >
struct Explicit
{
    Explicit( const Geometry& g, const Parameters& p );

    const Container& chi() const { return m_chi;}
    const Container& phi() const { return m_phi;}
    const Container& uE2() const { return m_uE2;}        

    dg::Elliptic<Geometry, Matrix, Container>& laplacianM( ) {
        return m_laplaceM;
    }

    unsigned ncalls() const{ return m_ncalls;}

    void operator()( double t, const Container & y,
            Container & yp);

    
  private:
    //use chi and omega as helpers to compute square velocity in uE2
    Container m_chi, m_uE2;
    Container m_phi, m_dxphi, m_dyphi, m_vx, m_vy;

    //matrices and solvers
    dg::Elliptic<Geometry, Matrix, Container> m_laplaceM;

    std::array<Matrix,2> m_centered;
    dg::Advection< Geometry, Matrix, Container> m_adv;
    dg::PCG<Container> m_pcg;
    dg::Extrapolation<Container> m_extra;
    dg::Helmholtz<Geometry, Matrix, Container> m_helmholtz;

    Parameters m_p;

    unsigned m_ncalls = 0;

};

template< class Geometry, class M, class Container>
Explicit< Geometry, M, Container>::Explicit( const Geometry& grid, const Parameters& p ):
    m_chi( evaluate( dg::zero, grid)), m_uE2(m_chi),
    m_phi( m_chi), m_dxphi(m_phi), m_dyphi( m_phi), m_vx(m_phi), m_vy(m_phi),
    m_laplaceM( grid,  p.diff_dir),
    m_adv( grid),
    m_pcg( m_phi, grid.size()),
    m_extra( 2, m_phi),
    m_helmholtz( -1., {grid, dg::centered}),
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
    
    //need to compute m_phi from y here!!! (y = phi - lap phi)
    m_extra.extrapolate( t, m_phi);
    m_pcg.solve(m_helmholtz, m_phi, y,
                m_helmholtz.precond(), m_helmholtz.weights(), m_p.eps_gamma);
    m_extra.update( t, m_phi);
    //dg::blas2::symv(m_laplaceM, m_phi, m_chi); //alternative variante für chi
    dg::blas1::axpby( -1., m_phi, 1., y, m_chi); //chi = - lap \phi
    //compute derivatives
    dg::blas2::symv( m_centered[0], m_phi, m_dxphi);
    dg::blas2::symv( m_centered[1], m_phi, m_dyphi);

    dg::blas1::axpby(-1., m_dyphi, 0., m_vx); //compute ExB velocities v = uE = (-dy phi, dx phi)
    dg::blas1::axpby(1., m_dxphi, 0., m_vy);
    m_adv.upwind( -1., m_vx, m_vy, m_chi, 0., yp); // yp = uE.nabla(lap phi)

    //gradient terms
    dg::blas1::axpby( -1./m_p.Ln, m_dyphi, 1., yp); //Ln is background gradient length in units of rho_s
    
    //compute uE2
    m_laplaceM.variation(m_phi, m_uE2);

    return;
}

}//namespace dg
