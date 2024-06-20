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
    const Container& Ni() const { return m_Ni;}
    const Container& uE2() const { return m_uE2;}        

    dg::Elliptic<Geometry, Matrix, Container>& laplacianM( ) {
        return m_laplaceM;
    }

    unsigned ncalls() const{ return m_ncalls;}

    void operator()( double t, const Container & y,
            Container & yp);

    
  private:
    //use chi and omega as helpers to compute square velocity in uE2
    Container m_chi, m_uE2, m_invgamma1phi, m_invgamma0ni, m_Ni;
    Container m_phi, m_psi, m_dxphi, m_dyphi, m_uex, m_uey, m_ustx, m_usty, m_nomi, m_sinedrive;

    //matrices and solvers
    dg::Elliptic<Geometry, Matrix, Container> m_laplaceM;

    std::array<Matrix,2> m_centered;
    dg::Advection< Geometry, Matrix, Container> m_adv;
    dg::PCG<Container> m_pcg;
    dg::Extrapolation<Container> m_extra;
    dg::Extrapolation<Container> m_extra2;
    dg::Extrapolation<Container> m_extra3;
    dg::Helmholtz<Geometry, Matrix, Container> m_helmholtz;
    dg::Helmholtz<Geometry, Matrix, Container> m_invgamma0;
    dg::Helmholtz<Geometry, Matrix, Container> m_invgamma1;
    dg::Helmholtz<Geometry, Matrix, Container> m_invgamma2; //invgamma2 = 1-(taui+1)lap
    dg::HelmholtzLN<Geometry, Matrix, Container> m_helmholtzLN;
    dg::SinY sinedrive;

    Parameters m_p;

    unsigned m_ncalls = 0;

};

template< class Geometry, class M, class Container>
Explicit< Geometry, M, Container>::Explicit( const Geometry& grid, const Parameters& p ):
    m_chi( evaluate( dg::zero, grid)), m_uE2(m_chi), m_invgamma1phi(m_chi), m_invgamma0ni(m_chi), m_Ni(m_chi),
    m_phi( m_chi), m_psi(m_chi), m_dxphi(m_phi), m_dyphi( m_phi), m_uex(m_phi), m_uey(m_phi), m_ustx(m_phi), m_usty(m_phi), m_nomi(m_phi),
    m_laplaceM( grid,  p.diff_dir),
    m_adv( grid),
    m_pcg( m_phi, grid.size()),
    m_extra( 2, m_phi),
    m_extra2( 2, m_invgamma1phi),
    m_extra3( 2, m_psi),
    m_helmholtz( -1., {grid, dg::centered}),
    m_invgamma0( -p.taui, {grid, dg::centered}),
    m_invgamma1( -p.taui/2., {grid, dg::centered}),
    m_invgamma2( -p.taui-1, {grid, dg::centered}),
    m_helmholtzLN( -1., {grid, dg::centered}),
    sinedrive(m_p.sinedrive_amp, 0., m_p.sinedrive_k*2.*M_PI/m_p.ly),
    m_p(p)   
{
    m_centered = {dg::create::dx( grid, m_p.bcx),
                  dg::create::dy( grid, m_p.bcy)};
    m_sinedrive = dg::evaluate(sinedrive, grid);
}


template< class G, class M, class Container>
void Explicit<G, M, Container>::operator()( double t,
        const Container & y, Container & yp)
{
    m_ncalls ++ ;
    if (m_p.model == "standardCHM"){
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

        dg::blas1::axpby(-1., m_dyphi, 0., m_uex); //compute ExB velocities v = uE = (-dy phi, dx phi)
        dg::blas1::axpby(1., m_dxphi, 0., m_uey);
        m_adv.upwind( -1., m_uex, m_uey, m_chi, 0., yp); // yp = uE.nabla(lap phi)

        //gradient terms
        dg::blas1::axpby( -1./m_p.Ln, m_dyphi, 1., yp); //Ln is background gradient length in units of rho_s
        
        //compute uE2
        m_laplaceM.variation(m_phi, m_uE2);
    }
    else if (m_p.model == "boussinesq"){
        ///############### invert y here...
        m_extra.extrapolate( t, m_phi);
        m_pcg.solve(m_helmholtzLN, m_phi, y,
                    m_helmholtzLN.precond(), m_helmholtzLN.weights(), m_p.eps_gamma);
        m_extra.update( t, m_phi);
        //dg::blas2::symv( m_laplaceM, m_phi, m_chi);
        dg::blas1::axpby( 1., m_phi, -1., y, m_chi); //chi = lap \phi (=v)

        //compute derivatives
        dg::blas2::symv( m_centered[0], m_phi, m_dxphi);
        dg::blas2::symv( m_centered[1], m_phi, m_dyphi);
        dg::blas1::axpby(-1., m_dyphi, 0., m_uex); //compute ExB velocities v = uE = (-dy phi, dx phi)
        dg::blas1::axpby(1., m_dxphi, 0., m_uey);
        m_laplaceM.variation(m_phi, m_uE2); //compute uE2
        dg::blas2::symv( 0.5, m_centered[1], m_uE2, 0., m_ustx); //compute advection velocities Ust
        dg::blas2::symv(-0.5, m_centered[0], m_uE2, 0., m_usty);

        //now start combining RHS terms:
        m_adv.upwind( -1, m_ustx, m_usty, m_phi, 0., yp);
        dg::blas1::axpby(1./m_p.Ln, m_ustx, 1., yp);

        dg::blas1::transform( m_chi, m_nomi, dg::PLUS<double>(1.)); //nomi=1 - lap phi
        dg::blas1::pointwiseDivide( 1., yp, m_nomi, 0., yp); //yp = yp/(1-lap phi)
        
        m_adv.upwind( -1., m_uex, m_uey, y, 1., yp);
        //m_adv.upwind( -1., m_uex, m_uey, m_phi, 1., yp);
        //m_adv.upwind( -1., m_uex, m_uey, m_chi, 1., yp);
        dg::blas1::axpby( -1./m_p.Ln, m_dyphi, 1., yp);

    }
    else if (m_p.model == "FLRapprox"){
         dg::blas1::copy(y, m_Ni);
        //need to compute m_phi from y here!!! (y = Ni)
        // Gamma2 Ni =: psi
        m_extra3.extrapolate( t, m_psi);
        m_pcg.solve(m_invgamma2, m_psi, y,
                    m_invgamma2.precond(), m_invgamma2.weights(), m_p.eps_gamma);
        m_extra.update( t, m_psi);

        dg::blas2::symv(m_invgamma1, m_psi, m_phi);

       

        //compute derivatives
        dg::blas2::symv( m_centered[0], m_psi, m_dxphi);
        dg::blas2::symv( m_centered[1], m_psi, m_dyphi);

        dg::blas1::axpby(-1., m_dyphi, 0., m_uex); //compute ExB velocities v = uE = (-dy phi, dx phi)
        dg::blas1::axpby(1., m_dxphi, 0., m_uey);
        m_adv.upwind( -1., m_uex, m_uey, y, 0., yp); // yp = uE.nabla(lap phi)

        //gradient terms
        dg::blas1::axpby( -1., m_dyphi, 1., yp); //Ln is background gradient length in units of rho_s
        
        //compute uE2 and vorticity=lap phi
        m_laplaceM.variation(m_phi, m_uE2);
        dg::blas2::symv( m_laplaceM, m_phi, m_chi);

    }
    else if(m_p.model == "FLR"){
        dg::blas1::copy(y, m_Ni);
        //need to compute m_phi from y here!!! (y = Ni)
        dg::blas2::symv(m_invgamma0, y, m_invgamma0ni);
        m_extra2.extrapolate(t, m_invgamma1phi);
        m_pcg.solve(m_invgamma2, m_invgamma1phi, m_invgamma0ni, m_invgamma2.precond(), m_invgamma2.weights(), m_p.eps_gamma); 
        m_extra2.update( t, m_invgamma1phi);

        //now we have to solve invgamma1phi for phi
        m_extra.extrapolate( t, m_phi);
        m_pcg.solve(m_invgamma1, m_phi, m_invgamma1phi,
                    m_invgamma1.precond(), m_invgamma1.weights(), m_p.eps_gamma);
        m_extra.update( t, m_phi);

        //now add a driving force to phi
        if (m_p.sinedrive_amp !=0.0){
            dg::blas1::axpby(1., m_sinedrive, 1., m_phi);
        }
        
        //now compute psi out of phi, psi = gamma1 phi
        m_extra3.extrapolate( t, m_psi);
        m_pcg.solve(m_invgamma1, m_psi, m_phi,
                    m_invgamma1.precond(), m_invgamma1.weights(), m_p.eps_gamma);
        m_extra.update( t, m_psi);


        //compute derivatives
        dg::blas2::symv( m_centered[0], m_psi, m_dxphi);
        dg::blas2::symv( m_centered[1], m_psi, m_dyphi);

        dg::blas1::axpby(-1., m_dyphi, 0., m_uex); //compute ExB velocities v = uE = (-dy phi, dx phi)
        dg::blas1::axpby(1., m_dxphi, 0., m_uey);
        m_adv.upwind( -1., m_uex, m_uey, y, 0., yp); // yp = uE.nabla(lap phi)

        //gradient terms
        dg::blas1::axpby( -1., m_dyphi, 1., yp); //Ln is background gradient length in units of rho_s
        
        //compute uE2 and vorticity=lap phi
        m_laplaceM.variation(m_phi, m_uE2);
        dg::blas2::symv( m_laplaceM, m_phi, m_chi);

    }

    return;
}

}//namespace dg
