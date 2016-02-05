#pragma once

#include "dg/backend/grid.h"
#include "dg/backend/gridX.h"
#include "dg/nullstelle.h"

#include "orthogonal.h"



namespace solovev
{

namespace detail
{

struct XPointer
{
    XPointer( const GeomParameters& gp): fieldRZtau_(gp), psip_(gp){
        HessianRZtau hessianRZtau(gp);
        R_X = gp.R_0-1.1*gp.triangularity*gp.a;
        Z_X = -1.1*gp.elongation*gp.a;
        thrust::host_vector<double> X(2,0), XN(X);
        X[0] = R_X, X[1] = Z_X;
        for( unsigned i=0; i<3; i++)
        {
            hessianRZtau.newton_iteration( X, XN);
            XN.swap(X);
        }
        R_X = X[0], Z_X = X[1];
        std::cout << "X-point set at "<<R_X<<" "<<Z_X<<"\n";
        R_i[0] = R_X + 1, Z_i[0] = Z_X;
        R_i[1] = R_X    , Z_i[1] = Z_X + 1;
        R_i[2] = R_X - 1, Z_i[2] = Z_X;
        R_i[3] = R_X    , Z_i[3] = Z_X - 1;
    }
    void set_quadrant( int quad){quad_ = quad;}
    double operator()( double x) const
    {
        thrust::host_vector<double> begin(2), end(2), end_old(2);
        begin[0] = R_i[quad_], begin[1] = Z_i[quad_];
        double eps = 1e10, eps_old = 2e10;
        unsigned N=10;
        if( quad_ == 0 || quad_ == 2) { begin[1] += x;}
        else if( quad_ == 1 || quad_ == 3) { begin[0] += x;}

        double psi0 = psip_(begin[0], begin[1]);
        while( (eps < eps_old || eps > 1e-4 ) && eps > 1e-7)
        {
            eps_old = eps; end_old = end;
            N*=2; dg::stepperRK17( fieldRZtau_, begin, end, psi0, 0, N);
            eps = sqrt( (end[0]-end_old[0])*(end[0]-end_old[0]) + (end[1]-end_old[1])*(end[1]-end_old[1]));
            if( isnan(eps)) { eps = eps_old/2.; end = end_old; }
        }
        //std::cout<< end_old[0]  - R_X << "\t";
        //std::cout<< end_old[1]  - Z_X << std::endl;
        if( quad_ == 0 || quad_ == 2){ return end_old[1] - Z_X;}
        return end_old[0] - R_X;
    }
    void point( double& R, double& Z, double x)
    {
        if( quad_ == 0 || quad_ == 2){ R = R_i[quad_], Z= Z_i[quad_] +x;}
        else if (quad_ == 1 || quad_ == 3) { R = R_i[quad_] + x, Z = Z_i[quad_];}
    }

    private:
    int quad_;
    FieldRZtau fieldRZtau_;
    Psip psip_;
    double R_X, Z_X;
    double R_i[4], Z_i[4];
};

//This leightweights struct and its methods finds the initial R and Z values and the coresponding f(\psi) as 
//good as it can, i.e. until machine precision is reached
struct FpsiX
{
    FpsiX( const GeomParameters& gp): 
        gp_(gp), fieldRZYT_(gp), fieldRZYZ_(gp), fieldRZtau_(gp), 
        xpointer_(gp), hessianRZtau_(gp), minimalCurve_(gp)

    {
        //////////////////////////////////////////////
        Psip psip(gp_);
        PsipR psipR(gp_);
        PsipZ psipZ(gp_);
        double R_ini = gp.R_0-1.1*gp.triangularity*gp.a;
        double Z_ini = -1.1*gp.elongation*gp.a;
        thrust::host_vector<double> X(2,0), XN(X);
        X[0] = R_ini, X[1] = Z_ini;
        for( unsigned i=0; i<3; i++)
        {
            hessianRZtau_.newton_iteration( X, XN);
            XN.swap(X);
        }
        R_ini = X[0], Z_ini = X[1];
        std::cout << "X-point set at "<<R_ini<<" "<<Z_ini<<"\n";
        //std::cout << "psi at X-point is "<<psip(R_ini, Z_ini)<<"\n";
        //std::cout << "gradient at X-point is "<<psipR(R_ini, Z_ini)<<" "<<psipZ(R_ini, Z_ini)<<"\n";
        //find four points; one in each quadrant
        hessianRZtau_.set_norm( false);
        minimalCurve_.set_norm( false);
        for( int i=0; i<4; i++)
        {
            hessianRZtau_.set_quadrant( i);
            unsigned N = 50;
            //thrust::host_vector<double> begin2d( 2, 0), end2d( begin2d), end2d_old(begin2d); 
            thrust::host_vector<double> begin2d( 4, 0), end2d( begin2d), end2d_old(begin2d); 
            begin2d[0] = end2d[0] = end2d_old[0] = R_ini;
            begin2d[1] = end2d[1] = end2d_old[1] = Z_ini;
            hessianRZtau_(begin2d, end2d); //find eigenvector
            begin2d[2] = end2d[0], begin2d[3] = end2d[1];
            double eps = 1e10, eps_old = 2e10;
            while( eps < eps_old && N<1e6 && eps > 1e-15)
            {
                //remember old values
                eps_old = eps; end2d_old = end2d;
                //compute new values
                N*=2;
                dg::stepperRK17( hessianRZtau_, begin2d, end2d, 0., 5., N);
                //dg::stepperRK17( minimalCurve_, begin2d, end2d, 0., 2., N);
                eps = sqrt( (end2d[0]-end2d_old[0])*(end2d[0]-end2d_old[0]) + (end2d[1]-end2d_old[1])*(end2d[1]-end2d_old[1]));
            }
            R_i[i] = end2d_old[0], Z_i[i] = end2d_old[1]; 
            vR_i[i] = end2d_old[2], vZ_i[i] = end2d_old[3];
            //std::cout << "Found the point "<<R_i[i]<<" "<<Z_i[i]<<" "<<psip(R_i[i], Z_i[i])<<"\n";
        }
        hessianRZtau_.set_norm( true);
        minimalCurve_.set_norm( true);
    }
    //finds the two starting points for the integration in y direction
    void find_initial( double psi, double* R_0, double* Z_0) 
    {
        solovev::Psip psip(gp_);
        thrust::host_vector<double> begin2d( 4, 0), end2d( begin2d), end2d_old(begin2d); 
        //std::cout << "In init function\n";
        //std::cout << "psi is "<<psi<<"\n";
        if( psi < 0)
        {
            for( unsigned i=0; i<2; i++)
            {
                unsigned N = 50;
                hessianRZtau_.set_quadrant( 1+2*i);
                begin2d[0] = end2d[0] = end2d_old[0] = R_i[1+2*i];
                begin2d[1] = end2d[1] = end2d_old[1] = Z_i[1+2*i];
                begin2d[2] = end2d[2] = end2d_old[2] = vR_i[1+2*i];
                begin2d[3] = end2d[3] = end2d_old[3] = vZ_i[1+2*i];
                double eps = 1e10, eps_old = 2e10;
                while( (eps < eps_old || eps > 1e-7) && N<1e6 && eps > 1e-11)
                {
                    //remember old values
                    eps_old = eps;
                    end2d_old = end2d;
                    //compute new values
                    N*=2;
                    dg::stepperRK17( hessianRZtau_, begin2d, end2d, psip(R_i[1+2*i], Z_i[1+2*i]), psi, N);
                    //dg::stepperRK17( minimalCurve_, begin2d, end2d, psip(R_i[1+2*i], Z_i[1+2*i]), psi, N);
                    eps = sqrt( (end2d[0]-end2d_old[0])*(end2d[0]-end2d_old[0]) + (end2d[1]-end2d_old[1])*(end2d[1]-end2d_old[1]));
                }
                //remember last call
                R_i[1+2*i] = R_0[i] = end2d_old[0], Z_i[1+2*i] = Z_0[i] = end2d_old[1];
                vR_i[1+2*i] = end2d_old[2], vZ_i[1+2*i] = end2d_old[3];
            }
        }
        else
        {
            for( unsigned i=0; i<2; i++)
            {
                unsigned N=50;
                hessianRZtau_.set_quadrant( 2*i);
                begin2d[0] = end2d[0] = end2d_old[0] = R_i[2*i];
                begin2d[1] = end2d[1] = end2d_old[1] = Z_i[2*i];
                begin2d[2] = end2d[2] = end2d_old[2] = vR_i[2*i];
                begin2d[3] = end2d[3] = end2d_old[3] = vZ_i[2*i];
                double eps = 1e10, eps_old = 2e10;
                while( (eps < eps_old || eps > 1e-7) && N<1e6 && eps > 1e-11)
                {
                    //remember old values
                    eps_old = eps;
                    end2d_old = end2d;
                    //compute new values
                    N*=2;
                    dg::stepperRK17( hessianRZtau_, begin2d, end2d, psip(R_i[2*i], Z_i[2*i]), psi, N);
                    //dg::stepperRK17( minimalCurve_, begin2d, end2d, psip(R_i[2*i], Z_i[2*i]), psi, N);
                    eps = sqrt( (end2d[0]-end2d_old[0])*(end2d[0]-end2d_old[0]) + (end2d[1]-end2d_old[1])*(end2d[1]-end2d_old[1]));
                }
                R_0[i] = end2d_old[0], Z_0[i] = end2d_old[1];
                //remember last call
                R_i[2*i] = R_0[i] = end2d_old[0], Z_i[2*i] = Z_0[i] = end2d_old[1];
                vR_i[2*i] = end2d_old[2], vZ_i[2*i] = end2d_old[3];
            }
        }

    }

    //compute f for a given psi between psi0 and psi1
    double construct_f( double psi, double* R_i, double* Z_i) 
    {
        dg::Timer t;
        t.tic();
        find_initial( psi, R_i, Z_i);
        t.toc();
        solovev::Psip psip( gp_);
        
        //std::cout << "find_initial took "<<t.diff()<< "s\n";
        //t.tic();
        //std::cout << "Begin error "<<eps_old<<" with "<<N<<" steps\n";
        //std::cout << "In Stepper function:\n";
        //double y_old=0;
        thrust::host_vector<double> begin( 3, 0), end(begin), end_old(begin);
        begin[0] = R_i[0], begin[1] = Z_i[0];
        //std::cout << begin[0]<<" "<<begin[1]<<" "<<begin[2]<<"\n";
        double eps = 1e10, eps_old = 2e10;
        unsigned N = 32; 
        //double y_eps;
        while( (eps < eps_old || eps > 1e-7) && N < 1e6)
        {
            //remember old values
            eps_old = eps, end_old = end; //y_old = end[2];
            //compute new values
            N*=2;
            if( psi < 0)
            {
                dg::stepperRK17( fieldRZYT_, begin, end, 0., 2.*M_PI, N);
                //std::cout << "result is "<<end[0]<<" "<<end[1]<<" "<<end[2]<<"\n";
                eps = sqrt( (end[0]-begin[0])*(end[0]-begin[0]) + (end[1]-begin[1])*(end[1]-begin[1]));
            }
            else
            {
                dg::stepperRK17( fieldRZYZ_, begin, end, begin[1], 0., N);
                thrust::host_vector<double> temp(end);
                dg::stepperRK17( fieldRZYT_, temp, end, 0., M_PI, N);
                temp = end; //temp[1] should be 0 now
                dg::stepperRK17( fieldRZYZ_, temp, end, temp[1], Z_i[1], N);
                eps = sqrt( (end[0]-R_i[1])*(end[0]-R_i[1]) + (end[1]-Z_i[1])*(end[1]-Z_i[1]));
            }
            if( isnan(eps)) { eps = eps_old/2.; end = end_old; 
                //std::cerr << "\t nan! error "<<eps<<"\n";
            } //near X-point integration can go wrong
            //y_eps = sqrt( (y_old - end[2])*(y_old-end[2]));
            //std::cout << "error "<<eps<<" with "<<N<<" steps| psip "<<psip(end[0], end[1])<<"\n";
            //std::cout <<"error in y is "<<y_eps<<"\n";
        }
        //std::cout << "\t error "<<eps<<" with "<<N<<" steps| err psip "<<fabs( psip(end[0], end[1]) - psi )/psi<<"\n";
        double f_psi = 2.*M_PI/end_old[2];
        //t.toc();
        //std::cout << "Finding f took "<<t.diff()<<"s\n";
        return f_psi;
        //return 1./f_psi;
    }
    double operator()( double psi)
    {
        double R_0[2], Z_0[2]; 
        return construct_f( psi, R_0, Z_0);
    }

    /**
     * @brief This function computes the integral x_0 = \int_{\psi}^{0} f(\psi) d\psi to machine precision
     *
     * @return x0
     */
    double find_x( double psi ) 
    {
        unsigned P=6;
        double x0 = 0, x0_old = 0;
        double eps=1e10, eps_old=2e10;
        //std::cout << "In x1 function\n";
        while( (eps < eps_old||eps>1e-7) && P < 20 )
        {
            eps_old = eps; x0_old = x0;
            P+=2;
            dg::Grid1d<double> grid( 0, 1, P, 1);
            if( psi>0)
            {
                dg::Grid1d<double> grid1( 0, psi, P, 1);
                grid = grid1;
            }
            else 
            {
                dg::Grid1d<double> grid2( psi, 0, P, 1);
                grid = grid2;
            }
            thrust::host_vector<double> psi_vec = dg::evaluate( dg::coo1, grid);
            thrust::host_vector<double> f_vec(grid.size(), 0);
            thrust::host_vector<double> w1d = dg::create::weights(grid);
            for( unsigned i=0; i<psi_vec.size(); i++)
            {
                f_vec[i] = this->operator()( psi_vec[i]);
                //std::cout << " "<<f_vec[i]<<"\n";
            }
            if( psi < 0)
                x0 = dg::blas1::dot( f_vec, w1d);
            else
                x0 = -dg::blas1::dot( f_vec, w1d);

            eps = fabs((x0 - x0_old)/x0);
            if( isnan(eps)) { std::cerr << "Attention!!\n"; eps = eps_old -1e-15; x0 = x0_old;} //near X-point integration can go wrong
            //std::cout << "X = "<<x0<<" rel. error "<<eps<<" with "<<P<<" polynomials\n";
        }
        return x0_old;

    }


    //compute the vector of r and z - values that form one psi surface
    //calls construct_f to find f and the starting point and then just integrates
    //the field-line and metric from 0 to 2pi in y
    void compute_rzy( double psi, unsigned n, unsigned N, double fy,
            thrust::host_vector<double>& r, 
            thrust::host_vector<double>& z, 
            thrust::host_vector<double>& yr, 
            thrust::host_vector<double>& yz,  
            double* R_0, double* Z_0, double& f ) 
    {
        assert( psi < 0);
        //find start points for first psi surface
        thrust::host_vector<double> begin( 2, 0), end(begin), temp(begin), end_old(end);
        double R_init[2], Z_init[2];
        for( unsigned i=0; i<2; i++)
        {
            xpointer_.set_quadrant( 2*i+1);
            double x_min = -1, x_max = 1;
            dg::bisection1d( xpointer_, x_min, x_max, 1e-6);
            xpointer_.point( R_init[i], Z_init[i], (x_min+x_max)/2.);
            std::cout << "Found initial point! "<<R_init[i]<<" "<<Z_init[i]<<"\n";

            begin[0] = R_init[i], begin[1] = Z_init[i]; end = begin;
            unsigned steps = 1;
            double eps = 1e10, eps_old=2e10;
            while( (eps < eps_old||eps > 1e-7) && eps > 1e-11)
            {
                eps_old = eps; end_old = end;
                dg::stepperRK17( fieldRZtau_, begin, end, Psip(gp_)(R_init[i], Z_init[i]), psi, steps);
                eps = sqrt( (end[0]-end_old[0])*(end[0]- end_old[0]) + (end[1]-end_old[1])*(end[1]-end_old[1]));
                //std::cout << "rel. error is "<<eps<<" with "<<steps<<" steps\n";
                if( isnan(eps)) { eps = eps_old/2.; end = end_old; }
                steps*=2;
            }
            std::cout << "Found initial point "<<end_old[0]<<" "<<end_old[1]<<"\n";
            R_init[i] = begin[0] = end_old[0], Z_init[i] = begin[1] = end_old[1];
        }
        ///////////////////////////now find y coordinate line//////////////
        dg::GridX1d g1d( -fy*2.*M_PI/(1.-2.*fy), 2*M_PI+fy*2.*M_PI/(1.-2.*fy), fy, n, N, dg::DIR);
        thrust::host_vector<double> y_vec = dg::evaluate( dg::coo1, g1d);
        thrust::host_vector<double> r_old(g1d.size(), 0), r_diff( r_old);
        thrust::host_vector<double> z_old(g1d.size(), 0), z_diff( z_old);
        const thrust::host_vector<double> w1d = dg::create::weights( g1d);
        r.resize( g1d.size()), z.resize(r.size()), yr.resize(r.size()), yz.resize(r.size());
        solovev::orthogonal::FieldRZY fieldRZY(gp_);
        //now compute f and starting values 
        double f_psi = construct_f( psi, R_0, Z_0);
        fieldRZY.set_f(f_psi);
        //fieldRZY.set_f(1./f_psi);
        unsigned steps = 1; double eps = 1e10, eps_old=2e10;
        while( (eps < eps_old||eps > 1e-7) && eps > 1e-11)
        {
            eps_old = eps, r_old = r, z_old = z;

            ////////////////////////bottom left region/////////////////////
            begin[0] = R_init[1], begin[1] = Z_init[1];
            dg::stepperRK17( fieldRZY, begin, end, 0, y_vec[n*g1d.outer_N()-1], steps);
            //r[n*g1d.outer_N()-1] = end[0], z[n*g1d.outer_N()-1] = end[1];
            r[n*g1d.outer_N()-1] = end[0], z[n*g1d.outer_N()-1] = end[1];
            for( int i=n*g1d.outer_N()-2; i>=0; i--)
            {
                temp = end;
                dg::stepperRK17( fieldRZY, temp, end, y_vec[i+1], y_vec[i], steps);
                //r[n*g1d.outer_N()-1-i] = end[0], z[n*g1d.outer_N()-1-i] = end[1];
                r[i] = end[0], z[i] = end[1];
            }
            ////////////////middle region///////////////////////////
            begin[0] = R_init[0], begin[1] = Z_init[0];
            dg::stepperRK17( fieldRZY, begin, end, 0, y_vec[n*g1d.outer_N()], steps);
            r[n*g1d.outer_N()] = end[0], z[n*g1d.outer_N()] = end[1];
            //std::cout <<end[0]<<" "<< end[1] <<"\n";
            for( unsigned i=n*g1d.outer_N()+1; i<n*(g1d.outer_N()+g1d.inner_N()); i++)
            {
                temp = end;
                dg::stepperRK17( fieldRZY, temp, end, y_vec[i-1], y_vec[i], steps);
                r[i] = end[0], z[i] = end[1];
            }
            temp = end;
            dg::stepperRK17( fieldRZY, temp, end, y_vec[n*(g1d.outer_N()+g1d.inner_N())-1], 2.*M_PI, steps);
            eps = sqrt( (end[0]-R_init[0])*(end[0]-R_init[0]) + (end[1]-Z_init[0])*(end[1]-Z_init[0]));
            std::cout << "abs. error is "<<eps<<" with "<<steps<<" steps\n";
            ////////////////////bottom right region

            begin[0] = R_init[1], begin[1] = Z_init[1];
            dg::stepperRK17( fieldRZY, begin, end, 2.*M_PI, y_vec[n*(g1d.outer_N()+g1d.inner_N())], steps);
            r[n*(g1d.outer_N()+g1d.inner_N())] = end[0], z[n*(g1d.outer_N()+g1d.inner_N())] = end[1];
            for( unsigned i=n*(g1d.outer_N()+g1d.inner_N())+1; i<n*g1d.N(); i++)
            {
                temp = end;
                dg::stepperRK17( fieldRZY, temp, end, y_vec[i-1], y_vec[i], steps);
                r[i] = end[0], z[i] = end[1];
            }
            //compute error in R,Z only
            dg::blas1::axpby( 1., r, -1., r_old, r_diff);
            dg::blas1::axpby( 1., z, -1., z_old, z_diff);
            double er = dg::blas2::dot( r_diff, w1d, r_diff);
            double ez = dg::blas2::dot( z_diff, w1d, z_diff);
            double ar = dg::blas2::dot( r, w1d, r);
            double az = dg::blas2::dot( z, w1d, z);
            eps =  sqrt( er + ez)/sqrt(ar+az);
            std::cout << "rel. error is "<<eps<<" with "<<steps<<" steps\n";
            steps*=2;
        }
        r = r_old, z = z_old;
        f = f_psi;
        PsipR psipR_(gp_);
        PsipZ psipZ_(gp_);
        for( unsigned i=0; i<r.size(); i++)
        {
            double psipR = psipR_( r[i], z[i]), psipZ = psipZ_( r[i], z[i]);
            //double psip2 = psipR*psipR+psipZ*psipZ;
            //yr[i] = psipZ*f/psip2;
            //yz[i] = -psipR*f/psip2;
            //yr[i] = psipZ*f/sqrt(psip2);
            //yz[i] = -psipR*f/sqrt(psip2);
            //yr[i] = psipZ/f/sqrt(psip2);
            //yz[i] = -psipR/f/sqrt(psip2);
            yr[i] = psipZ*f;
            yz[i] = -psipR*f;
        }

    }
    private:
    const GeomParameters gp_;
    const solovev::orthogonal::FieldRZYT fieldRZYT_;
    const solovev::orthogonal::FieldRZYZ fieldRZYZ_;
    const FieldRZtau fieldRZtau_;
    XPointer xpointer_;
    HessianRZtau hessianRZtau_;
    MinimalCurve minimalCurve_;
    double R_i[4], Z_i[4], vR_i[4], vZ_i[4];

};
//This struct computes -2pi/f with a fixed number of steps for all psi
struct XFieldFinv
{
    XFieldFinv( const GeomParameters& gp, unsigned N_steps = 500): 
        fpsi_(gp), fieldRZYT_(gp), fieldRZYZ_(gp) , N_steps(N_steps),
        R_0_(gp.R_0), psipR_(gp), psipZ_(gp),
        psipRR_(gp), psipZZ_(gp), psipRZ_(gp)
            { xAtOne_ = fpsi_.find_x(0.1); }
    void operator()(const std::vector<thrust::host_vector<double> >& y, std::vector<thrust::host_vector<double> >& yp) 
    { 
        //y[0] = R, y[1] = Z , y[2] = g, y[3] = yr, y[4] = yz
        //this->operator()( y[2], yp[2]);

        unsigned size = y[0].size();
        double psipR, psipZ, psipRR, psipRZ, psipZZ, psip2;
        for( unsigned i=0; i<size; i++)
        {
            psipR = psipR_(y[0][i], y[1][i]), psipZ = psipZ_(y[0][i], y[1][i]);
            psipRR = psipRR_(y[0][i], y[1][i]), psipRZ = psipRZ_(y[0][i], y[1][i]), psipZZ = psipZZ_(y[0][i], y[1][i]);
            psip2 = psipR*psipR+psipZ*psipZ;
            yp[0][i] = psipR/psip2;
            yp[1][i] = psipZ/psip2;
            //yp[2][i] = y[2][i]/psip2*( 2./psip2*( psipR*psipR*psipRR +psipZ*psipZ*psipZZ+2.*psipZ*psipR*psipRZ )  -(psipRR+psipZZ) );
            //yp[2][i] = y[2][i]/psip2*( 1./psip2/sqrt(psip2)*( psipR*psipR*psipRR +psipZ*psipZ*psipZZ+2.*psipZ*psipR*psipRZ )  -(psipRR+psipZZ) );//g/gradpsi^1/2
            yp[2][i] = y[2][i]/psip2*( -(psipRR+psipZZ) );//g
            yp[3][i] = 1./psip2 *( -psipRR*y[3][i] - psipRZ*y[4][i]);
            yp[4][i] = 1./psip2 *( -psipRZ*y[3][i] - psipZZ*y[4][i]);
        }
    }
    void operator()(const thrust::host_vector<double>& psi, thrust::host_vector<double>& fpsiM) 
    { 
        thrust::host_vector<double> begin( 3, 0), end(begin), end_old(begin);
        double R_i[2], Z_i[2];
        dg::Timer t;
        t.tic();
        fpsi_.find_initial( psi[0], R_i, Z_i);
        t.toc();
        //std::cout << "find_initial took "<<t.diff()<< "s\n";
        t.tic();
        begin[0] = R_i[0], begin[1] = Z_i[0];
        unsigned N = N_steps;
        if( psi[0] < -1. && psi[0] > -2.) N*=2;
        if( psi[0] < 0 && psi[0] > -1.) N*=10;
        if( psi[0] <0  )
            dg::stepperRK17( fieldRZYT_, begin, end, 0., 2.*M_PI, N);
        else
        {
            dg::stepperRK17( fieldRZYZ_, begin, end, begin[1], 0., N);
            thrust::host_vector<double> temp(end);
            dg::stepperRK17( fieldRZYT_, temp, end, 0., M_PI, N/2);
            temp = end; //temp[1] should be 0 now
            dg::stepperRK17( fieldRZYZ_, temp, end, temp[1], Z_i[1], N);
        }
        //eps = sqrt( (end[0]-begin[0])*(end[0]-begin[0]) + (end[1]-begin[1])*(end[1]-begin[1]));
        fpsiM[0] = - end[2]/2./M_PI;
        //fpsiM[0] = - 2.*M_PI/end[2];
        t.toc();
        //std::cout << "Finding f took "<<t.diff()<<"s\n";
        //std::cout <<"fpsiMinverse is "<<fpsiM[0]<<" "<<-1./fpsi_(psi[0])<<" "<<eps<<"\n";
    }
    double find_psi( double x)
    {
        assert( x > 0);
        //integrate from x0 to x, with psi(x0) = 1;
        double x0 = xAtOne_; 
        thrust::host_vector<double> begin( 1, 0.1), end(begin), end_old(begin);
        double eps = 1e10, eps_old = 2e10;
        unsigned N = 1;
        while( eps < eps_old && N < 1e6 &&  eps > 1e-9)
        {
            eps_old = eps, end_old = end; 
            N*=2; dg::stepperRK17( *this, begin, end, x0, x, N);
            eps = fabs( end[0]- end_old[0]);
            //std::cout << "\t error "<<eps<<" with "<<N<<" steps\n";
        }
        return end_old[0];
    }

    private:
    FpsiX fpsi_;
    solovev::orthogonal::FieldRZYT fieldRZYT_;
    solovev::orthogonal::FieldRZYZ fieldRZYZ_;
    thrust::host_vector<double> fpsi_neg_inv;
    unsigned N_steps;
    double xAtOne_;
    double R_0_;
    PsipR psipR_;
    PsipZ psipZ_;
    PsipRR psipRR_;
    PsipZZ psipZZ_;
    PsipRZ psipRZ_;
};
} //namespace detail

template< class container>
struct OrthogonalXGrid2d; 

/**
 * @brief A three-dimensional grid based on "almost-conformal" coordinates by Ribeiro and Scott 2010
 */
template< class container>
struct OrthogonalXGrid3d : public dg::GridX3d
{
    typedef dg::CurvilinearCylindricalTag metric_category;
    typedef OrthogonalXGrid2d<container> perpendicular_grid;

    /**
     * @brief Construct 
     *
     * @param gp The geometric parameters define the magnetic field
     * @param psi_0 lower boundary for psi
     * @param fx factor in x-direction
     * @param fy factor in y-direction
     * @param n The dG number of polynomials
     * @param Nx The number of points in x-direction
     * @param Ny The number of points in y-direction
     * @param Nz The number of points in z-direction
     * @param bcx The boundary condition in x (z is periodic)
     * @param bcy The boundary condition in y (z is periodic)
     */
    OrthogonalXGrid3d( GeomParameters gp, double psi_0, double fx, double fy, unsigned n, unsigned Nx, unsigned Ny, unsigned Nz, dg::bc bcx, dg::bc bcy): 
        dg::GridX3d( 0,1, -2.*M_PI*fy/(1.-2.*fy), 2.*M_PI*(1.+fy/(1.-2.*fy)), 0., 2*M_PI, fx, fy, n, Nx, Ny, Nz, bcx, bcy, dg::PER)
    { 
        assert( psi_0 < 0 );
        assert( gp.c[10] != 0);
        //construct x-grid in two parts
        solovev::detail::FpsiX fpsi(gp);
        std::cout << "FIND X FOR PSI_0\n";
        const double x_0 = fpsi.find_x(psi_0);
        const double x_1 = -fx/(1.-fx)*x_0;
        init_X_boundaries( x_0, x_1);
        //compute psi(x) for a grid on x 
        dg::Grid1d<double> g1d_( this->x0(), this->x1(), n, Nx, bcx);
        std::cout << "X0 is "<<x_0<<" and X1 is "<<x_1<<"\n";
        g1d_.display();
        thrust::host_vector<double> x_vec = dg::evaluate( dg::coo1, g1d_);
        thrust::host_vector<double> psi_x(n*Nx, 0), psi_old(psi_x), psi_diff( psi_old);
        f_x_.resize( psi_x.size());
        thrust::host_vector<double> w1d = dg::create::weights( g1d_);
        unsigned N = 1;
        std::cout << "In psi function:\n";
        double x0=this->x0(), x1 = x_vec[0];
        detail::XFieldFinv fpsiMinv_(gp, 500);
        const unsigned idx = inner_Nx()*this->n();
        const double psi_const = fpsiMinv_.find_psi( x_vec[idx]);
        double eps = 1e10;//, eps_old=2e10;
        //while( eps <  eps_old && N < 1e6)
        while( eps >  1e-8 && N < 1e6 )
        {
           // eps_old = eps; 
            psi_old = psi_x; 
            x0 = this->x0(), x1 = x_vec[0];

            thrust::host_vector<double> begin(1,psi_0), end(begin), temp(begin);
            dg::stepperRK6( fpsiMinv_, begin, end, x0, x1, N);
            psi_x[0] = end[0]; fpsiMinv_(end,temp); f_x_[0] = temp[0];
            for( unsigned i=1; i<idx; i++)
            {
                temp = end;
                x0 = x_vec[i-1], x1 = x_vec[i];
                dg::stepperRK6( fpsiMinv_, temp, end, x0, x1, N);
                psi_x[i] = end[0]; fpsiMinv_(end,temp); f_x_[i] = temp[0];
                //std::cout << "FOUND PSI "<<end[0]<<"\n";
            }
            end[0] = psi_const;
            //std::cout << "FOUND PSI "<<end[0]<<"\n";
            psi_x[idx] = end[0]; fpsiMinv_(end,temp); f_x_[idx] = temp[0];
            for( unsigned i=idx+1; i<g1d_.size(); i++)
            {
                temp = end;
                x0 = x_vec[i-1], x1 = x_vec[i];
                dg::stepperRK6( fpsiMinv_, temp, end, x0, x1, N);
                psi_x[i] = end[0]; fpsiMinv_(end,temp); f_x_[i] = temp[0];
                //std::cout << "FOUND PSI "<<end[0]<<"\n";
            }
            dg::blas1::axpby( 1., psi_x, -1., psi_old, psi_diff);
            eps = sqrt( dg::blas2::dot( psi_diff, w1d, psi_diff)/ dg::blas2::dot( psi_x, w1d, psi_x));
            psi_1_numerical_ = psi_0 + dg::blas1::dot( f_x_, w1d);

            //eps = fabs( psi_1_numerical-psi_1); 
            //std::cout << "Effective absolute Psi error is "<<psi_1_numerical-psi_1<<" with "<<N<<" steps\n"; 
            std::cout << "Effective Psi error is "<<eps<<" with "<<N<<" steps\n"; 
            std::cout << "psi 1               is "<<psi_1_numerical_<<"\n"; 
            N*=2;
        }
        construct_rz( gp, psi_0, psi_x);
        construct_metric();
    }
    const thrust::host_vector<double>& f()const{return f_;}
    const thrust::host_vector<double>& g()const{return g_;}
    const thrust::host_vector<double>& r()const{return r_;}
    const thrust::host_vector<double>& z()const{return z_;}
    const thrust::host_vector<double>& xr()const{return xr_;}
    const thrust::host_vector<double>& yr()const{return yr_;}
    const thrust::host_vector<double>& xz()const{return xz_;}
    const thrust::host_vector<double>& yz()const{return yz_;}
    const thrust::host_vector<double>& f_x()const{return f_x_;}
    thrust::host_vector<double> x()const{
        dg::Grid1d<double> gx( x0(), x1(), n(), Nx());
        return dg::create::abscissas(gx);}
    const container& g_xx()const{return g_xx_;}
    const container& g_yy()const{return g_yy_;}
    const container& g_xy()const{return g_xy_;}
    const container& g_pp()const{return g_pp_;}
    const container& vol()const{return vol_;}
    const container& perpVol()const{return vol2d_;}
    perpendicular_grid perp_grid() const { return OrthogonalXGrid2d<container>(*this);}
    const thrust::host_vector<double>& rx0()const{return r_x0;}
    const thrust::host_vector<double>& zx0()const{return z_x0;}
    const thrust::host_vector<double>& rx1()const{return r_x1;}
    const thrust::host_vector<double>& zx1()const{return z_x1;}
    double psi1()const{return psi_1_numerical_;}
    private:
    //call the construct_rzy function for all psi_x and lift to 3d grid
    //construct r,z,xr,xz,yr,yz,f_x
    void construct_rz( const GeomParameters& gp, double psi_0, thrust::host_vector<double>& x_vec)
    {
        //x_vec contains psi(x) values
        // convergence utilities
        dg::Grid1d<double> g1d_( this->x0(), this->x1(), this->n(), this->Nx(), this->bcx());
        unsigned Nx = this->n()*this->Nx(), Ny = this->n()*this->Ny();
        thrust::host_vector<double> psi_x(Nx, 0), psi_old(psi_x), psi_diff( psi_old);
        thrust::host_vector<double> w1d = dg::create::weights( g1d_);
        //fpsiMinv utilities
        detail::XFieldFinv fpsiMinv_(gp, 500);
        thrust::host_vector<double> rvec( Ny), zvec(Ny), yrvec(Ny), yzvec(Ny);
        std::vector<thrust::host_vector<double> > begin(5);
        double R0[2], Z0[2], f0;
        detail::FpsiX fpsi(gp);
        fpsi.compute_rzy( psi_0, this->n(), this->Ny(), this->fy(), rvec, zvec, yrvec, yzvec, R0, Z0, f0);
        thrust::host_vector<double> gvec(Ny, f0);
        //thrust::host_vector<double> gvec(Ny, 1./f0);
        begin[0] = rvec, begin[1] = zvec;
        begin[2] = gvec, begin[3] = yrvec, begin[4] = yzvec;
        //now we have the starting values of r, z, psi
        std::vector<thrust::host_vector<double> > end(begin), temp(begin);
        r_.resize(size()), z_.resize(size()), f_.resize(size()), g_.resize(size());
        yr_ = r_, yz_ = z_, xr_ = r_, xz_ = r_ ;
        std::cout << "In RZ  function:\n";
        double x0=this->x0(), x1 = x_vec[0];
        const unsigned idx = inner_Nx()*this->n();
        const double psi_const = fpsiMinv_.find_psi( x_vec[idx]);
        //while( eps <  eps_old && N < 1e6)
        PsipR psipR_(gp);
        PsipZ psipZ_(gp);
        double eps = 1e10;
        unsigned N=1; 
        while( eps >  1e-8 && N < 1e6 )
        {
           // eps_old = eps; 
            psi_old = psi_x; 
            //x0 = this->x0(), x1 = x_vec[0];
            x0 = psi_0, x1 = x_vec[0];
            //////////////////////////////////////////////////
            dg::stepperRK6( fpsiMinv_, begin, end, x0, x1, N);
            psi_x[0] = end[0][0]; //R
            for( unsigned j=0; j<Ny; j++)
            {
                 r_[j*Nx+0] = end[0][j],  z_[j*Nx+0]  = end[1][j];
                yr_[j*Nx+0] = end[3][j], yz_[j*Nx+0] = end[4][j];
                xr_[j*Nx+0] = psipR_(end[0][j], end[1][j])/f_x_[0]; 
                xz_[j*Nx+0] = psipZ_(end[0][j], end[1][j])/f_x_[0]; 
                 g_[j*Nx+0] = end[2][j], f_[j*Nx + 0] = -1./f_x_[0];
            }

            //////////////////////////////////////////////////
            for( unsigned i=1; i<g1d_.size(); i++)
            {
                temp = end;
                x0 = x_vec[i-1], x1 = x_vec[i];
                //////////////////////////////////////////////////
                dg::stepperRK6( fpsiMinv_, temp, end, x0, x1, N);
                psi_x[idx] = end[0][i]; //R
                for( unsigned j=0; j<Ny; j++)
                {
                     r_[j*Nx+i] = end[0][j],  z_[j*Nx+i] = end[1][j];
                    yr_[j*Nx+i] = end[3][j], yz_[j*Nx+i] = end[4][j];
                    xr_[j*Nx+i] = psipR_(end[0][j], end[1][j])/f_x_[i]; 
                    xz_[j*Nx+i] = psipZ_(end[0][j], end[1][j])/f_x_[i]; 
                     g_[j*Nx+i] = end[2][j],  f_[j*Nx + i] = -1./f_x_[i];
                }
                //////////////////////////////////////////////////
            }
            dg::blas1::axpby( 1., psi_x, -1., psi_old, psi_diff);
            eps = sqrt( dg::blas2::dot( psi_diff, w1d, psi_diff)/ dg::blas2::dot( psi_x, w1d, psi_x));
            //psi_1_numerical_ = psi_0 + dg::blas1::dot( f_x_, w1d);

            //eps = fabs( psi_1_numerical-psi_1); 
            //std::cout << "Effective absolute Psi error is "<<psi_1_numerical-psi_1<<" with "<<N<<" steps\n"; 
            std::cout << "Effective R  error is "<<eps<<" with "<<N<<" steps\n"; 
            //std::cout << "R 1                is "<<psi_1_numerical_<<"\n"; 
            N*=2;
        }
        //r_x1 = r_x0, z_x1 = z_x0; //periodic boundaries
        //now lift to 3D grid
        for( unsigned i=0; i<Nx; i++)
            f_x_[i] = -1./f_x_[i];
        for( unsigned k=1; k<this->Nz(); k++)
            for( unsigned i=0; i<Nx*Ny; i++)
            {
                f_[k*Nx*Ny+i] = f_[(k-1)*Nx*Ny+i];
                g_[k*Nx*Ny+i] = g_[(k-1)*Nx*Ny+i];
                r_[k*Nx*Ny+i] = r_[(k-1)*Nx*Ny+i];
                z_[k*Nx*Ny+i] = z_[(k-1)*Nx*Ny+i];
                yr_[k*Nx*Ny+i] = yr_[(k-1)*Nx*Ny+i];
                yz_[k*Nx*Ny+i] = yz_[(k-1)*Nx*Ny+i];
                xr_[k*Nx*Ny+i] = xr_[(k-1)*Nx*Ny+i];
                xz_[k*Nx*Ny+i] = xz_[(k-1)*Nx*Ny+i];
            }
    }
    //compute metric elements from xr, xz, yr, yz, r and z
    void construct_metric()
    {
        std::cout << "CONSTRUCTING METRIC\n";
        thrust::host_vector<double> tempxx( r_), tempxy(r_), tempyy(r_), tempvol(r_);
        unsigned Nx = this->n()*this->Nx(), Ny = this->n()*this->Ny();
        for( unsigned k=0; k<this->Nz(); k++)
            for( unsigned i=0; i<Ny; i++)
                for( unsigned j=0; j<Nx; j++)
                {
                    unsigned idx = k*Ny*Nx+i*Nx+j;
                    tempxx[idx] = (xr_[idx]*xr_[idx]+xz_[idx]*xz_[idx]);
                    tempxy[idx] = (yr_[idx]*xr_[idx]+yz_[idx]*xz_[idx]);
                    tempyy[idx] = (yr_[idx]*yr_[idx]+yz_[idx]*yz_[idx]);
                    //tempvol[idx] = r_[idx]/(1.0*f_[idx]*f_[idx] + 0.001*tempxx[idx]);
                    tempvol[idx] = r_[idx]/sqrt(tempxx[idx]*tempyy[idx]-tempxy[idx]*tempxy[idx]);
                    //tempvol[idx] = r_[idx]/fabs(f_[idx]*g_[idx]);
                }
        g_xx_=tempxx, g_xy_=tempxy, g_yy_=tempyy, vol_=tempvol;
        dg::blas1::pointwiseDivide( tempvol, r_, tempvol);
        vol2d_ = tempvol;
        thrust::host_vector<double> ones = dg::evaluate( dg::one, *this);
        dg::blas1::pointwiseDivide( ones, r_, tempxx);
        dg::blas1::pointwiseDivide( tempxx, r_, tempxx); //1/R^2
        g_pp_=tempxx;
    }
    thrust::host_vector<double> f_x_; //1d vector
    thrust::host_vector<double> f_, g_, r_, z_, xr_, xz_, yr_, yz_; //3d vector
    container g_xx_, g_xy_, g_yy_, g_pp_, vol_, vol2d_;
    thrust::host_vector<double> r_x0, r_x1, z_x0, z_x1; //boundary points in y
    double psi_1_numerical_;
};

/**
 * @brief A three-dimensional grid based on "almost-conformal" coordinates by Ribeiro and Scott 2010
 */
template< class container>
struct OrthogonalXGrid2d : public dg::GridX2d
{
    typedef dg::CurvilinearCylindricalTag metric_category;
    OrthogonalXGrid2d( const GeomParameters gp, double psi_0, double fx, double fy, double y, unsigned n, unsigned Nx, unsigned Ny, dg::bc bcx, dg::bc bcy): 
        dg::GridX2d( 0, 1,-fy*2.*M_PI/(1.-2.*fy), 2*M_PI+fy*2.*M_PI/(1.-2.*fy), fx, fy, n, Nx, Ny, bcx, bcy)
    {
        solovev::detail::FpsiX fpsi(gp);
        const double x0 = fpsi.find_x(psi_0);
        const double x1 = -fx/(1.-fx)*x0;
        //const double psi1 = fpsi.find_x(1);
        init_X_boundaries( x0,x1);
        OrthogonalXGrid3d<container> g( gp, psi_0, fx,fy, n,Nx,Ny,1,bcx,bcy);
        f_x_ = g.f_x();
        f_ = g.f(), g_ = g.g(), r_=g.r(), z_=g.z(), xr_=g.xr(), xz_=g.xz(), yr_=g.yr(), yz_=g.yz();
        g_xx_=g.g_xx(), g_xy_=g.g_xy(), g_yy_=g.g_yy();
        vol2d_=g.perpVol();
    }
    OrthogonalXGrid2d( const OrthogonalXGrid3d<container>& g):
        dg::GridX2d( g.x0(), g.x1(), g.y0(), g.y1(), g.fx(), g.fy(), g.n(), g.Nx(), g.Ny(), g.bcx(), g.bcy())
    {
        f_x_ = g.f_x();
        unsigned s = this->size();
        f_.resize(s), g_.resize(s), r_.resize( s), z_.resize(s), xr_.resize(s), xz_.resize(s), yr_.resize(s), yz_.resize(s);
        g_xx_.resize( s), g_xy_.resize(s), g_yy_.resize(s), vol2d_.resize(s);
        for( unsigned i=0; i<s; i++)
        {f_[i] = g.f()[i], g_[i] = g.g()[i], r_[i]=g.r()[i], z_[i]=g.z()[i], xr_[i]=g.xr()[i], xz_[i]=g.xz()[i], yr_[i]=g.yr()[i], yz_[i]=g.yz()[i];}
        thrust::copy( g.g_xx().begin(), g.g_xx().begin()+s, g_xx_.begin());
        thrust::copy( g.g_xy().begin(), g.g_xy().begin()+s, g_xy_.begin());
        thrust::copy( g.g_yy().begin(), g.g_yy().begin()+s, g_yy_.begin());
        thrust::copy( g.perpVol().begin(), g.perpVol().begin()+s, vol2d_.begin());
    }
    const thrust::host_vector<double>& f()const{return f_;}
    const thrust::host_vector<double>& g()const{return g_;}
    const thrust::host_vector<double>& r()const{return r_;}
    const thrust::host_vector<double>& z()const{return z_;}
    const thrust::host_vector<double>& xr()const{return xr_;}
    const thrust::host_vector<double>& yr()const{return yr_;}
    const thrust::host_vector<double>& xz()const{return xz_;}
    const thrust::host_vector<double>& yz()const{return yz_;}
    thrust::host_vector<double> x()const{
        dg::Grid1d<double> gx( x0(), x1(), n(), Nx());
        return dg::create::abscissas(gx);}
    const thrust::host_vector<double>& f_x()const{return f_x_;}
    const container& g_xx()const{return g_xx_;}
    const container& g_yy()const{return g_yy_;}
    const container& g_xy()const{return g_xy_;}
    const container& vol()const{return vol2d_;}
    const container& perpVol()const{return vol2d_;}
    private:
    thrust::host_vector<double> f_x_; //1d vector
    thrust::host_vector<double> f_, g_, r_, z_, xr_, xz_, yr_, yz_; //2d vector
    container g_xx_, g_xy_, g_yy_, vol2d_;
};

/**
 * @brief Integrates the equations for a field line and 1/B
 */ 
struct OrthogonalXField
{
    OrthogonalXField( solovev::GeomParameters gp,const dg::GridX2d& gXY, const thrust::host_vector<double>& g):
        gp_(gp),
        psipR_(gp), psipZ_(gp),
        ipol_(gp), invB_(gp), gXY_(gXY), g_(dg::create::forward_transform(g, gXY)) 
    { 
        HessianRZtau hessianRZtau(gp);
        R_X = gp.R_0-1.1*gp.triangularity*gp.a;
        Z_X = -1.1*gp.elongation*gp.a;
        thrust::host_vector<double> X(2,0), XN(X);
        X[0] = R_X, X[1] = Z_X;
        for( unsigned i=0; i<3; i++)
        {
            hessianRZtau.newton_iteration( X, XN);
            XN.swap(X);
        }
        R_X = X[0], Z_X = X[1];
    
    }

    /**
     * @brief \f[ \frac{d \hat{R} }{ d \varphi}  = \frac{\hat{R}}{\hat{I}} \frac{\partial\hat{\psi}_p}{\partial \hat{Z}}, \hspace {3 mm}
     \frac{d \hat{Z} }{ d \varphi}  =- \frac{\hat{R}}{\hat{I}} \frac{\partial \hat{\psi}_p}{\partial \hat{R}} , \hspace {3 mm}
     \frac{d \hat{l} }{ d \varphi}  =\frac{\hat{R}^2 \hat{B}}{\hat{I}  \hat{R}_0}  \f]
     */ 
    void operator()( const dg::HVec& y, dg::HVec& yp)
    {
        //x,y,s,R,Z
        double psipR = psipR_(y[3],y[4]), psipZ = psipZ_(y[3],y[4]), ipol = ipol_( y[3],y[4]);
        double xs = y[0],ys=y[1];
        if( y[4] > Z_X) //oberhalb vom X-Punkt
            gXY_.shift_topologic( y[0], M_PI, xs,ys);
        else 
        {
            if( y[1] > M_PI) //Startpunkt vermutlich in der rechten Hälfte
                gXY_.shift_topologic( y[0], gXY_.y1()-1e-10, xs,ys);
            else
                gXY_.shift_topologic( y[0], gXY_.y0()+1e-10, xs,ys);
        }
        if( !gXY_.contains(xs,ys))
        {
            if(y[0] > R_X) ys = gXY_.y1()-1e-10;
            else           ys = gXY_.y0()+1e-10;
        }
        double g = dg::interpolate( xs,  ys, g_, gXY_);
        yp[0] =  0;
        yp[1] =  y[3]*g*(psipR*psipR+psipZ*psipZ)/ipol;
        yp[2] =  y[3]*y[3]/invB_(y[3],y[4])/ipol/gp_.R_0; //ds/dphi =  R^2 B/I/R_0_hat
        yp[3] =  y[3]*psipZ/ipol;              //dR/dphi =  R/I Psip_Z
        yp[4] = -y[3]*psipR/ipol;             //dZ/dphi = -R/I Psip_R

    }
    /**
     * @brief \f[   \frac{1}{\hat{B}} = 
      \frac{\hat{R}}{\hat{R}_0}\frac{1}{ \sqrt{ \hat{I}^2  + \left(\frac{\partial \hat{\psi}_p }{ \partial \hat{R}}\right)^2
      + \left(\frac{\partial \hat{\psi}_p }{ \partial \hat{Z}}\right)^2}}  \f]
     */ 
    double operator()( double R, double Z) const { return invB_(R,Z); }
    /**
     * @brief == operator()(R,Z)
     */ 
    double operator()( double R, double Z, double phi) const { return invB_(R,Z,phi); }
    
    private:
    GeomParameters gp_;
    PsipR  psipR_;
    PsipZ  psipZ_;
    Ipol   ipol_;
    InvB   invB_;
    const dg::GridX2d gXY_;
    thrust::host_vector<double> g_;
    double R_X, Z_X;
   
};

}//namespace solovev

namespace dg{

/**
 * @brief This function pulls back a function defined in cartesian coordinates R,Z to the conformal coordinates x,y,\phi
 *
 * i.e. F(x,y) = f(R(x,y), Z(x,y))
 * @tparam BinaryOp The function object 
 * @param f The function defined on R,Z
 * @param g The grid
 *
 * @return A set of points representing F(x,y)
 */
template< class BinaryOp, class container>
thrust::host_vector<double> pullback( BinaryOp f, const solovev::OrthogonalXGrid2d<container>& g)
{
    thrust::host_vector<double> vec( g.size());
    for( unsigned i=0; i<g.size(); i++)
        vec[i] = f( g.r()[i], g.z()[i]);
    return vec;
}

///@cond
template<class container>
thrust::host_vector<double> pullback( double(f)(double,double), const solovev::OrthogonalXGrid2d<container>& g)
{
    return pullback<double(double,double),container>( f, g);
}
///@endcond

/**
 * @brief This function pulls back a function defined in cylindrical coordinates R,Z,\phi to the conformal coordinates x,y,\phi
 *
 * i.e. F(x,y,\phi) = f(R(x,y), Z(x,y), \phi)
 * @tparam TernaryOp The function object 
 * @param f The function defined on R,Z,\phi
 * @param g The grid
 *
 * @return A set of points representing F(x,y,\phi)
 */
template< class TernaryOp, class container>
thrust::host_vector<double> pullback( TernaryOp f, const solovev::OrthogonalXGrid3d<container>& g)
{
    thrust::host_vector<double> vec( g.size());
    unsigned size2d = g.n()*g.n()*g.Nx()*g.Ny();
    Grid1d<double> gz( g.z0(), g.z1(), 1, g.Nz());
    thrust::host_vector<double> absz = create::abscissas( gz);
    for( unsigned k=0; k<g.Nz(); k++)
        for( unsigned i=0; i<size2d; i++)
            vec[k*size2d+i] = f( g.r()[k*size2d+i], g.z()[k*size2d+i], absz[k]);
    return vec;
}
///@cond
template<class container>
thrust::host_vector<double> pullback( double(f)(double,double,double), const solovev::OrthogonalXGrid3d<container>& g)
{
    return pullback<double(double,double,double),container>( f, g);
}
///@endcond

}//namespace dg

