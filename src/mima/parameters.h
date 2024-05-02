#pragma once
#include <string>
#include "dg/algorithm.h"
#include "dg/file/json_utilities.h"
#include "json/json.h"

namespace mima{

struct Parameters
{
    unsigned n, Nx, Ny, state;
    double lx, ly;
    dg::bc bcx, bcy;

    double eps_gamma;
    enum dg::direction diff_dir;

    double amp, sigma, posX, posY, R;

    unsigned N_kR, N_kZ;
    double R_min, Z_min, bath_gamma, L_E, bath_amp;

    std::string model;
    std::string init_cond;
    double Ln, taui;
    Parameters() = default;

    Parameters( const dg::file::WrappedJsonValue& js) {
        n  = js["grid"]["n"].asUInt();
        Nx = js["grid"]["Nx"].asUInt();
        Ny = js["grid"]["Ny"].asUInt();
        lx = js["grid"]["lx"].asDouble();
        ly = js["grid"]["ly"].asDouble();

        state = js["init"]["state"].asUInt();

        eps_gamma = js["elliptic"]["eps_gamma"].asDouble();
        
        diff_dir = dg::centered;
        sigma = (double(std::min(lx,ly))/45.);
        amp = js["init"]["amp"].asDouble();
        //amp = (double(std::min(lx,ly))/320.);
        posX = js["init"]["posX"].asDouble();
        posY = js["init"]["posY"].asDouble();
        bcx = dg::str2bc(js["bc"][0].asString());
        bcy = dg::str2bc(js["bc"][1].asString());
        model = js["model"].get("type", "standardCHM").asString();
        taui = 0;
        init_cond = js["model"]["init_cond"].asString();


        N_kR = js["bath"]["N_kR"].asUInt();
        N_kZ = js["bath"]["N_kZ"].asUInt();
        R_min = js["bath"]["R_min"].asDouble();
        Z_min = js["bath"]["Z_min"].asDouble();
        bath_gamma = js["bath"]["gamma"].asDouble();
        L_E = js["bath"]["L_E"].asDouble();
        bath_amp = js["bath"]["bath_amp"].asDouble();


        if( "standardCHM" == model)
        {
            Ln = js["model"]["Ln"].asDouble();
        }
        else if( "boussinesq" == model)
        {
            Ln = js["model"]["Ln"].asDouble();
        }
        else if( "boussinesq2" == model)
        {
            Ln = js["model"]["Ln"].asDouble();
        }
        else if( "FLR" == model)
        {
            Ln = js["model"]["Ln"].asDouble();
            taui = js["model"]["taui"].asDouble();
        }
        else
            throw dg::Error( dg::Message(_ping_) << "Model : type `"<<model<<"` not recognized!\n");
    }

    void display( std::ostream& os = std::cout ) const
    {
        os << "Physical parameters are: \n"
            <<"    Background gradient length     = "<<Ln<<"\n"
            <<"    Ion temperature taui           = "<<taui<<"\n";
        os << "Equation parameters are: \n"
            <<"    "<<model<<"\n";
        os << "Boundary parameters are: \n"
            <<"    lx = "<<lx<<"\n"
            <<"    ly = "<<ly<<"\n";
        os << "Boundary conditions in x are: \n"
            <<"    "<<bc2str(bcx)<<"\n";  //Curious! dg:: is not needed due to ADL!
        os << "Boundary conditions in y are: \n"
            <<"    "<<bc2str(bcy)<<"\n";
        os << "Algorithmic parameters are: \n"
            <<"    n  = "<<n<<"\n"
            <<"    Nx = "<<Nx<<"\n"
            <<"    Ny = "<<Ny<<"\n";
        os  <<"Initial conditions are: \n"
            << "    Initialized function:            "<<init_cond<<"\n"
            << "    Blob width / Vortex radius:        "<<sigma<<"\n"
            << "    Blob / Sine amplitude / u_dipole:    "<<amp<<"\n"
            << "    posX:         "<<posX<<"\n"
            << "    posY / Sine_k_y:         "<<posY<<"\n"
            << "    vortex state:      "<<state<<"\n";
        os  <<"Stopping for Gamma CG:   "<<eps_gamma<<std::endl;
    }
};
}//namespace mima
