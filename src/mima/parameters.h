#pragma once
#include <string>
#include "dg/algorithm.h"
#include "dg/file/json_utilities.h"
#include "json/json.h"

namespace mima{

struct Parameters
{
    unsigned n, Nx, Ny;
    double lx, ly;
    dg::bc bcx, bcy;

    double eps_gamma;
    enum dg::direction diff_dir;

    double amp, sigma, posX, posY;

    std::string model;
    double Ln;
    Parameters() = default;

    Parameters( const dg::file::WrappedJsonValue& js) {
        n  = js["grid"]["n"].asUInt();
        Nx = js["grid"]["Nx"].asUInt();
        Ny = js["grid"]["Ny"].asUInt();
        lx = js["grid"]["lx"].asDouble();
        ly = js["grid"]["ly"].asDouble();

        eps_gamma = js["elliptic"]["eps_gamma"].asDouble();
        
        diff_dir = dg::centered;
        sigma = (double(std::min(lx,ly))/32.);
        amp = js["init"]["amp"].asDouble();
        //amp = (double(std::min(lx,ly))/320.);
        posX = js["init"]["posX"].asDouble();
        posY = js["init"]["posY"].asDouble();
        bcx = dg::str2bc(js["bc"][0].asString());
        bcy = dg::str2bc(js["bc"][1].asString());
        model = js["model"].get("type", "standardCHM").asString();
        
        if( "standardCHM" == model)
        {
            Ln = js["model"]["Ln"].asDouble();
        }
        else if( "global" == model)
        {
            Ln = js["model"]["Ln"].asDouble();
        }
        else
            throw dg::Error( dg::Message(_ping_) << "Model : type `"<<model<<"` not recognized!\n");
    }

    void display( std::ostream& os = std::cout ) const
    {
        os << "Physical parameters are: \n"
            <<"    Background gradient length:     = "<<Ln<<"\n";
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
        os  <<"Blob parameters are: \n"
            << "    width:        "<<sigma<<"\n"
            << "    amplitude:    "<<amp<<"\n"
            << "    posX:         "<<posX<<"\n"
            << "    posY:         "<<posY<<"\n";
        os  <<"Stopping for Gamma CG:   "<<eps_gamma<<std::endl;
    }
};
}//namespace mima
