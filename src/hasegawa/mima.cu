#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>

#include "draw/host_window.h"

//#include "draw/device_window.cuh"


#include "mima.cuh"
#include "../toefl/parameters.h"
#include "dg/file/json_utilities.h"

int main( int argc, char* argv[])
{
    ////Parameter initialisation ////////////////////////////////////////////
    std::stringstream title;
    Json::Value js;
    if( argc == 1)
        dg::file::file2Json( "input/default.json", js, dg::file::comments::are_discarded, dg::file::error::is_throw);
    else if( argc == 2)
        dg::file::file2Json( argv[1], js, dg::file::comments::are_discarded, dg::file::error::is_throw);
    else
    {
        std::cerr << "ERROR: Too many arguments!\nUsage: "<< argv[0]<<" [filename]\n";
        return -1;
    }
    std::cout<<"Point0\n";
    toefl::Parameters p( js);std::cout<<"Point1\n";
    p.display( std::cout); std::cout<<"Point2\n";
    /////////glfw initialisation ////////////////////////////////////////////
    dg::file::file2Json( "window_params.json", js, dg::file::comments::are_discarded, dg::file::error::is_throw); std::cout<<"Point2\n";
    GLFWwindow* w = draw::glfwInitAndCreateWindow( js["width"].asDouble(), js["height"].asDouble(), ""); std::cout<<"Point3\n";
    draw::RenderHostData render(js["rows"].asDouble(), js["cols"].asDouble());
    /////////////////////////////////////////////////////////////////////////
    std::cout<<"Point10\n";
    dg::x::CartesianGrid2d grid( 0, p.lx, 0, p.ly, p.n, p.Nx, p.Ny, p.bcx, p.bcy);
    //create RHS
    
    bool mhw = ( p.model == "fullF");
    mima::Mima< dg::x::CartesianGrid2d, dg::x::DMatrix, dg::x::DVec > mima( grid, p.kappa, p.tau, p.eps_pol[0], p.nu, mhw);
    dg::x::DVec one( grid.size(), 1.);
    //create initial vector
    dg::Gaussian gaussian( p.posX*grid.lx(), p.posY*grid.ly(), p.sigma, p.sigma, p.amp); //gaussian width is in absolute values
    dg::Vortex vortex( p.posX*grid.lx(), p.posY*grid.ly(), 0, p.sigma, p.amp);

//     dg::DVec phi = dg::evaluate( vortex, grid), omega( phi), y0(phi), y1(phi);
    dg::x::DVec phi = dg::evaluate( gaussian, grid), omega( phi), y0(phi), y1(phi);
    dg::Elliptic<dg::x::CartesianGrid2d, dg::x::DMatrix, dg::x::DVec> laplaceM( grid,  dg::centered);
    dg::blas2::gemv( laplaceM, phi, omega);
    dg::blas1::axpby( 1., phi, 1., omega, y0);

    dg::DVec w2d( dg::create::weights( grid));
    if( p.bcx == dg::PER && p.bcy == dg::PER)
    {
        double meanMass = dg::blas2::dot( y0, w2d, one)/(double)(p.lx*p.ly);
        std::cout << "Mean Mass is "<<meanMass<<"\n";
        dg::blas1::axpby( -meanMass, one, 1., y0);
    }
    std::string tableau;
    double rtol, atol, time = 0.;
    try{
        rtol = js["timestepper"].get("rtol", 1e-5).asDouble();
        atol = js["timestepper"].get("atol", 1e-5).asDouble();
        tableau = js[ "timestepper"].get( "tableau",
                "Bogacki-Shampine-4-2-3").asString();
    }catch ( std::exception& error){
        DG_RANK0 std::cerr << "Error in input file " << argv[1]<< std::endl;
        DG_RANK0 std::cerr << error.what() << std::endl;
        dg::abort_program();
    }
    DG_RANK0 std::cout<< "Construct timeloop ...\n";
    dg::Adaptive< dg::ERKStep< dg::DVec>> adapt(tableau, y0);

    dg::DVec dvisual( grid.size(), 0.);
    dg::HVec hvisual( grid.size(), 0.), visual(hvisual);
    dg::IHMatrix equi = dg::create::backscatter( grid);
    draw::ColorMapRedBlueExt colors( 1.);
    //create timer
    dg::Timer t;
    double dt = 1e-5;
    unsigned itstp = js["output"]["itstp"].asUInt();
    std::cout << "Begin computation \n";
    std::cout << std::scientific << std::setprecision( 2);
    unsigned step = 0;
    while ( !glfwWindowShouldClose( w ))
    {
        if( p.bcx == dg::PER && p.bcy == dg::PER)
        {
            double meanMass = dg::blas2::dot( y0, w2d, one)/(double)(p.lx*p.ly);
            std::cout << "Mean Mass is "<<meanMass<<"\n";
        }
        //transform field to an equidistant grid
        dvisual = mima.potential();

        dg::assign( dvisual, hvisual);
        dg::blas2::gemv( equi, hvisual, visual);
        //compute the color scale
        colors.scale() =  (float)thrust::reduce( visual.begin(), visual.end(), 0., dg::AbsMax<double>() );
        //draw ions
        title << std::setprecision(2) << std::scientific;
        title <<"ne / "<<colors.scale()<<"\t";
        render.renderQuad( visual, grid.n()*grid.Nx(), grid.n()*grid.Ny(), colors);

        //transform phi
        dg::blas2::gemv( laplaceM, mima.potential(), y1);
        dg::assign( y1, hvisual);
        dg::blas2::gemv( equi, hvisual, visual);
        //compute the color scale
        colors.scale() =  (float)thrust::reduce( visual.begin(), visual.end(), 0., dg::AbsMax<double>() );
        //draw phi and swap buffers
        title <<"omega / "<<colors.scale()<<"\t";
        title << std::fixed;
        title << " &&   time = "<<time;
        render.renderQuad( visual, grid.n()*grid.Nx(), grid.n()*grid.Ny(), colors);
        glfwSetWindowTitle(w,title.str().c_str());
        title.str("");
        glfwPollEvents();
        glfwSwapBuffers( w);

        //step 
#ifdef DG_BENCHMARK
        t.tic();
#endif//DG_BENCHMARK
        for( unsigned i=0; i<itstp; i++)
        {
            step++;
            if( p.bcx == dg::PER && p.bcy == dg::PER)
            {
                double meanMass = dg::blas2::dot( y0, w2d, one)/(double)(p.lx*p.ly);
                dg::blas1::axpby( -meanMass, one, 1., y0);
                meanMass = dg::blas2::dot( y0, w2d, one)/(double)(p.lx*p.ly);
                dg::blas1::axpby( -meanMass, one, 1., y0);
            }

            try{
                adapt.step( mima, time, y0, time, y0, dt, dg::pid_control,
                        dg::l2norm, rtol, atol);
            }
            catch( std::exception& fail) {
                std::cerr << "ERROR in Timestepper\n";
                std::cerr << fail.what() << std::endl;
                std::cerr << "Does Simulation respect CFL condition?\n";
                glfwSetWindowShouldClose( w, GL_TRUE);
                break;
            }
        }
#ifdef DG_BENCHMARK
        t.toc();
        std::cout << "\n\t Step "<<step;
        std::cout << "\n\t Average time for one step: "<<t.diff()/(double)itstp<<"s\n\n";
#endif//DG_BENCHMARK
    }
    glfwTerminate();
    ////////////////////////////////////////////////////////////////////

    return 0;

}
