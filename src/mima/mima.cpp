#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#ifdef WITH_MPI
#include <mpi.h>
#endif // WITH_MPI

#ifdef WITH_GLFW
#include "draw/host_window.h"
#endif // WITH_GLFW

#include "dg/algorithm.h"
#include "dg/file/file.h"

#include "mima.h"
#include "diag.h"

int main( int argc, char* argv[])
{
#ifdef WITH_MPI
    dg::mpi_init( argc, argv);
    MPI_Comm comm;
    dg::mpi_init2d( dg::DIR, dg::PER, comm, std::cin, true);
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
#endif //WITH_MPI
    dg::file::WrappedJsonValue js( dg::file::error::is_throw);
    mima::Parameters p;
    try{
        std::string inputfile = "input/default.json";
        if( argc != 1) inputfile = argv[1];
        dg::file::file2Json( inputfile.c_str(), js.asJson(),
                dg::file::comments::are_discarded, dg::file::error::is_throw);
        p = { js};
    } catch( std::exception& e) {
        DG_RANK0 std::cerr << "ERROR in input file "<<argv[1]<<std::endl;
        DG_RANK0 std::cerr << e.what()<<std::endl;
        dg::abort_program();
    }
    DG_RANK0 std::cout << js.asJson() << std::endl;
    DG_RANK0 p.display(std::cout);

    //Construct grid
    dg::x::CartesianGrid2d grid( 0, p.lx, 0., p.ly, p.n, p.Nx, p.Ny,
            p.bcx, p.bcy
        #ifdef WITH_MPI
        , comm
        #endif //WITH_MPI
        );
    //create RHS
    mima::Explicit<dg::x::CartesianGrid2d, dg::x::DMatrix, dg::x::DVec>
        rhs( grid, p);
    //////////////////create initial vector///////////////////////////////////////
    dg::Gaussian g( p.posX*p.lx, p.posY*p.ly, p.sigma, p.sigma, p.amp); //gaussian initial condition
    dg::Gaussian g1( 0.45*p.lx, p.posY*p.ly, p.sigma, p.sigma, p.amp);
    dg::Gaussian g2( 0.55*p.lx, p.posY*p.ly, p.sigma, p.sigma, p.amp);
    dg::SinY siny(p.amp, 0., p.posY*2.*M_PI/p.ly); // siny initial condition where posY denotes ky
    dg::Vortex vortex(p.posX*p.lx, p.posY*p.ly, p.state, p.sigma, p.amp); //
    dg::Lamb lamb(p.posX*p.lx, p.posY*p.ly, p.sigma, p.amp);
    dg::BathRZ bath(p.N_kR, p.N_kZ, p.R_min, p.Z_min, p.bath_gamma, p.L_E, p.bath_amp);

    dg::x::DVec y0;
    dg::x::DVec firstgauss;
    dg::x::DVec secondgauss;
    
    dg::Helmholtz<dg::x::CartesianGrid2d, dg::x::DMatrix, dg::x::DVec> m_helmholtz( -1., {grid, dg::centered});
    dg::Helmholtz<dg::x::CartesianGrid2d, dg::x::DMatrix, dg::x::DVec> m_invgamma0( -p.taui, {grid, dg::centered});
    dg::Helmholtz<dg::x::CartesianGrid2d, dg::x::DMatrix, dg::x::DVec> m_invgamma1( -p.taui/2., {grid, dg::centered});
    dg::Helmholtz<dg::x::CartesianGrid2d, dg::x::DMatrix, dg::x::DVec> m_invgamma2( -p.taui-1., {grid, dg::centered}); //invgamma2 = 1-(taui+1)lap */
    if (p.init_cond == "siny"){ 
        y0 = dg::evaluate(siny, grid); // n_i = sin(y)
        }
    else if (p.init_cond == "vortex"){
        y0 = dg::evaluate(vortex, grid);
    }
    else if (p.init_cond == "lamb"){
        y0 = dg::evaluate(lamb, grid);
    }
    else if (p.init_cond == "bath"){
        y0 = dg::evaluate(bath, grid);
    }
    else if (p.init_cond == "doublegauss"){
        y0 = dg::evaluate(g, grid);
        firstgauss = dg::evaluate(g1,grid);
        secondgauss = dg::evaluate(g2,grid);
        dg::blas1::axpby(1.,firstgauss, 1., secondgauss, y0);
    }
    else{
        y0 = dg::evaluate(g, grid);
    }
    auto temp1 = y0;
    
    if (p.model == "FLRapprox"){
        dg::blas2::symv(m_invgamma2, y0, temp1);
        dg::PCG<dg::x::DVec> m_pcg( y0, grid.size());
        m_pcg.solve(m_invgamma1, y0, temp1, m_invgamma0.precond(), m_invgamma0.weights(), p.eps_gamma);
    }
    else{
        dg::blas2::symv(m_invgamma1, temp1, y0);
        dg::blas2::symv(m_invgamma2, y0, temp1);
        dg::PCG<dg::x::DVec> m_pcg( y0, grid.size());
        m_pcg.solve(m_invgamma0, y0, temp1, m_invgamma0.precond(), m_invgamma0.weights(), p.eps_gamma);
    }
    
    

    
    
/*
    if( p.model == "gravity_local" || p.model == "gravity_global" ||
            p.model == "drift_global"){
        y0[1] = dg::evaluate( dg::zero, grid);
    }*/
    //////////////////////////////////////////////////////////////////////
    // Construct timestepper
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
    using Vec = dg::x::DVec;
    dg::Adaptive< dg::ERKStep< Vec>> adapt(tableau, y0);
    dg::AdaptiveTimeloop<Vec> timeloop( adapt, rhs,
                        dg::pid_control, dg::l2norm, rtol, atol);

    ////////////////////////////////////////////////////////////////////

    mima::Variables var = { rhs, grid, p};
    // trigger first computation of potential
    {
        DG_RANK0 std::cout<< "First potential\n";
        auto temp = y0;
        rhs( 0., y0, temp);
        DG_RANK0 std::cout<< "Done\n";
    }
    std::string output = js[ "output"]["type"].asString("glfw");
    dg::Timer t;
    t.tic();
#ifdef WITH_GLFW
    if( "glfw" == output)
    {
        double dt = 1e-5;
        /////////glfw initialisation ////////////////////////////////////////////
        dg::file::WrappedJsonValue ws;
        dg::file::file2Json( "window_params.json", ws.asJson(),
                dg::file::comments::are_discarded);
        GLFWwindow* w = draw::glfwInitAndCreateWindow( ws["width"].asDouble(),
                ws["height"].asDouble(), "");
        draw::RenderHostData render(ws["rows"].asDouble(), ws["cols"].asDouble());
        //create visualisation vectors
        dg::DVec visual( grid.size()), temp(visual);
        dg::HVec hvisual( grid.size());
        //transform vector to an equidistant grid
        std::stringstream title;
        draw::ColorMapRedBlueExtMinMax colors( -1.0, 1.0);
        dg::IDMatrix equidistant = dg::create::backscatter( grid );
        // the things to plot:
        std::map<std::string, const dg::DVec* > v2d;
        v2d["ne / "] = &rhs.phi();
        v2d["Vor / "] = &rhs.chi();
        unsigned itstp = js["output"]["itstp"].asUInt();
        unsigned step = 0;

        while ( !glfwWindowShouldClose( w ))
        {
            for( auto pair : v2d)
            {
                if( pair.first == "Vor / " )
                {
                    //dg::blas1::axpby(1., rhs.phi(), -1., y0, 0.,temp);
                    //dg::blas2::gemv( rhs.laplacianM(), *pair.second, temp);
                    dg::blas2::gemv( equidistant, rhs.chi(), visual);
                }
                else
                    dg::blas2::gemv( equidistant, *pair.second, visual);
                dg::assign( visual, hvisual);
                colors.scalemax() = dg::blas1::reduce(
                    hvisual, 0., dg::AbsMax<double>() );
                colors.scalemin() = -colors.scalemax();
                title << std::setprecision(2) << std::scientific;
                title <<pair.first << colors.scalemax()<<"   ";
                render.renderQuad( hvisual, grid.n()*grid.Nx(),
                        grid.n()*grid.Ny(), colors);
            }
            title << std::fixed;
            title << " &&   time = "<<time;
            glfwSetWindowTitle(w,title.str().c_str());
            title.str("");
            glfwPollEvents();
            glfwSwapBuffers( w);

            //step
            dg::Timer ti;
            ti.tic();
            std::cout << "\n\t dt "<<dt<<"\n";
            for( unsigned i=0; i<itstp; i++)
            {
                try{
                    adapt.step( rhs, time, y0, time, y0, dt, dg::pid_control,
                            dg::l2norm, rtol, atol);
                }
                catch( std::exception& fail) {
                    std::cerr << "ERROR in Timestepper\n";
                    std::cerr << fail.what() << std::endl;
                    std::cerr << "Does Simulation respect CFL condition?\n";
                    glfwSetWindowShouldClose( w, GL_TRUE);
                    break;
                }
                step++;
            }
            ti.toc();
            std::cout << "\n\t Step "<<step;
            std::cout << "\n\t Time "<<time<<" dt "<<dt;
            std::cout << "\n\t Average time for one step: "<<ti.diff()/(double)itstp<<"s\n\n";
        }
        glfwTerminate();
    }
#endif //WITH_GLFW
    if( "netcdf" == output)
    {
        std::string outputfile;
        if( argc == 1 || argc == 2)
            outputfile = "mima.nc";
        else
            outputfile = argv[2];
        // Create netcdf file
        dg::file::NC_Error_Handle err;
        int ncid=-1;
        try{
            DG_RANK0 err = nc_create(outputfile.c_str(), NC_NETCDF4|NC_CLOBBER, &ncid);
        }catch( std::exception& e)
        {
            DG_RANK0 std::cerr << "ERROR creating file "<<argv[1]<<std::endl;
            DG_RANK0 std::cerr << e.what() << std::endl;
            dg::abort_program();
        }
        std::map<std::string, std::string> att;
        att["title"] = "Output file of feltor/src/mima/mima.cpp";
        att["Conventions"] = "CF-1.8";
        ///Get local time and begin file history
        auto ttt = std::time(nullptr);

        std::ostringstream oss;
        ///time string  + program-name + args
        oss << std::put_time(std::localtime(&ttt), "%F %T %Z");
        for( int i=0; i<argc; i++) oss << " "<<argv[i];
        att["history"] = oss.str();
        att["comment"] = "Find more info in feltor/src/mima/mima.tex";
        att["source"] = "FELTOR";
        att["git-hash"] = GIT_HASH;
        att["git-branch"] = GIT_BRANCH;
        att["compile-time"] = COMPILE_TIME;
        att["references"] = "https://github.com/feltor-dev/feltor";
        // Here we put the inputfile as a string without comments so that it can be read later by another parser
        att["inputfile"] = js.asJson().toStyledString();
        for( auto pair : att)
            DG_RANK0 err = nc_put_att_text( ncid, NC_GLOBAL,
                pair.first.data(), pair.second.size(), pair.second.data());

        unsigned n_out     = js[ "output"]["n"].asUInt( 3);
        unsigned Nx_out    = js[ "output"]["Nx"].asUInt( 48);
        unsigned Ny_out    = js[ "output"]["Ny"].asUInt( 48);

        dg::x::CartesianGrid2d grid_out( 0., p.lx, 0., p.ly,
                    n_out, Nx_out, Ny_out, p.bcx, p.bcy
                    #ifdef WITH_MPI
                    , comm
                    #endif //WITH_MPI
                    );
        dg::x::IHMatrix projection = dg::create::interpolation( grid_out, grid);
        int dim_ids[3], tvarID;
        // the dimensions are the ones of grid_out!
        err = dg::file::define_dimensions( ncid, dim_ids, &tvarID, grid_out,
                        {"time", "y", "x"});

        std::map<std::string, int> id3d;
        for( auto& record : mima::diagnostics2d_list.at( p.model))
        {
            std::string name = record.name;
            std::string long_name = record.long_name;
            id3d[name] = 0;
            DG_RANK0 err = nc_def_var( ncid, name.data(), NC_DOUBLE, 3, dim_ids,
                    &id3d.at(name));
            DG_RANK0 err = nc_put_att_text( ncid, id3d.at(name), "long_name",
                    long_name.size(), long_name.data());
        }
        dg::x::HVec resultH = dg::evaluate( dg::zero, grid);
        dg::x::HVec transferH = dg::evaluate( dg::zero, grid_out);
        dg::x::DVec resultD = transferH; // transfer to device
        for( auto& record : mima::diagnostics2d_static_list)
        {
            std::string name = record.name;
            std::string long_name = record.long_name;
            int staticID = 0;
            DG_RANK0 err = nc_def_var( ncid, name.data(), NC_DOUBLE, 2, &dim_ids[1],
                &staticID);
            DG_RANK0 err = nc_put_att_text( ncid, staticID, "long_name",
                                           long_name.size(), long_name.data());
            record.function( resultD, var);
            dg::assign( resultD, resultH);
            dg::blas2::gemv( projection, resultH, transferH);
            dg::file::put_var_double( ncid, staticID, grid_out, transferH);
        }
        dg::x::DVec volume = dg::create::volume( grid);
        size_t start = {0};
        size_t count = {1};
        for( auto& record : mima::diagnostics2d_list.at( p.model))
        {
            record.function( resultD, var);
            dg::assign( resultD, resultH);
            dg::blas2::gemv( projection, resultH, transferH);
            // note that all processes call this function (for MPI)
            dg::file::put_vara_double( ncid, id3d.at(record.name), start,
                    grid_out, transferH);
        }
        DG_RANK0 err = nc_put_vara_double( ncid, tvarID, &start, &count, &time);
        DG_RANK0 err = nc_close( ncid);
        double Tend = js["output"].get("tend", 1.0).asDouble();
        unsigned maxout = js["output"].get("maxout", 10).asUInt();
        double deltaT = Tend/(double)maxout;
        bool abort = false;
        unsigned ncalls = rhs.ncalls();
        for( unsigned u=1; u<=maxout; u++)
        {
            dg::Timer ti;
            ti.tic();
            try{
                timeloop.integrate( time, y0, u*deltaT, y0,
                                  u < maxout ? dg::to::at_least : dg::to::exact);
            }catch ( std::exception& fail)
            {
                DG_RANK0 std::cerr << "ERROR in Timestepper\n";
                DG_RANK0 std::cerr << fail.what() << std::endl;
                DG_RANK0 std::cerr << "Writing last output and exit ..."<<std::endl;
                abort = true;
            }
            unsigned delta_ncalls = rhs.ncalls() - ncalls;
            ncalls = rhs.ncalls();
            ti.toc();
            DG_RANK0 std::cout << "\n\t Time "<<time <<" of "<<Tend <<" with current timestep "<<timeloop.get_dt();
            DG_RANK0 std::cout << "\n\t # of rhs calls since last output "<<delta_ncalls;
            DG_RANK0 std::cout << "\n\t Average time for one step: "<<ti.diff()/(double)delta_ncalls<<"s\n\n"<<std::flush;
            start = u;
            DG_RANK0 err = nc_open(outputfile.c_str(), NC_WRITE, &ncid);
            // First write the time variable
            DG_RANK0 err = nc_put_vara_double( ncid, tvarID, &start, &count, &time);
            for( auto& record : mima::diagnostics2d_list.at(p.model))
            {
                record.function( resultD, var);
                dg::assign( resultD, resultH);
                dg::blas2::gemv( projection, resultH, transferH);
                dg::file::put_vara_double( ncid, id3d.at(record.name),
                                          start, grid_out, transferH);
            }
            DG_RANK0 err = nc_close( ncid);
            if( abort) break;
        }
    }
    if( !("netcdf" == output)
#ifdef WITH_GLFW
            && !("glfw" == output)
#endif
            )
    {
        DG_RANK0 std::cerr <<"Error: Wrong value for output type `"<<output<<"` Must be netcdf!";
#ifdef WITH_GLFW
        DG_RANK0 std::cerr <<" Or glfw!\n";
#endif
        DG_RANK0 std::cerr <<" Exit now!\n";
        dg::abort_program();
    }
    ////////////////////////////////////////////////////////////////////
    t.toc();
    unsigned hour = (unsigned)floor(t.diff()/3600);
    unsigned minute = (unsigned)floor( (t.diff() - hour*3600)/60);
    double second = t.diff() - hour*3600 - minute*60;
    DG_RANK0 std::cout << std::fixed << std::setprecision(2) <<std::setfill('0');
    DG_RANK0 std::cout <<"Computation Time \t"<<hour<<":"<<std::setw(2)<<minute<<":"<<second<<" for "<<rhs.ncalls()<<"\n";
    DG_RANK0 std::cout <<"which is         \t"<<t.diff()/rhs.ncalls()<<"s/rhs call\n";

#ifdef WITH_MPI
    MPI_Finalize();
#endif //WITH_MPI

    return 0;

}
