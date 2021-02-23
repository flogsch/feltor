#pragma once
#include <string>
#include <vector>
#ifdef JSONCPP_VERSION_STRING
#include <dg/file/json_utilities.h>
#endif
/*!@file
 *
 * Geometry parameters
 */
namespace dg
{
namespace geo
{
namespace solovev
{
/*! @class hide_solovev_json
 * @code
// Solovev (and Taylor) geometry parameters
{
    "equilibrium": "solovev",
    // Note that for the taylor field you need to include boost before the geometries header!
    // "equilibrium" : "taylor",
    "A": 0,
    "R_0": 213.36,
    "PP": 1,
    "PI": 1,
    "c":[
        0.072597888572520090,
        -0.14926096478076946,
        // ... 12 coefficients in total
    ],
    "description" : "standardX",
    "inverseaspectratio": 0.3211009174311926,
    "triangularity": 0.3,
    "elongation": 1.44
}
@endcode
*/
/**
 * @brief Constructs and display geometric parameters for the solovev and taylor fields
 * @ingroup solovev
 * @note include \c json/json.h before \c geometries.h in order to activate json functionality
 */
struct Parameters
{
    double A, //!< A coefficient
           R_0, //!< major tokamak radius
           pp, //!< prefactor for Psi_p
           pi, //!< prefactor for current I
           a,  //!<  little tokamak radius
           elongation, //!< elongation of the magnetic surfaces
           triangularity; //!< triangularity of the magnetic surfaces
    std::vector<double> c;  //!< 12 coefficients for the solovev equilibrium;
    std::string description;
#ifdef JSONCPP_VERSION_STRING
    /**
     * @brief Construct from Json dataset
     * @copydoc hide_solovev_json
     * @sa \c dg::geo::description to see valid values for the %description field
     * @note the \c dg::geo::taylor field is chosen by setting "taylor" in the equilibrium field (but also note that you need to include boost for the taylor field)
     * @param js valid Json object (see code above to see the valid key : value pairs)
     * @param mode determine what happens when a key is missing
     * @note the default values in brackets are taken if the variables are not found in the input file
     * @attention This Constructor is only defined if \c json/json.h is included before \c dg/geometries/geometries.h
     */
    Parameters( const Json::Value& js, dg::file::error mode = dg::file::error::is_silent) {
        A  = dg::file::get( mode, js, "A", 0).asDouble();
        pp  = dg::file::get( mode, js, "PP", 1).asDouble();
        pi  = dg::file::get( mode, js, "PI", 1).asDouble();
        c.resize(12);
        for (unsigned i=0;i<12;i++)
            c[i] = dg::file::get_idx( mode, js, "c", i, 0.).asDouble();

        R_0  = dg::file::get( mode, js, "R_0", 0.).asDouble();
        a  = R_0*dg::file::get( mode, js, "inverseaspectratio", 0.).asDouble();
        elongation=dg::file::get( mode, js, "elongation", 1.).asDouble();
        triangularity=dg::file::get( mode, js, "triangularity", 0.).asDouble();
        try{
            description = dg::file::get( dg::file::error::is_throw, js, "description", "standardX").asString();
        } catch ( std::exception& err)
        {
            if( isToroidal())
                description = "none";
            else if( !hasXpoint())
                description = "standardO";
            else
                description = "standardX";
        }
    }
    /**
     * @brief Put values into a json string
     *
     * @return Json value
     * @attention This member is only defined if \c json/json.h is included before \c dg/geometries/geometries.h
     */
    Json::Value dump( ) const
    {
        Json::Value js;
        js["A"] = A;
        js["PP"] = pp;
        js["PI"] = pi;
        for (unsigned i=0;i<12;i++) js["c"][i] = c[i];
        js["R_0"] = R_0;
        js["inverseaspectratio"] = a/R_0;
        js["elongation"] = elongation;
        js["triangularity"] = triangularity;
        js[ "equilibrium"] = "solovev";
        js[ "description"] = description;
        return js;
    }
#endif // JSONCPP_VERSION_STRING
    /**
    * @brief True if any coefficient \c c_i!=0 with \c 7<=i<12
    *
    * The Xpoint is situated close to
     <tt> R_X = R_0-1.1*triangularity*a</tt>
     <tt> Z_X = -1.1*elongation*a</tt>
    *
    * @return \c true if Psip has an Xpoint, \c false else
    */
    bool hasXpoint( ) const{
        bool Xpoint = false;
        for( int i=7; i<12; i++)
            if( fabs(c[i]) >= 1e-10)
                Xpoint = true;
        return Xpoint;
    }
    /**
    * @brief True if \c pp==0
    *
    * @return \c true if the flux function is a constant
    */
    bool isToroidal() const{
        if( pp == 0)
            return true;
        return false;
    }
    ///Write variables as a formatted string
    void display( std::ostream& os = std::cout ) const
    {
        os << "Solovev Geometrical parameters are: \n"
            <<" A               = "<<A<<"\n"
            <<" Prefactor Psi   = "<<pp<<"\n"
            <<" Prefactor I     = "<<pi<<"\n";
        for( unsigned i=0; i<12; i++)
            os<<" c"<<i+1<<"\t\t = "<<c[i]<<"\n";

        os  <<" R0            = "<<R_0<<"\n"
            <<" a             = "<<a<<"\n"
            <<" epsilon_a     = "<<a/R_0<<"\n"
            <<" description   = "<<description<<"\n"
            <<" elongation    = "<<elongation<<"\n"
            <<" triangularity = "<<triangularity<<"\n";
        os << std::flush;

    }
};
} //namespace solovev
} //namespace geo
} //namespace dg
