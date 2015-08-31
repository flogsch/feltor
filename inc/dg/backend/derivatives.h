#ifndef _DG_DERIVATIVES_CUH_
#define _DG_DERIVATIVES_CUH_

#include "grid.h"
#include "dx.h"

/*! @file 
  
  Convenience functions to create 2D derivatives
  */
namespace dg{


/**
 * @brief Contains functions used for matrix creation
 */
namespace create{

///@addtogroup highlevel
///@{

//dx, dy, jumpX, jumpY

/**
 * @brief Create 2d derivative in x-direction
 *
 * @param g The grid on which to create dx
 * @param bcx The boundary condition
 * @param dir The direction of the first derivative
 *
 * @return A host matrix 
 */
EllSparseBlockMat dx( const Grid2d<double>& g, bc bcx, direction dir = centered)
{
    EllSparseBlockMat dx;
    dx = dx_normed( g.n(), g.Nx(), g.hx(), bcx, dir);
    dx.left = g.n()*g.Ny();
    return dx;
}

/**
 * @brief Create 2d derivative in x-direction
 *
 * @param g The grid on which to create dx (boundary condition is taken from here)
 * @param dir The direction of the first derivative
 *
 * @return A host matrix
 */
EllSparseBlockMat dx( const Grid2d<double>& g, direction dir = centered) { return dx( g, g.bcx(), dir);}

/**
 * @brief Create 2d derivative in y-direction
 *
 * @param g The grid on which to create dy
 * @param bcy The boundary condition
 * @param dir The direction of the first derivative
 *
 * @return A host matrix
 */
EllSparseBlockMat dy( const Grid2d<double>& g, bc bcy, direction dir = centered)
{
    EllSparseBlockMat dy;
    dy = dx_normed( g.n(), g.Ny(), g.hy(), bcy, dir);
    dy.right = g.n()*g.Nx();
    return dy;
}

/**
 * @brief Create 2d derivative in y-direction
 *
 * @param g The grid on which to create dy (boundary condition is taken from here)
 * @param dir The direction of the first derivative
 *
 * @return A host matrix 
 */
EllSparseBlockMat dy( const Grid2d<double>& g, direction dir = centered){ return dy( g, g.bcy(), dir);}

/**
 * @brief Matrix that contains 2d jump terms in X direction
 *
 * @param g grid
 * @param bcx boundary condition in x
 *
 * @return A host matrix 
 */
EllSparseBlockMat jumpX( const Grid2d<double>& g, bc bcx)
{
    EllSparseBlockMat jx;
    jx = jump( g.n(), g.Nx(), g.hx(), bcx);
    jx.left = g.n()*g.Ny();
    return jx;
}

/**
 * @brief Matrix that contains 2d jump terms in Y direction
 *
 * @param g grid
 * @param bcy boundary condition in y
 *
 * @return A host matrix 
 */
EllSparseBlockMat jumpY( const Grid2d<double>& g, bc bcy)
{
    EllSparseBlockMat jy;
    jy = jump( g.n(), g.Ny(), g.hy(), bcy);
    jy.right = g.n()*g.Nx();
    return jy;
}

/**
 * @brief Matrix that contains 2d jump terms in X direction taking boundary conditions from the grid
 *
 * @param g grid
 *
 * @return A host matrix 
 */
EllSparseBlockMat jumpX( const Grid2d<double>& g)
{
    return jumpX( g, g.bcx());
}

/**
 * @brief Matrix that contains 2d jump terms in Y direction taking boundary conditions from the grid
 *
 * @param g grid
 *
 * @return A host matrix 
 */
EllSparseBlockMat jumpY( const Grid2d<double>& g)
{
    return jumpY( g, g.bcy());
}

///////////////////////////////////////////3D VERSIONS//////////////////////
//jumpX, jumpY, jumpZ, dx, dy, dz
/**
 * @brief Matrix that contains jump terms in X direction in 3D
 *
 * @param g The 3D grid
 * @param bcx boundary condition in x
 *
 * @return A host matrix 
 */
EllSparseBlockMat jumpX( const Grid3d<double>& g, bc bcx)
{
    EllSparseBlockMat jx;
    jx = jump( g.n(), g.Nx(), g.hx(), bcx);
    jx.left = g.n()*g.Ny()*g.Nz();
    return jx;
}

/**
 * @brief Matrix that contains jump terms in Y direction in 3D
 *
 * @param g The 3D grid
 * @param bcy boundary condition in y
 *
 * @return A host matrix 
 */
EllSparseBlockMat jumpY( const Grid3d<double>& g, bc bcy)
{
    EllSparseBlockMat jy;
    jy = jump( g.n(), g.Ny(), g.hy(), bcy);
    jy.right = g.n()*g.Nx();
    jy.left = g.Nz();
    return jy;
}

/**
 * @brief Matrix that contains jump terms in Z direction in 3D
 *
 * @param g The 3D grid
 * @param bcz boundary condition in z
 *
 * @return A host matrix 
 */
EllSparseBlockMat jumpZ( const Grid3d<double>& g, bc bcz)
{
    EllSparseBlockMat jz;
    jz = jump( 1, g.Nz(), g.hz(), bcz);
    jz.right = g.n()*g.Nx()*g.n()*g.Ny();
    return jz;
}

/**
 * @brief Matrix that contains 3d jump terms in X direction taking boundary conditions from the grid
 *
 * @param g grid
 *
 * @return A host matrix
 */
EllSparseBlockMat jumpX( const Grid3d<double>& g)
{
    return jumpX( g, g.bcx());
}

/**
 * @brief Matrix that contains 3d jump terms in Y direction taking boundary conditions from the grid
 *
 * @param g grid
 *
 * @return A host matrix
 */
EllSparseBlockMat jumpY( const Grid3d<double>& g)
{
    return jumpY( g, g.bcy());
}

/**
 * @brief Matrix that contains 3d jump terms in Z direction taking boundary conditions from the grid
 *
 * @param g grid
 *
 * @return A host matrix
 */
EllSparseBlockMat jumpZ( const Grid3d<double>& g)
{
    return jumpZ( g, g.bcz());
}


/**
 * @brief Create 3d derivative in x-direction
 *
 * @param g The grid on which to create dx
 * @param bcx The boundary condition
 * @param dir The direction of the first derivative
 *
 * @return A host matrix 
 */
EllSparseBlockMat dx( const Grid3d<double>& g, bc bcx, direction dir = centered)
{
    EllSparseBlockMat dx;
    dx = dx_normed( g.n(), g.Nx(), g.hx(), bcx, dir);
    dx.left = g.n()*g.Ny()*g.Nz();
    return dx;
}

/**
 * @brief Create 3d derivative in x-direction
 *
 * @param g The grid on which to create dx (boundary condition is taken from here)
 * @param dir The direction of the first derivative
 *
 * @return A host matrix 
 */
EllSparseBlockMat dx( const Grid3d<double>& g, direction dir = centered) { return dx( g, g.bcx(), dir);}

/**
 * @brief Create 3d derivative in y-direction
 *
 * @param g The grid on which to create dy
 * @param bcy The boundary condition
 * @param dir The direction of the first derivative
 *
 * @return A host matrix 
 */
EllSparseBlockMat dy( const Grid3d<double>& g, bc bcy, direction dir = centered)
{
    EllSparseBlockMat dy;
    dy = dx_normed( g.n(), g.Ny(), g.hy(), bcy, dir);
    dy.right = g.n()*g.Nx();
    dy.left = g.Nz();
    return dy;
}

/**
 * @brief Create 3d derivative in y-direction
 *
 * @param g The grid on which to create dy (boundary condition is taken from here)
 * @param dir The direction of the first derivative
 *
 * @return A host matrix 
 */
EllSparseBlockMat dy( const Grid3d<double>& g, direction dir = centered){ return dy( g, g.bcy(), dir);}

/**
 * @brief Create 3d derivative in z-direction
 *
 * @param g The grid on which to create dz
 * @param bcz The boundary condition
 * @param dir The direction of the stencil
 *
 * @return A host matrix 
 */
EllSparseBlockMat dz( const Grid3d<double>& g, bc bcz, direction dir = centered)
{
    EllSparseBlockMat dz;
    dz = dx_normed( 1, g.Nz(), g.hz(), bcz, dir);
    dz.right = g.n()*g.n()*g.Nx()*g.Ny();
    return dz;

}

/**
 * @brief Create 3d derivative in z-direction
 *
 * @param g The grid on which to create dz (boundary condition is taken from here)
 * @param dir The direction of the stencil
 *
 * @return A host matrix 
 */
EllSparseBlockMat dz( const Grid3d<double>& g, direction dir = centered){ return dz( g, g.bcz(), dir);}



///@}

} //namespace create

} //namespace dg

#endif//_DG_DERIVATIVES_CUH_
