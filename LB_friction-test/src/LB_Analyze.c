/**
 * @file LB_Analyze.c
 * @author Kent S. Ridl
 * @date 7 July 2016 (last modified: 26 November 2017)
 *
 * The module _LB_Analyze.c_ contains the code used to set theoretical expectations for the lattice Boltzmann simulations and to process simulation results.
 * It includes everything from the free energy minimizer, gradient stencils, logging functions, etc.
 */

#include "LB_friction-test.h"


/**
 * @brief Function _gradient_ returns a spatial gradient at a specified location.
 *
 * This function calculates a spatial gradient for both periodic and solid boundaries of a simulation domain.
 * It implements 2nd order forward, backward, and central differences.
 * It does not implement a halo to account for gradients at edges of the domain.
 *
 * The global flag _useBoundaryConditionsPeriodic_ toggles between the periodic and solid boundary calculations.
 * - 1: uses periodic boundary conditions (central difference only)
 * - 0: uses solid walls at domain boundaries (forward, backward, central differences allowed)
 *
 * @param [in] *var An input array (pointer) representing values across space.
 * @param [in] i The array index of the location for which the spatial gradient is to be calculated.
 *
 * @return The calculated spatial gradient.
 */
double gradient(double *var, int i) {
	double result = 0;
	int ipp = i+2, ip = i+1; /// use indices for 2 cells before and after cell i
	int im = i-1, imm = i-2;

	if (useBoundaryConditionsPeriodic) { // gradient for a lattice with periodic BCs
		if (ip == XDIM) ip = 0;
		else if (im == -1) im = XDIM - 1;
		result = 1./2. * (var[ip]-var[im]);
	}
	else { // gradient for a lattice with walls at the boundaries; adjust for existence of the walls with forward/backward differencing
		if (ip == XDIM)  result = 1./2. * (var[imm] - 4.*var[im] + 3.*var[i]); // backward difference
		else if (im == -1) result = 1./2. * (-3.*var[i] + 4.*var[ip] - var[ipp]); // forward difference
		else result = 1./2. * (var[ip]-var[im]); // central difference... normal case
	}
	return result;
} // end function gradient()


/**
 * @brief Function _laplace_ returns a spatial laplacian at a specified location.
 *
 * This function calculates a laplacian for both periodic and solid boundaries of a simulation domain.
 * It implements 2nd order forward, backward, and central differences.
 * It does not implement a halo to account for the laplacian at edges of the domain.
 *
 * The global flag _useBoundaryConditionsPeriodic_ toggles between the periodic and solid boundary calculations.
 * - 1: uses periodic boundary conditions (central difference only)
 * - 0: uses solid walls at domain boundaries (forward, backward, central differences allowed)
 *
 * @param [in] *var An input array (pointer) representing values across space.
 * @param [in] i The array index of the location for which the laplacian is to be calculated.
 *
 * @return The calculated laplacian.
 */
double laplace(double *var, int i) {
	double result = 0;
	int ippp = i+3, ipp = i+2, ip = i+1; /// use indices for 3 cells before and after cell i
	int im = i-1, imm = i-2, immm = i-3;

	if (useBoundaryConditionsPeriodic) { // gradient for lattice with periodic BCs; instead of a halo, reset the out-of-bounds index to opposite end the lattice
		if (ip == XDIM) ip = 0;
		else if (im == -1) im = XDIM - 1;
		result = var[ip] - 2.*var[i] + var[im];
	}
	else { // gradient for a lattice with walls at the boundaries; adjust for existence of the walls with forward/backward differencing
		if (ip == XDIM) result = -1.*var[immm] + 4.*var[imm] -5.*var[im] + 2.*var[i]; // backward difference
		else if (im == -1) result = 2.*var[i] - 5.*var[ip] + 4.*var[ipp] - var[ippp]; // forward difference
		else result = var[ip] - 2.*var[i] + var[im]; // central difference... normal case
	}
	return result;
} // end function laplace()

