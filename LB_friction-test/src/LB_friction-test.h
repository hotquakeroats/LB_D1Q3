/**
 * @file LB_D1Q3_2-components.h
 * @author Kent S. Ridl
 * @date 2 December 2017
 *
 * The module _LB_D1Q3_2-components.h_ contains the declarations for all global variables, data structures, and shared function headers.
 */

#ifndef LB_D1Q3_2_COMPONENTS_H_
#define LB_D1Q3_2_COMPONENTS_H_

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>
#include <sys/stat.h>

#include "mygraph.h"

#define XDIM 401 ///< The length of the D1Q3 lattice.
#define isEqual(a, b, threshold) ( (fabs((a)-(b)) < (threshold)) ? 1 : 0 ) ///< Equality tester for all the doubles out there.

// D1Q3 equilibrium distribution of each component velocities
double f1_0[XDIM]; ///< Component A, zero velocity
double f1_1[XDIM]; ///< Component A, +1 velocity
double f1_2[XDIM]; ///< Component A, -1 velocity
double f2_0[XDIM]; ///< Component B, zero velocity
double f2_1[XDIM]; ///< Component B, +1 velocity
double f2_2[XDIM]; ///< Component B, -1 velocity
double psi1[XDIM]; ///< Component A, 4th order correction term
double psi2[XDIM]; ///< Component B, 4th order correction term

// Number densities of each component
double n1[XDIM]; 		///< Component A mass/density
double n2[XDIM]; 		///< Component B mass/density
double n[XDIM];  		///< Total mass/density (\f$\rho\f$)
double nReduced[XDIM];	///< Mass term for friction force calculations \f$\left(\frac{n1*n2}{n1+n2}\right)\f$

// Component and bulk velocities
double u1[XDIM]; 	///< Component A velocity (not a true hydrodynamic variable)
double u2[XDIM]; 	///< Component B velocity (not a true hydrodynamic variable)
double uHat1[XDIM];	///< Component A corrected velocity
double uHat2[XDIM];	///< Component B corrected velocity
double u[XDIM]; 	///< Mean fluid velocity

// Component thermodynamic properties
double mu1[XDIM]; 						///< Component A, chemical potential recovered by simulation
double muNonIdeal1[XDIM]; 				///< Component A, non-ideal part of the chemical potential (\f$\mu - \rho ln(\rho)\f$)
double mu2[XDIM]; 						///< Component B, chemical potential recovered by simulation
double muNonIdeal2[XDIM]; 				///< Component B, non-ideal part of the chemical potential (\f$\mu - \rho ln(\rho)\f$)
double pressure[XDIM];					///< pressure of the mixture recovered by simulation

// Forces on each component; potential and friction
double F1[XDIM];			///< Component A total force
double F2[XDIM];			///< Component B total force
double friction1[XDIM];		///< Component A friction (momentum transferred to the B component)
double friction2[XDIM];		///< Component B friction (momentum transferred to the A component)
double gradMuForce1[XDIM];  ///< Component A conservative force from the chemical potential gradient
double gradMuForce2[XDIM];  ///< Component B conservative force from the chemical potential gradient


//
// Externs
//

// Physical properties common to minimization and LB simulations
extern double nA0;					///< Total A-component material (per lattice site for LB simulations)
extern double nB0;					///< Total B-component material (per lattice site for LB simulations)
extern double uA0;
extern double uB0;
extern double theta;				///< Lattice temperature (k_B*T)
extern double aA;					///< van der Waals constant for the A component self-attraction
extern double aB;					///< van der Waals constant for the B component self-attraction
extern double aAB;					///< van der Waals constant for the A-B interaction
extern double bA;					///< van der Waals constant for the A component excluded volume
extern double bB;					///< van der Waals constant for the B component excluded volume
extern double vdwInteractionFactor;	///< Coefficient of the free energy A-B interaction (1.0 is a neutral interaction, < 1.0 is repulsive, > 1.0 is attractive)
extern double Amp;					///< Amplitude of the random term for the mixed density profile initialization
extern double tcA;					///< Critical temperature of the A component
extern double ncA;					///< Critical density of the A component
extern double pcA;					///< Critical pressure of the A component
extern double tcB;					///< Critical temperature of the B component
extern double ncB;					///< Critical density of the B component
extern double pcB;					///< Critical pressure of the B component
extern double tau;					///< Lattice relaxation constant
extern double oneOverTau;   		///< Inverse of the lattice relaxation constant
extern double lambda;     			///< Friction coefficient

// LB simulation initialization and runtime control
extern double kappa;								///< Interfacial free energy term
extern double gammaP; 								///< Coefficient to control the pressure applied per time step
extern double gammaMu;								///< Coefficient to control the chemical potential applied (i.e. rate of forcing) per time step
extern int useChemicalPotentialForcingMethod;		///< Choose between "nid" and "log" chemical potential forcing methods
extern int useBoundaryConditionsPeriodic;			///< Flag to toggle between periodic or wall boundaries
extern int usePressurePartitionFunction;			///< Flag to calculate pressure based on a partition function or a constructed tensor
extern int initializeRandomComponents;

// GUI
extern int next;					///< Step a LB simulation a specified number of iterations (_Repeat_); does not update plots during simulation execution
extern int Pause;					///< Start/stop a LB simulation with no bound on the number of iterations to run
extern int run;						///< Step a LB simulation a specified number of iterations (_phase_iterations_); updates plots during simulation execution
extern int done;					///< Program exit flag
extern int Repeat;					///< Number of iterations of a LB simulation to run before updating plots
extern int iterations;				///< Total accumulated number of iterations of a LB simulation
extern int tmp_phase_iterations;	///< Iteration number at which to stop the LB simulation
extern int phase_iterations;		///< Number of iterations to run a LB simulation


//
// Function Prototypes
//

// Shared from module LB_Analyze.c
double gradient(double *, int);
double laplace(double *,int);

// Shared from module LB_D1Q3_2-components.c
void iteration();
void setCollisionForcingNewChemicalPotentialGradient();

// Shared from module LB_GUI.c
void GUI();

// Shared from module LB_Initialize.c
void init();
void initialize();
void setInitializeRandom();
void setInitializeTheoreticalTwoPhases();


#endif /* LB_D1Q3_2_COMPONENTS_H_ */
