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

//#define DEBUG_MINIMIZATION_ON
//#define DEBUG_TIELINES_ON
#define XDIM 401 ///< The length of the D1Q3 lattice.
#define isEqual(a, b, threshold) ( (fabs((a)-(b)) < (threshold)) ? 1 : 0 ) ///< Equality tester for all the doubles out there.
#define MAXSTRLEN 255 ///< Maximum length of a file name

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
double theoreticalMuAArray[XDIM];		///< Component A, theoretical chemical potential from free energy minimization
double theoreticalMuBArray[XDIM];		///< Component B, theoretical chemical potential from free energy minimization
double pressure[XDIM];					///< pressure of the mixture recovered by simulation
double theoreticalPressureArray[XDIM];	///< theoretical pressure from free energy minimization

// Forces on each component; potential and friction
double F1[XDIM];			///< Component A total force
double F2[XDIM];			///< Component B total force
double friction1[XDIM];		///< Component A friction (momentum transferred to the B component)
double friction2[XDIM];		///< Component B friction (momentum transferred to the A component)
double gradMuForce1[XDIM];  ///< Component A conservative force from the chemical potential gradient
double gradMuForce2[XDIM];  ///< Component B conservative force from the chemical potential gradient

double threePhaseRegion[6]; ///< coordinates of 3-phase region from free energy minimization

/// Free energy values used during minimization
struct FreeEnergy {
	double freeEnergyTotal;	///< Total free energy
	double freeEnergy1;		///< Phase 1 free energy
	double freeEnergy2;		///< Phase 2 free energy
	double freeEnergy3;		///< Phase 3 free energy
};

/// Binodal point (A,B) coordinates
struct BinodalPoint {
	double rhoA;
	double rhoB;
};


//
// Externs
//

extern char *dataDirectory;	///< Directory for log/data files and gnuplot scripts

// Domain and phase information
extern double theoreticalDensityA1;			///< The expected A-component density in the 1st phase from a minimized free energy
extern double theoreticalDensityA2;			///< The expected A-component density in the 2nd phase from a minimized free energy
extern double theoreticalDensityA3;			///< The expected A-component density in the 3rd phase from a minimized free energy
extern double theoreticalDensityB1;			///< The expected B-component density in the 1st phase from a minimized free energy
extern double theoreticalDensityB2;			///< The expected B-component density in the 2nd phase from a minimized free energy
extern double theoreticalDensityB3;			///< The expected B-component density in the 3rd phase from a minimized free energy
extern double theoreticalVolume1;			///< The expected volume occupied by the 1st phase from a minimized free energy
extern double theoreticalVolume2;			///< The expected volume occupied by the 2nd phase from a minimized free energy
extern double theoreticalVolume3;			///< The expected volume occupied by the 3rd phase from a minimized free energy
extern double theoreticalPressure;			///< The expected bulk pressure calculated with densities from a minimized free energy
extern double theoreticalMuA;				///< The expected A-component bulk chemical potential calculated with densities from a minimized free energy
extern double theoreticalMuB;				///< The expected B-component bulk chemical potential calculated with densities from a minimized free energy
extern double interfaceWidth;				///< The interface width of the initialized density profile
extern double interfaceWidthForGivenKappa;	///< The interface width of the optimal tanh profile fit to a near-equilbrium LB simulation
extern int theoreticalPhases;				///< The number of phases (2 or 3) in the minimized free energy for an (A,B) test point
extern int initializeRandomComponents; 		///< Chooser for which component (0 = both, 1 = A, 2 = B) to initialize with a random profile
extern int domain1; 						///< Lattice point that is the center of the A-rich phase of the initialized density profile
extern int domain1Width;					///< Width in lattice units of the A-rich phase
extern int domain2; 						///< Lattice point that is the center of the first vapor phase separating A-rich and B-rich phases
extern int domain2Width;					///< Width in lattice units of the first vapor phase
extern int domain3; 						///< Lattice point that is the center of the B-rich phase of the initialized density profile
extern int domain3Width;					///< Width in lattice units of the B-rich phase
extern int domain4;							///< Lattice point that is the center of the second vapor phase separating A-rich and B-rich phases
extern int domain4Width;					///< Width in lattice units of the second vapor phase
extern int phase1Index;						///< Lattice space from which to read the A- and B-component densities of the first phase in a 2-phase LB simulation
extern int phase2Index;						///< Lattice space from which to read the A- and B-component densities of the second phase in a 2-phase LB simulation

// Physical properties common to minimization and LB simulations
extern double nA0;					///< Total A-component material (per lattice site for LB simulations)
extern double nB0;					///< Total B-component material (per lattice site for LB simulations)
extern double nAIntegrated;			///< Amount of A-component per lattice site during a LB simulation
extern double nBIntegrated;			///< Amount of B-component per lattice site during a LB simulation
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
extern double g;          			///< Gravitational acceleration term
extern double lambda;     			///< Friction coefficient

// Free energy minimization control
extern double minimizationParticlesStepSize;	///< Resolution of the density landscape for free energy minimization
extern double minimizationStepFactor;			///< Initial step size for free energy minimization density trials (fraction of the density landscape resolution)
extern int sortABinodalByIncreasingB;			///< Control sorting logic for constructing A-component liquid-vapor binodal
extern int sortBBinodalByIncreasingA;			///< Control sorting logic for constructing B-component liquid-vapor binodal
extern int threePhaseRegionExists;				///< Flag to indicate if a phase diagram has 3-phase behavior

// 3-phase region properties and control
extern double rhoAThreePhaseRegion;			///< A-component coordinate of the point used to define the three-phase region
extern double rhoBThreePhaseRegion;			///< B-component coordinate of the point used to define the three-phase region
extern double maxArea;						///< Metric used to define the 3-phase region of a phase diagram
extern double metastableThreshold;			///< Threshold for metric used to separate unstable 2-phase from metastable or 3-phase behavior
extern double metastableThresholdFactor;	///< Metric coefficient used to separate 3-phase behavior from metastable 2-phase behavior
extern int setThreePhaseRegion;				///< Flag to indicate if 3-phase data needs to be read in
extern int setLineAPoint;					///< Flag to indicate point 2 of the line separating binary liquid behavior needs to be set
extern int setLineBPoint;					///< Flag to indicate point 1 of the line separating binary liquid behavior needs to be set

// LB simulation initialization and runtime control
extern double kappa;								///< Interfacial free energy term
extern double kappaFactor;							///< Correction factor to fit kappa to an interface width
extern double gammaP; 								///< Coefficient to control the pressure applied per time step
extern double gammaMu;								///< Coefficient to control the chemical potential applied (i.e. rate of forcing) per time step
extern double gammaFactor;							///< Correction factor to change gammaMu to stabilize a simulation
extern double endGammaFactor;						///< Maximum gamma factor for the search space to find a stabilizing gammaMu
extern int useTwoOrThreePhaseInitialization;		///< Choose 2 or 3 phases for tanh profile initialization
extern int useStepOrRandomInitialization;			///< Flag to choose random or tanh profile for initialization
extern int useTheoreticalDensities;					///< Flag to choose to use densities from free energy minimization or manually input densities
extern int useTheoreticalVolumina;					///< Flag to choose to use volumina from free energy minimization or equal volumina for each phase
extern int useTheoreticalPhase1;					///< Flag to choose theoretical phase 1 when initializing a 2-phase density profile
extern int useTheoreticalPhase2;					///< Flag to choose theoretical phase 2 when initializing a 2-phase density profile
extern int useTheoreticalPhase3;					///< Flag to choose theoretical phase 3 when initializing a 2-phase density profile
extern int suppressTheoreticalPhase1;				///< Flag to suppress the first liquid phase when initializing a 3-phase density profile
extern int suppressTheoreticalPhase2;				///< Flag to suppress the second liquid phase when initializing a 3-phase density profile
extern int isolateFourthPhase;						///< Flag to control the arrangement of domains when initializing the 4-phase behavior attempt
extern int overrideMinimumInterfaceWidth;			///< Flag to disregard minimum interface width of 2 for density profile initialization
extern int goodInterfaceFit;						///< Flag for initialization process to indicate if an optimal interface width was found for the given kappa
extern int fitInterfaceWidthToKappa;				///< State machine control for the stages of fitting an interface width to a given kappa
extern int determineInterfaceWidthAutomatically;	///< Flag to use an interface width that is manually or automatically defined
extern int determineKappaAutomatically;				///< Flag to use a kappa that is manually defined or automatically determined from single-component theory
extern int kappaFactorDetermined;					///< Flag to indicate if an interface width has been fit for the given kappa or if a new fit is needed
extern int useChemicalPotentialForcingMethod;		///< Choose between "nid" and "log" chemical potential forcing methods
extern int useBoundaryConditionsPeriodic;			///< Flag to toggle between periodic or wall boundaries
extern int usePressurePartitionFunction;			///< Flag to calculate pressure based on a partition function or a constructed tensor
extern int useMuVdwInteraction;						///< Flag to use _vdwInteractionFactor_ as a kappa coefficient for cross-term (A-B) forcing
extern int setPhaseDiagramTwoPhaseRegion;			///< Flag to include/exclude points from the 2-phase region when building a phase diagram from LB simulations
extern int setPhaseDiagramThreePhaseRegion;			///< Flag to include/exclude points from the 3-phase region when building a phase diagram from LB simulations
extern int checkTwoPhaseMetastableLBPoints;			///< Flag to include/exclude metastable points from the 3-phase region in the LB phase diagram
extern int checkThreePhaseLBPoints;					///< Flag to include/exclude points with anticipated 3-phase behavior in the LB phase diagram
extern int applyMomentumCorrection;					///< Flag to apply force corrections to account for discrete errors that may accumulate

extern double childRhoA;		///< A-component coordinate of a point for which the parent minimization or LB simulation point is desired
extern double childRhoB;		///< B-component coordinate of a point for which the parent minimization or LB simulation point is desired
extern double singlePointRhoA;	///< A-component coordinate of a specific, single point to minimize
extern double singlePointRhoB;	///< B-component coordinate of a specific, single point to minimize

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
int calculateMinimizationPath(double, double, double *, double *);
int stableSimulation(double *);
double gradient(double *, int);
double laplace(double *,int);
void generatePhaseDiagramBinodalLines();
void generateMinimizedPressureAndMuDeviations();
void generatePhaseDiagramRhoAVsRhoB();
void sortMetastableMinimizationResults();
void sortMetastableLBResults();
void generateEigenvalueMap();
void generateInterfaceWidthFitValues();
void printComponentDensities();
void printComponentMaxMin();
void minimizeSinglePoint();
//void calculatePhaseDiagramRhoAVsRhoBTwoPhasesTheoretical();
void calculatePhaseDiagramRhoAVsRhoBThreePhasesTheoretical();
void exploreVDWShieldRegion();
void calculateKandSParameters();
void logVDWParameters();
void printMinimizationParentPoint();
void printLBParentPoint();
void generateFreeEnergyData();
void generateFreeEnergySlices();

// Shared from module LB_D1Q3_2-components.c
void iteration();
void setCollisionForcingNewChemicalPotentialGradient();

// Shared from module LB_Files.c
FILE * openFile(char *, char *, char *);
int setDataDirectory();
void gnuplotTwoComponentTwoPhase();
void gnuplotTwoComponentThreePhase();
void gnuplotEigenstuffData();
void gnuplotInterfaceWidthFit();
void gnuplotPublicationGraphics();
void logVDWParameters();
void writeThreePhaseRegionCoordinates();
void readThreePhaseRegionCoordinates();
void getDensityProfile();
void getPressureProfile();
void getChemicalPotentialProfile();

// Shared from module LB_GUI.c
void GUI();

// Shared from module LB_Initialize.c
void init();
void initialize();
void setInitializeRandom();
void setInitializeTheoreticalTwoPhases();
void setInitializeTheoreticalThreePhases();
void setInitializeTheoreticalFourPhases();
void setLBInitializationProfile();
int fitGammaFactor();


#endif /* LB_D1Q3_2_COMPONENTS_H_ */


//struct SpinodalPoint {
//	int spinodalLevel;
//	double particlesATotal;
//	double particlesBTotal;
//	double particlesA1;
//	double particlesB1;
//	double volume1;
//	double particlesA2;
//	double particlesB2;
//	double volume2;
//	double eigenvalue1;
//	double eigenvalue2;
//};
