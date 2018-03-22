/**
 * @file LB_Initialize.c
 * @author Kent S. Ridl
 * @date 7 July 2016 (last modified: 28 November 2017)
 *
 * The module _LB_Initialize.c_ contains the code used to initialize lattice Boltzmann simulations.  This includes fetching the theoretical results from
 * the free energy minimization for a test point, setting the near-equilibrium phases for 2-phase and 3-phase behavior, and optimizing simulation
 * parameters for the anticipated interface width.
 */

#include "LB_friction-test.h"

//TODO: add ability to initialize 3-phase point with only 3 interfaces (a-vapor, b-vapor, a-b) instead of 4
//TODO: may need to break kappa symmetry
//TODO: implement a choice to do the full golden search for an optimal value or to simply find one that works and go (to save time)

void (*initializeProfile)(int); ///< function pointer to choose the method to initialize the lattice density profile (random, 2-phase, 3-phase)


/**
 * @brief Function _initializeRandom_ initializes a fully mixed 2-component simulation (i.e. equal densities across the lattice) with an additional noise term.
 *
 * @note The global value _Amp_ controls the amplitude of the noise term added to the density profile.
 * @note The global values _nA0_ and _nB0_ set the total mass/densities to be initialized.
 *
 * @param [in] i The lattice site to be initialized.
 */
void initializeRandom(int i) {
	if (nA0 < Amp || nB0 < Amp) Amp = 0.1 * (nA0 < nB0 ? nA0 : nB0); // make sure fluctuations don't give negative particle counts

	if (initializeRandomComponents == 1 || initializeRandomComponents == 2) Amp = 0.0; // no noise if only 1 component is random (only done if spoofing zero)
	if (initializeRandomComponents == 1 || initializeRandomComponents == 0) n1[i] = nA0 + Amp*(1.0*rand()/RAND_MAX-0.5);
	if (initializeRandomComponents == 2 || initializeRandomComponents == 0) n2[i] = nB0 + Amp*(1.0*rand()/RAND_MAX-0.5);
} // end function initializeRandom()



/**
 * @brief Function _initializeTheoreticalTwoPhases_ initializes two density domains for a lattice Boltzmann simulation.
 *
 * This function initializes two density domains for a simulation, one A-rich domain and one B-rich domain, using tanh profiles as an approximation for the
 * domain shape.  The bulk density values are taken from the free energy minimization of the (A,B) pair that is to be simulated.  It's nominal use is for
 * simulations where the expected behavior is phase separation into 2 phases.
 *
 * @note The global values _nA0_ and _nB0_ set the total mass/densities to be initialized.
 *
 * @param [in] i The lattice site to be initialized.
 */
//void initializeTheoreticalTwoPhases(int i) {
//	int interface = 0.5 * XDIM;
//	double transition = 0.0;
//	double rhoA1 = theoreticalDensityA1, rhoB1 = theoreticalDensityB1;
//	double rhoA2 = theoreticalDensityA2, rhoB2 = theoreticalDensityB2;
//
//	domain4 = XDIM;
//	if (i < interface-0.25*XDIM) {
//		transition = 0.5 + 0.5*tanh((double)(i%XDIM)/interfaceWidth);
//		n1[i] = (1.0-transition)*rhoA2 + transition*rhoA1;
//		n2[i] = (1.0-transition)*rhoB2 + transition*rhoB1;
//	}
//	else if (i >= interface-0.25*XDIM && i <= interface+0.25*XDIM) {
//		transition = 0.5 + 0.5*tanh((double)(i-interface)/interfaceWidth);
//		n1[i] = (1.0-transition)*rhoA1 + transition*rhoA2;
//		n2[i] = (1.0-transition)*rhoB1 + transition*rhoB2;
//	}
//	else {
//		transition = 0.5 + 0.5*tanh((double)((i-domain4)%XDIM)/interfaceWidth);
//		n1[i] = (1.0-transition)*rhoA2 + transition*rhoA1;
//		n2[i] = (1.0-transition)*rhoB2 + transition*rhoB1;
//	}
//
//	// If one component is practically zero, use the random initializer on only that component
//	if (nA0 < 1e-10) {
//		initializeRandomComponents = 1;
//		initializeRandom(i);
//		initializeRandomComponents = 0;
//	}
//	if (nB0 < 1e-10) {
//		initializeRandomComponents = 2;
//		initializeRandom(i);
//		initializeRandomComponents = 0;
//	}
//
//	// Set indices from which to read equilibrium densities in 2-phase simulations
//	phase1Index = interface-0.25*XDIM;
//	phase2Index = interface+0.25*XDIM;
//} // end function initializeTheoreticalTwoPhases()


/**
 * @brief Function _setInitializeRandom_ sets the simulation to initialize a uniform mixture.
 */
void setInitializeRandom() {
	initializeProfile = initializeRandom;
	initialize();
} // end function setInitializeRandom()


/**
 * @brief Function _setInitializeTheoreticalTwoPhases_ sets the simulation to initialize 2-phase behavior.
 */
//void setInitializeTheoreticalTwoPhases() {
//	initializeProfile = initializeTheoreticalTwoPhases;
//	initialize();
//} // end function setInitializeTheoreticalTwoPhases()


/**
 * @brief Function _initialize_ is the controlling function to initialize a lattice Boltzmann simulation.
 *
 * This function is the controlling entity for all lattice Boltzmann simulations.  When initializing 2-phase or 3-phase behavior, the default is to use bulk
 * densities from a free energy minimization along with a manually set value for kappa and an interface width that is automatically fit to that kappa.
 *
 * @note The global flag _useTheoreticalDensities_ controls if the minimization results for an (A,B) test pair are to be used to initialize bulk phases.
 * @note The global flag _determineInterfaceWidthAutomatically_ controls if the initialized interface width is automatically or manually set.  If a manual
 * interface width is chosen, the initialization forces kappa to use a manual value as well.
 * @note The global flag _determineKappaAutomatically_ controls if the initialized kappa value is automatically or manually set.
 * @note Once the interface width is automatically fit to a specific manual kappa value, subsequent initializations use the same value of the global parameter
 * _kappaFactor_ to change the interface with when kappa is changed (width = kappaFactor * sqrt(kappa)).  If this results unstable simulations for kappas other
 * than the one used for the interface fit, setting the global flag _kappaFactorDetermined_ to "off" will force a new interface fit for a new kappa value.
 * @note The value for gammaMu is always automatically determined.  However, it may be manually set using a new value for the global parameter _gammaFactor_.
 */
void initialize() {

	// tcA, tcB, and ncB vary as the degrees of freedom
	pcA = 3.*tcA/8.; // determine pcA
	ncA = pcA / ((3./8.)*tcA); // fix ncA to be 1
	pcB = (3./8.)*tcB*ncB; // determine pcB

	// Reset the values of the VDW constants for each component
	aA = (27./64.)*(tcA*tcA/pcA);
	aB = (27./64.)*(tcB*tcB/pcB);
	aAB = sqrt(aA*aB) * vdwInteractionFactor;
	bA = tcA/(8.*pcA);
	bB = tcB/(8.*pcB);

	iterations = 0;
	tmp_phase_iterations = phase_iterations;

	// Loop to initialize lattice values
	for (int i = 0; i < XDIM; i++) {
		initializeProfile(i);
		n[i] = n1[i] + n2[i];
		nReduced[i] = n1[i]*n2[i] / n[i];
		u1[i] = uA0; //0.0;
		u2[i] = uB0; //0.0;
		uHat1[i] = u1[i]; //0.0;
		uHat2[i] = u2[i]; //0.0;
		u[i] = (n1[i]*u1[i] + n2[i]*u2[i]) / n[i]; //0.0;

		// Initialize 1st component
		f1_0[i] = n1[i] - n1[i]*theta - n1[i]*u1[i]*u1[i];            		// zero velocity
		f1_1[i] = 0.5 * (n1[i]*u1[i]*u1[i] + n1[i]*u1[i] + n1[i]*theta); 	// +1 velocity (moving right)
		f1_2[i] = 0.5 * (n1[i]*u1[i]*u1[i] - n1[i]*u1[i] + n1[i]*theta); 	// -1 velocity (moving left)
		psi1[i] = 0.0;
		mu1[i] = 0.0;
		muNonIdeal1[i] = 0.0;

		// Initialize 2nd component
		f2_0[i] = n2[i] - n2[i]*theta - n2[i]*u2[i]*u2[i];            		// zero velocity
		f2_1[i] = 0.5 * (n2[i]*u2[i]*u2[i] + n2[i]*u2[i] + n2[i]*theta); 	// +1 velocity (moving right)
		f2_2[i] = 0.5 * (n2[i]*u2[i]*u2[i] - n2[i]*u2[i] + n2[i]*theta); 	// -1 velocity (moving left)
		psi2[i] = 0.0;
		mu2[i] = 0.0;
		muNonIdeal2[i] = 0.0;

		pressure[i] = 0.0;
		gradMuForce1[i] = 0.0;
		gradMuForce2[i] = 0.0;
		friction1[i] = 0.0;
		friction2[i] = 0.0;
		F1[i] = 0.0;
		F2[i] = 0.0;
	} // end for

} // end function initialize()


/**
 * @brief Function _init_ is a wrapper function for the very first program/simulation initialization.  It sets the default forcing method to the "log"
 * implementation of chemical potential forcing and the random profile for the density initialization.
 *
 * @note This function is only called once by the main() function on program start-up.
 */
void init() {
	setCollisionForcingNewChemicalPotentialGradient();
	initializeProfile = initializeRandom;
	initialize();
} // end function init()

