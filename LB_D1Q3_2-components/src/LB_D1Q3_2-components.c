/**
 * @file LB_D1Q3_2-components.c
 * @author Kent S. Ridl
 * @date 2 December 2017
 *
 * The module _LB_D1Q3_2-components.c_ contains the code to run a lattice Boltzmann simulation.  It also contains the definitions of all externs from the header
 * file _LB_D1Q3_2-components.h_.
 */


#include "LB_D1Q3_2-components.h"


//
// Define the externs declared in the .h file
//

// Directory to read/write log files and scripts
char *dataDirectory = "";

// Domain and phase information
double theoreticalDensityA1 = 0.0;					// theoretical densities and volumina/interface widths
double theoreticalDensityA2 = 0.0;
double theoreticalDensityA3 = 0.0;
double theoreticalDensityB1 = 0.0;
double theoreticalDensityB2 = 0.0;
double theoreticalDensityB3 = 0.0;
double theoreticalVolume1 = 0.0;
double theoreticalVolume2 = 0.0;
double theoreticalVolume3 = 0.0;
double theoreticalPressure = 0;
double theoreticalMuA = 0;
double theoreticalMuB = 0;
double interfaceWidth = XDIM / 20;					// default to 5% of the lattice size
double interfaceWidthForGivenKappa = XDIM / 20;
int theoreticalPhases = 0;
int initializeRandomComponents = 0;					// default both components when random init
int domain1 = XDIM / 4;								// locations and widths of 4 possible domains
int domain1Width = XDIM / 8;
int domain2 = XDIM / 2;
int domain2Width = XDIM / 8;
int domain3 = 3 * XDIM / 4;
int domain3Width = XDIM / 8;
int domain4 = XDIM;
int domain4Width = XDIM / 8;
int phase1Index = 0;
int phase2Index = 0;

// Physical properties common to minimization and LB simulations
double nA0 = 1.0;
double nB0 = 1.0;
double nAIntegrated = 0.0;  		// track total mass/densities
double nBIntegrated = 0.0;
double theta = 1./3.;
double aA = 0.1;					// VDW constants for each component, interaction
double aB = 0.1;
double aAB = 0.1;
double bA = 1./3.;
double bB = 1./3.;
double vdwInteractionFactor = 0.5;	// nu in paper
double Amp=0.01;
double tcA = 0.4;					// critical point traits
double ncA = 1.0;
double pcA = 1.0;
double tcB = 0.4;
double ncB = 1.0;
double pcB = 1.0;
double oneOverTau = 1;				// 1/tau from write-up (relaxation constant)
double tau = 1;
double g = 0;          				// gravitational acceleration term
double lambda = 1.0;				// friction (F12) coefficient

// Free energy minimization control
double minimizationParticlesStepSize = 0.01;
double minimizationStepFactor = 0.5; //0.1;
double invalidFreeEnergy = 1000.0;
int sortABinodalByIncreasingB = 1;
int sortBBinodalByIncreasingA = 1;
int threePhaseRegionExists = 0; 					// assume default is no 3-phase region exist

// 3-phase region properties and control
double threePhaseRegion[6] = {0, 0, 0, 0, 0, 0}; 	// default to no 3-phase region in a phase diagram
double rhoAThreePhaseRegion = 0.0;
double rhoBThreePhaseRegion = 0.0;
double maxArea = 0.0;
double metastableThreshold = 0.01; //0.009;
double metastableThresholdFactor = 20.0;
int setThreePhaseRegion = 1; 						// need to set once program initialized
int setLineAPoint = 1;								// divide 3-phase and binary liquid regions
int setLineBPoint = 1;

// LB simulation initialization and runtime control
double kappa = 0.1;									// interfacial free energy
double kappaFactor = 1.0;
double gammaP = 1;									// pressure coefficient for forcing rate
double gammaMu = 0.1;								// chemical potential coefficient for forcing rate
double gammaFactor = 1.0;
double endGammaFactor = 1000.0;
int useTwoOrThreePhaseInitialization = 2;			// default to 2-phase step profile
int useStepOrRandomInitialization = 1;				// default to step/tanh profile
int useTheoreticalDensities = 1; 					// default to read in theory minimization densities
int useTheoreticalVolumina = 0;						// default to equal volumina for each domain
int useTheoreticalPhase1 = 1;
int useTheoreticalPhase2 = 1;
int useTheoreticalPhase3 = 0;
int suppressTheoreticalPhase1 = 0;
int suppressTheoreticalPhase2 = 0;
int isolateFourthPhase = 0;
int overrideMinimumInterfaceWidth = 0;
int goodInterfaceFit = 1;
int fitInterfaceWidthToKappa = 0;
int determineInterfaceWidthAutomatically = 1;
int determineKappaAutomatically = 1;
int kappaFactorDetermined = 0;
int useChemicalPotentialForcingMethod = 2;
int useBoundaryConditionsPeriodic = 1; 				// default periodic BCs
int usePressurePartitionFunction = 1; 				// default to pressure derived from PF (not from construction)
int useMuVdwInteraction = 1;
int setPhaseDiagramTwoPhaseRegion = 1;
int setPhaseDiagramThreePhaseRegion = 0;
int checkTwoPhaseMetastableLBPoints = 1;			// sub-choices for simulating 3-phase region
int checkThreePhaseLBPoints = 1;
int applyMomentumCorrection = 0;

// For search tool to find parent from either minimization or LB sim
double childRhoA = 1.0;
double childRhoB = 1.0;

// To specify a single test point to minimize (instead of a whole landscape of test points)
double singlePointRhoA = 0.4; //1.0;
double singlePointRhoB = 0.4; //1.0;

// GUI control flags
int next=0;
int Pause=1;
int run = 0;
int done=0;
int Repeat=100;
int iterations;
int tmp_phase_iterations = 0;
int phase_iterations = 50000;

void (*collision)();  ///< function pointer to choose the collision method for a LB simulation


/**
 * @brief Function _calculateMass_ calculates the mass each cell over the lattice.  It also performs a numerical integration of the mass over the lattice to help verify
 * the average densities of each component stay constant throughout a simulation.
 */
void calculateMass() {
	nAIntegrated = 0;
	nBIntegrated = 0;

	// Sweep across the lattice to conserve mass and determine the pressures at each cell
	#pragma omp parallel for reduction(+:nAIntegrated,nBIntegrated)
	for (int i = 0; i < XDIM; i++) {
		n1[i] = f1_0[i] + f1_1[i] + f1_2[i]; 		// 1st component mass density for this step
		n2[i] = f2_0[i] + f2_1[i] + f2_2[i]; 		// 2st component mass density for this step
		n[i] = n1[i] + n2[i];          				// conservation of particles
		nAIntegrated += n1[i];
		nBIntegrated += n2[i];

		if (n[i] != 0) nReduced[i] = n1[i]*n2[i] / n[i];
		else nReduced[i] = 0;
	} // end for

	nAIntegrated /= XDIM;
	nBIntegrated /= XDIM;
} // end function calculateMass()


/**
 * @brief Function _calculateVelocities_ calculates the mixture's mean velocity (a hydrodynamic variable) and the velocity for each component (non-hydrodynamic variables)
 * for each cell over the lattice.
 */
void calculateVelocities() {
	#pragma omp parallel for
	for (int i = 0; i < XDIM; i++) {
		if (n1[i] != 0) u1[i] = (f1_1[i]-f1_2[i]) / n1[i];
		else u1[i] = 0.0;
		
		if (n2[i] != 0) u2[i] = (f2_1[i]-f2_2[i]) / n2[i];
		else u2[i] = 0.0;

		if (n[i] != 0) u[i] = (n1[i]*u1[i] + n2[i]*u2[i]) / n[i];	// bulk average velocity
		else u[i] = 0;
	} // end for
} // end function calculateVelocities()


/**
 * @brief Function _correctVelocities_ applies the (0.5/rho*Force) correction factor to each component's velocity.  Only the thermodynamic force from a
 * chemical potential gradient is used to correct the velocities.
 */
void correctVelocities() {
	#pragma omp parallel for
	for (int i = 0; i < XDIM; i++) {
		if (n1[i] != 0) uHat1[i] = u1[i] + 0.5/n1[i]*gradMuForce1[i]; //*F1[i]; // velocity correction needed for the forcing methods
		else uHat1[i] = 0.0;

		if (n2[i] != 0) uHat2[i] = u2[i] + 0.5/n2[i]*gradMuForce2[i]; //*F2[i];  // velocity correction needed for the forcing methods
		else uHat2[i] = 0.0;
	} // end for
} // end function correctVelocities()


/**
 * @brief Function _correctExcessMomentum_ calculates a force per particle correction term.  The correction is applied to each lattice cell by
 * weighting the correction according to the density in that cell and subtracting from the actual force for each component applied to the cell.
 * This adjustment helps to correct for velocity errors that arise due to the discrete nature of the gradients used in force calculations.
 *
 * @note The global flag _applyMomentumCorrection_ controls whether or not this correction is applied. (default is off)
 */
void correctExcessMomentum() {
	double totalParticles = 0.0, totalForce = 0.0;
	double correction = 0.0;

	#pragma omp parallel for reduction(+:totalParticles,totalForce)
	for (int i = 0; i < XDIM; i++) {
		totalParticles += n1[i] + n2[i];
		totalForce += F1[i] + F2[i];
	}

	correction = totalForce / totalParticles;

	#pragma omp parallel for
	for (int i = 0; i < XDIM; i++) {
		F1[i] -= n1[i] * correction;
		F2[i] -= n2[i] * correction;
	}

//	static double totalF1 = 0.0, totalF2 = 0.0;
//
//	// Total force per lattice site is correction applied
//	#pragma omp parallel for
//	for (int i = 0; i < XDIM; i++) {
//		F1[i] -= totalF1 / XDIM; //F1correction;
//		F2[i] -= totalF2 / XDIM; //F2correction;
//	}
//
//	// Sum up new total forces for next iteration's correction application
//	totalF1 = 0;
//	totalF2 = 0;
//	#pragma omp parallel for reduction(+:totalF1,totalF2)
//	for (int i = 0; i < XDIM; i++) {
//		totalF1 += F1[i];
//		totalF2 += F2[i];
//	}
} // end function correctExcessMomentum()


/**
 * @brief Function _calculateLBPressure_ calculates pressure for the full 2-component mixture for each cell over the lattice.
 *
 * The pressure of the mixture is calculated in this function.  The pressure includes gradient terms, and the components are coupled together and cannot be
 * separated into meaningful partial pressures.
 *
 * @note The global parameter _gammaP_ is a "filter" to modulate the pressure that is applied each time step and aid
 * in stabilizing the simulation.
 * @note The global parameter _usePressurePartitionFunction_ allows the user to select from two pressure formulations: one is a constructed pressure tensor
 * and one is derived from a partition.  This feature was used in development and is preserved for "gee whiz" purposes.
 */
void calculateLBPressure() {
	#pragma omp parallel for
	for (int i = 0; i < XDIM; i++) {
		if (usePressurePartitionFunction) {
			pressure[i] = n[i]*theta*(1 + (bA*n1[i]+bB*n2[i])/(1.-bA*n1[i]-bB*n2[i])) - aA*n1[i]*n1[i] - 2*aAB*n1[i]*n2[i] - aB*n2[i]*n2[i];
		}
		else {
			pressure[i] = n1[i]*theta/(1.-bA*n1[i]-bB*n2[i]) + n2[i]*theta/(1.-bA*n1[i]-bB*n2[i]) - aA*n1[i]*n1[i] - 2*aAB*n1[i]*n2[i] - aB*n2[i]*n2[i];
		}

		// Gradient corrections for each single component, including self-interactions
		pressure[i] += -kappa*( n1[i]*laplace(n1,i) + 0.5*gradient(n1,i)*gradient(n1,i) + n2[i]*laplace(n2,i) + 0.5*gradient(n2,i)*gradient(n2,i) );
		pressure[i] += kappa*( gradient(n1,i)*gradient(n1,i) + gradient(n2,i)*gradient(n2,i) );

		// Gradient corrections for the cross terms (i.e. cross interactions)
		pressure[i] += -kappa*( n1[i]*laplace(n2,i) + n2[i]*laplace(n1,i) + gradient(n1,i)*gradient(n2,i) );
		pressure[i] += kappa*( 2.*gradient(n1,i)*gradient(n2,i) );

		pressure[i] *= gammaP;
	}
} // end function calculateLBPressure()


/**
 * @brief Function _calculateLBChemicalPotentials_ calculates the chemical potential for each component for each cell over the lattice.
 *
 * This function calculates the chemical potentials including gradient terms that are the basis for the LB forcing terms.  It calculates both the full
 * chemical potentials and the non-ideal parts of the chemical potentials.
 */
void calculateLBChemicalPotentials() {
	#pragma omp parallel for
	for (int i = 0; i < XDIM; i++) {
		mu1[i] = gammaMu * ( theta*log(n1[i]/(1.-bA*n1[i]-bB*n2[i])) + theta*bA*(n1[i]+n2[i])/(1.-bA*n1[i]-bB*n2[i])
				- 2*aA*n1[i] - 2*aAB*n2[i] - kappa*laplace(n1,i) - (useMuVdwInteraction ? vdwInteractionFactor : 1)*kappa*laplace(n2,i) );
		mu2[i] = gammaMu * ( theta*log(n2[i]/(1.-bA*n1[i]-bB*n2[i])) + theta*bB*(n1[i]+n2[i])/(1.-bA*n1[i]-bB*n2[i])
				- 2*aB*n2[i] - 2*aAB*n1[i] - kappa*laplace(n2,i) - (useMuVdwInteraction ? vdwInteractionFactor : 1)*kappa*laplace(n1,i) );

		muNonIdeal1[i] = mu1[i] - theta*log(n1[i]);
		muNonIdeal2[i] = mu2[i] - theta*log(n2[i]);
	}
} // end function calculateLBChemicalPotentials()


/**
 * @brief Function _calculateFriction_ calculates the average momentum transfer (friction force) between components A and B.  It then adds the
 * friction force to the conservative force from a chemical potential gradient to give the total force for each lattice cell.
 */
void calculateFriction() {
	#pragma omp parallel for
	for (int i = 0; i < XDIM; i++) {
		friction1[i] = nReduced[i] * (uHat1[i]-uHat2[i]); // friction forces on 1st component
		friction2[i] = nReduced[i] * (uHat2[i]-uHat1[i]); // friction forces on 2nd component

		// Final forcing for each component including friction
		F1[i] = gradMuForce1[i] - lambda*friction1[i];
		F2[i] = gradMuForce2[i] - lambda*friction2[i];
	} // end for
} // end function calculateFriction()


/**
 * @brief Function _collisionForcingNewChemicalPotentialGradient_ uses a chemical potential gradient forcing method for the lattice Boltzmann collision step.
 *
 * This function uses the chemical potentials of each component as the basis of the forcing in the LB collision step of the algorithm.  There are 2 chemical
 * potential gradient formulations:
 * 1. The "nid" formulation is a gradient of the non-ideal portion of the chemical potential.
 * 2. The "log" forumlation is a gradient of the full chemical potential minus the gradient of an ideal pressure.
 * The global parameter _useChemicalPotentialForcingMethod_ is used to choose the formulation (1=nid, 2=log).
 * This function also applies a 4th order force correction to help insure thermodynamic consistency of the LB simulation (see reference).
 *
 * @see A. J. Wagner, Phys. Rev. E 74, 056703 (2006).
 */
void collisionForcingNewChemicalPotentialGradient() {
	static int lastMuForcingMethod = 0; // default to the first case below

	iterations++;

	calculateMass();
	calculateVelocities();
	calculateLBPressure();
	calculateLBChemicalPotentials();

	//
	// Forcing derived from chemical potential gradients... a la Gibbs-Duhem (sum over components of rho*gradMu equals pressure gradient)
	//
	if (useChemicalPotentialForcingMethod != lastMuForcingMethod) { // print out a statement when the forcing method changes
		lastMuForcingMethod = useChemicalPotentialForcingMethod;
		switch (useChemicalPotentialForcingMethod) {
		case 1: // grad non-ideal mu
			printf("\"nid\" grad-Mu forcing: -nx*grad(MuNidx)\n");
			break;
		case 2: // gradient of mu minus ideal pressure
			printf("\"log\" grad-Mu forcing: -nx*grad(Mux)-theta*grad(nx)\n");
			break;
		default:
			useChemicalPotentialForcingMethod = 2;
			lastMuForcingMethod = 2;
			printf("Invalid Selection!  Defaulting to \"log\" grad-Mu forcing method 2: -nx*grad(Mux)-theta*grad(nx)\n");
		} // end switch
	} // end if
	switch (useChemicalPotentialForcingMethod) {
	case 1: // grad non-ideal mu (nid)
		#pragma omp parallel for
		for (int i = 0; i < XDIM; i++) {
			gradMuForce1[i] = -1.*n1[i]*gradient(muNonIdeal1,i) + n1[i]*g;
			gradMuForce2[i] = -1.*n2[i]*gradient(muNonIdeal2,i) + n2[i]*g;
		}
		break;
	case 2: // gradient of mu minus ideal pressure/chemical potential gradient (log)
		#pragma omp parallel for
		for (int i = 0; i < XDIM; i++) {
			gradMuForce1[i] = -1. * ( n1[i]*gradient(mu1,i)-theta*gradient(n1,i) ) + n1[i]*g;
			gradMuForce2[i] = -1. * ( n2[i]*gradient(mu2,i)-theta*gradient(n2,i) ) + n2[i]*g;
		}
		break;
	//default: // do nothing...
	} // end switch

	correctVelocities(); // velocities are corrected by the conservative rho*gradMu force only (friction not included)
	calculateFriction();
	if (applyMomentumCorrection) correctExcessMomentum();

	#pragma omp parallel for
	for (int i = 0; i < XDIM; i++) {

		// Correction to the equilibrium distribution that alters the actual forcing
		if (n1[i] !=0) {
			psi1[i] = -oneOverTau * ( (tau-0.25)*F1[i]*F1[i]/n1[i] + (1./12.)*laplace(n1,i) );	// subtract psi, so minus sign relative to paper
		}
		else psi1[i] = 0;

		// Calculate particle densities at current lattice spot with forcing included
		f1_0[i] += oneOverTau * ( (n1[i] - n1[i]*theta - n1[i]*u1[i]*u1[i]) - f1_0[i] ) - ( 2.*F1[i]*u1[i] - psi1[i] );
		f1_1[i] += oneOverTau * ( 0.5*(n1[i]*u1[i]*u1[i]+n1[i]*u1[i]+n1[i]*theta)-f1_1[i] ) - ( -F1[i]*u1[i] - 0.5*F1[i] + 0.5*psi1[i] );
		f1_2[i] += oneOverTau * ( 0.5*(n1[i]*u1[i]*u1[i]-n1[i]*u1[i]+n1[i]*theta)-f1_2[i] ) - ( -F1[i]*u1[i] + 0.5*F1[i] + 0.5*psi1[i] );

		// Correction to the equilibrium distribution that alters the actual PGF to pressure is constant in equilibirum
		if (n2[i] != 0) {
			psi2[i] = -oneOverTau * ( (tau-0.25)*F2[i]*F2[i]/n2[i] + (1./12.)*laplace(n2,i) );	// subtract psi, so minus sign relative to paper
		}
		else psi2[i] = 0;

		f2_0[i] += oneOverTau * ( (n2[i] - n2[i]*theta - n2[i]*u2[i]*u2[i]) - f2_0[i] ) - ( 2.*F2[i]*u2[i] - psi2[i] );
		f2_1[i] += oneOverTau * ( 0.5*(n2[i]*u2[i]*u2[i] + n2[i]*u2[i] + n2[i]*theta) - f2_1[i] ) - ( -F2[i]*u2[i] - 0.5*F2[i] + 0.5*psi2[i] );
		f2_2[i] += oneOverTau * ( 0.5*(n2[i]*u2[i]*u2[i] - n2[i]*u2[i] + n2[i]*theta) - f2_2[i] ) - ( -F2[i]*u2[i] + 0.5*F2[i] + 0.5*psi2[i] );
	} // end for
} // end function collisionForcingNewChemicalPotentialGradient()


/**
 * @brief Function _setCollisionForcingNewChemicalPotentialGradient_ sets the forcing method used by the function _collision_ to be the gradient of a
 * chemical potential.
 */
void setCollisionForcingNewChemicalPotentialGradient() {
	collision = collisionForcingNewChemicalPotentialGradient;
	printf("Using gradMu forcing - new...\n");
} // end function setCollisionForcingNewChemicalPotentialGradient()


/**
 * @brief Function _streaming_ is the streaming step of the lattice Boltzmann algorithm.
 *
 * The streaming step of the LB algorithm shifts the particle distributions for each component to "move" the particles.  This function is specific to the
 * 1-D lattice and defaults to use periodic boundary conditions.  If desired, the global flag _useBoundaryConditionsPeriodic_ may be toggled to change the
 * ends of the lattice to impose solid end points with bounce-back boundaries.
 */
void streaming() {
	double tmp;

	/* Original wrap-around end points */
	tmp=f1_1[XDIM-1];                                   // save right end point
	memmove(&f1_1[1],&f1_1[0],(XDIM-1)*sizeof(double)); // shift all cells +1
	f1_1[0]=tmp;                                        // rotate former end to first lattice cell
	tmp=f1_2[0];                                        // save left end point
	memmove(&f1_2[0],&f1_2[1],(XDIM-1)*sizeof(double)); // shift all cells -1
	f1_2[XDIM-1]=tmp;                                   // rotate former first lattice cell to end

	tmp=f2_1[XDIM-1];                                   // save right end point
	memmove(&f2_1[1],&f2_1[0],(XDIM-1)*sizeof(double)); // shift all cells +1
	f2_1[0]=tmp;                                        // rotate former end to first lattice cell
	tmp=f2_2[0];                                        // save left end point
	memmove(&f2_2[0],&f2_2[1],(XDIM-1)*sizeof(double)); // shift all cells -1
	f2_2[XDIM-1]=tmp;                                   // rotate former first lattice cell to end

	// Walls at the end points; bounce from lattice origin (0)
	if (!useBoundaryConditionsPeriodic) {
		tmp = f1_1[0];
		f1_1[0] = f1_2[XDIM-1];
		f1_2[XDIM-1]=tmp;

		tmp = f2_1[0];
		f2_1[0] = f2_2[XDIM-1];
		f2_2[XDIM-1]=tmp;
	}
} // end function streaming()


/**
 * @brief Function _iteration_ is the governing function for each iteration of the lattice Boltzmann algorithm, executing the collision and streaming steps
 * in succession.
 */
void iteration(){

	// Need to reset the critical and VDW constants each iteration
	// Keeps them all in sync if one is changed during a simulation
	pcA = 3.*tcA/8.; // determine pcA
	ncA = pcA / ((3./8.)*tcA); // fix ncA to be 1
	pcB = (3./8.)*tcB*ncB; // determine pcB

	aA = (27./64.)*(tcA*tcA/pcA);
	aB = (27./64.)*(tcB*tcB/pcB);
	aAB = sqrt(aA*aB) * vdwInteractionFactor;
	bA = tcA/(8.*pcA);
	bB = tcB/(8.*pcB);

	collision();
	streaming();

	// When running manual simulations, automatically stop the simulation if this iteration blows up
	if ((!Pause || run) && (!stableSimulation(n1) || !stableSimulation(n2))) {
		Pause = 1;
		run = 0;
	}
} // end function iteration()


//void calculateFreeEnergyLattice(){
//	int i = 0;
//	double excludedVolume = 0;
//
//	// Sum over the lattice to determine total free energy given free energy densities at each lattice site
////	#pragma omp parallel for private(excludedVolume)
//	for (i = 0; i < XDIM; i++) {
//		excludedVolume = n1[i]*bA + n2[i]*bB;
//		if ((n1[i] < 0) || (n2[i] < 0) || (excludedVolume > volumeTotal)) {
//			freeEnergy[i] = invalidFreeEnergy;
//		}
//		else if ((n1[i] == 0) && (n2[i] == 0)) {
//			freeEnergy[i] = 0;
//		}
//		else if (n1[i] == 0) {
//			freeEnergy[i] = n2[i]*theta*log(n2[i]/(volumeTotal-excludedVolume)) - aB*n2[i]*n2[i]/volumeTotal - theta*n2[i];
//		}
//		else if (n2[i] == 0) {
//			freeEnergy[i] = n1[i]*theta*log(n1[i]/(volumeTotal-excludedVolume)) - aA*n1[i]*n1[i]/volumeTotal - theta*n1[i];
//		}
//		else {
//			freeEnergy[i] = n1[i]*theta*log(n1[i]/(volumeTotal-excludedVolume)) + n2[i]*theta*log(n2[i]/(volumeTotal-excludedVolume))
//			- aA*n1[i]*n1[i]/volumeTotal - 2*aAB*n2[i]*n1[i]/volumeTotal - aB*n2[i]*n2[i]/volumeTotal - theta*(n1[i]+n2[i]);
//		}
//	} // end for
//} // end function calculateFreeEnergy()
//
//
//void correctPressure() {
//	int i = 0;
//
//	for (i = 0; i < XDIM; i++) {
//		pressureCorrection[i] = -(tau-0.25)*(F1[i]*F1[i]+F2[i]*F2[i])/n[i] + 0.25*(n1[i]*gradient(F1,i)+n2[i]*gradient(F2,i))/n[i] - (1./12.)*(n1[i]*laplace(n1,i)+n2[i]*laplace(n2,i))/n[i]; // -(tau-0.25)*(F1[i]*F1[i]+F2[i]*F2[i])/n[i] + 0.25*(n1[i]*F1[i]*F1[i]+n2[i]*F2[i]*F2[i])/n[i]
//		correctedPressure[i] = pressure[i] - pressureCorrection[i];
//	}
//}
