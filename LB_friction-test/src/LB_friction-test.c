/**
 * @file LB_D1Q3_2-components.c
 * @author Kent S. Ridl
 * @date 2 December 2017
 *
 * The module _LB_D1Q3_2-components.c_ contains the code to run a lattice Boltzmann simulation.  It also contains the definitions of all externs from the header
 * file _LB_D1Q3_2-components.h_.
 */


#include "LB_friction-test.h"


//
// Define the externs declared in the .h file
//

int initializeRandomComponents = 0;					// default both components when random init

// Physical properties common to minimization and LB simulations
double nA0 = 1.0;
double nB0 = 1.0;
double uA0 = 0.0;
double uB0 = 0.0;
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
double lambda = 0.01;				// friction (F12) coefficient

// LB simulation initialization and runtime control
double kappa = 0.1;									// interfacial free energy
double gammaP = 1;									// pressure coefficient for forcing rate
double gammaMu = 0.1;								// chemical potential coefficient for forcing rate
int useChemicalPotentialForcingMethod = 2;
int useBoundaryConditionsPeriodic = 1; 				// default periodic BCs
int usePressurePartitionFunction = 1; 				// default to pressure derived from PF (not from construction)

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
	#pragma omp parallel for
	for (int i = 0; i < XDIM; i++) {
		n1[i] = f1_0[i] + f1_1[i] + f1_2[i]; 		// 1st component mass density for this step
		n2[i] = f2_0[i] + f2_1[i] + f2_2[i]; 		// 2st component mass density for this step
		n[i] = n1[i] + n2[i];          				// conservation of particles

		if (n[i] != 0) nReduced[i] = n1[i]*n2[i] / n[i];
		else nReduced[i] = 0;
	} // end for
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
 * @brief Function _correctVelocities_ applies the (0.5/rho*Force) correction factor to each component's velocity.
 */
void correctVelocities() {
	#pragma omp parallel for
	for (int i = 0; i < XDIM; i++) {
		if (n1[i] != 0) uHat1[i] = u1[i] + 0.5/n1[i]*gradMuForce1[i]; //*F1[i]; // velocity correction needed for the forcing methods
		else uHat1[i] = 0.0;

		if (n2[i] != 0) uHat2[i] = u2[i] + 0.5/n2[i]*gradMuForce2[i]; //F2[i];  // velocity correction needed for the forcing methods
		else uHat2[i] = 0.0;
//		uHat1[i] = u1[i];
//		uHat2[i] = u2[i];
	} // end for
} // end function correctVelocities()


///**
// * @brief Function _calculatePressureCoupled_ calculates pressure for the full 2-component mixture for each cell over the lattice.
// *
// * The pressure of the mixture is calculated in this function.  The pressure includes gradient terms, and the components are coupled together and cannot be
// * separated into meaningful partial pressures.
// *
// * @note The global parameter _gammaP_ is a "filter" to modulate the pressure that is applied each time step and aid
// * in stabilizing the simulation.
// * @note The global parameter _usePressurePartitionFunction_ allows the user to select from two pressure formulations: one is a constructed pressure tensor
// * and one is derived from a partition.  This feature was used in development and is preserved for "gee whiz" purposes.
// */
//void calculatePressureCoupled() {
//	#pragma omp parallel for
//	for (int i = 0; i < XDIM; i++) {
//		if (usePressurePartitionFunction) {
//			pressure[i] = n[i]*theta*(1 + (bA*n1[i]+bB*n2[i])/(1.-bA*n1[i]-bB*n2[i])) - aA*n1[i]*n1[i] - 2*aAB*n1[i]*n2[i] - aB*n2[i]*n2[i];
//		}
//		else {
//			pressure[i] = n1[i]*theta/(1.-bA*n1[i]-bB*n2[i]) + n2[i]*theta/(1.-bA*n1[i]-bB*n2[i]) - aA*n1[i]*n1[i] - 2*aAB*n1[i]*n2[i] - aB*n2[i]*n2[i];
//		}
//
//		// Gradient corrections for each single component, including self-interactions
//		pressure[i] += -kappa*( n1[i]*laplace(n1,i) + 0.5*gradient(n1,i)*gradient(n1,i) + n2[i]*laplace(n2,i) + 0.5*gradient(n2,i)*gradient(n2,i) );
//		pressure[i] += kappa*( gradient(n1,i)*gradient(n1,i) + gradient(n2,i)*gradient(n2,i) );
//
//		// Gradient corrections for the cross terms (i.e. cross interactions)
//		pressure[i] += -kappa*( n1[i]*laplace(n2,i) + n2[i]*laplace(n1,i) + gradient(n1,i)*gradient(n2,i) );
//		pressure[i] += kappa*( 2.*gradient(n1,i)*gradient(n2,i) );
//
//		pressure[i] *= gammaP;
//	}
//} // end function calculatePressureCoupled()


/**
 * @brief Function _calculateChemicalPotentials_ calculates the chemical potential for each component for each cell over the lattice.
 *
 * This function calculates the chemical potentials including gradient terms that are the basis for the LB forcing terms.  It calculates both the full
 * chemical potentials and the non-ideal parts of the chemical potentials.
 */
void calculateChemicalPotentials() {
	#pragma omp parallel for
	for (int i = 0; i < XDIM; i++) {
//		mu1[i] = gammaMu * ( theta*log(n1[i]/(1.0-bA*n1[i])) + theta*bA*n1[i]/(1.0-bA*n1[i]) - 2.0*aA*n1[i] - kappa*laplace(n1,i) );
//		mu2[i] = gammaMu * ( theta*log(n2[i]/(1.0-bB*n2[i])) + theta*bB*n2[i]/(1.0-bB*n2[i]) - 2.0*aB*n2[i] - kappa*laplace(n2,i) );

		mu1[i] = gammaMu * ( theta*log(n1[i]/(1.-bA*n1[i]-bB*n2[i])) + theta*bA*(n1[i]+n2[i])/(1.-bA*n1[i]-bB*n2[i])
				- 2*aA*n1[i] - 2*aAB*n2[i] - kappa*laplace(n1,i) - 0.5*kappa*laplace(n2,i) );
		mu2[i] = gammaMu * ( theta*log(n2[i]/(1.-bA*n1[i]-bB*n2[i])) + theta*bB*(n1[i]+n2[i])/(1.-bA*n1[i]-bB*n2[i])
				- 2*aB*n2[i] - 2*aAB*n1[i] - kappa*laplace(n2,i) - 0.5*kappa*laplace(n1,i) );

		muNonIdeal1[i] = mu1[i] - theta*log(n1[i]);
		muNonIdeal2[i] = mu2[i] - theta*log(n2[i]);
	}
} // end function calculateChemicalPotentials()


/**
 * @brief Function _calculateFriction_ calculates the average momentum transfer (friction force) between components A and B.  It then adds the
 * friction force to the conservative force from a chemical potential gradient to give the total force for each lattice cell.
 */
void calculateFriction() {

	// TODO: Verify this by systematically modeling mixture of ideal fluids
	// 1. model where friction is the only force; initialize both fluids with small mean velocity and expect them to slow down
	// 2. model where friction and conservative mu force both exist; should see phase separation where the phases move but slow down
	// Should be able to apply very large values of lambda and stop all fluid motion in a single time step
	#pragma omp parallel for
	for (int i = 0; i < XDIM; i++) {
//		// 1st run - order: friction, correct velocities
//		// lambda=0.0001: uHat barely 10^-4 at ~50K iterations when error starts increasing... unstable ~225K iterations
//		//
//		// 2nd run - same state as the last stable Amp=0.01, lambda=0.01 below
//		// Amp=0.01
//		// lambda=0.01: unstable by 1600 iterations - uHat, friction, force magnitudes increase every step
//		// lambda=0.001: everything fluctuates around same initial magnitudes for ~8500 iterations, then slowly increases...
//		//				 ... first uHat noise by 14,500 iterations, unstable by 16,000 iterations
//		// lambda=0.0001: same as 1st run
//		friction1[i] = nReduced[i] * (uHat1[i]-uHat2[i]); // friction forces on 1st component
//		friction2[i] = nReduced[i] * (uHat2[i]-uHat1[i]); // friction forces on 2nd component
//		F1[i] = gradMuForce1[i] + lambda*friction1[i];
//		F2[i] = gradMuForce2[i] + lambda*friction2[i];

		//////////////////////////////////////////////////
		//					Best Run					//
		//////////////////////////////////////////////////

		// 1st run - order: friction, correct velocities
		// lambda=0.0001: uHat 10^-6, ~100K iterations error+noise increases until unstable ~300K iterations
		//
		// 2nd run - work down towards same state as stable Amp=0.01, lambda=0.01 below
		// Amp=0.01
		// lambda=0.0001: same as 1st run
		//
		// subsequent runs: - troubleshot with Alexander
		friction1[i] = nReduced[i] * (uHat1[i]-uHat2[i]); // friction forces on 1st component
		friction2[i] = nReduced[i] * (uHat2[i]-uHat1[i]); // friction forces on 2nd component
		F1[i] = gradMuForce1[i] - lambda*friction1[i];
		F2[i] = gradMuForce2[i] - lambda*friction2[i];

//		// 1st run - order: friction, correct velocities
//		// lambda=0.0001: both uHat 10^-15, can't recall seeing acceleration
//		// lambda=0.001: ~500K iterations, both uHat 10^-11/-12 start accelerating - direction (can't see at 10^-10 scale)
//		//
//		// 2nd run - no code changes
//		// lambda=0.0001: uHat1 10^-15, uHat2 noise 10^-13
//		//				  init-rand: both uHat 10^-15 (difference 10^-14) and accelerating - direction
//		//				  init-rand: both uHat 10^-13 noise and accelerating + direction (?!!!)
//		//				  init-rand: both uHat 10^-15, no acceleration
//		//				  reinit: both uHat 10^-15 (difference 10^-14) and accelerating - direction
//		//				  reinit: uHat1 10^-15, uHat2 noise 10^-13, no acceleration
//		//				  reinit: uHat1 10^-15, uHat2 noise 10^-13, no acceleration
//		//				  reinit: both uHat 10^-14, - accel, 10^-13 noise, + accel, noise dis/reappear, accel reverse (cycles)
//		//				  reinit: both uHat 10^-15, - accel, ~1.85M iterations accel stop/start with error+noise 10^-14/10^-13
//		//
//		// 3rd run - init nReduced with actual formula instead of setting = n[i]
//		// lambda=0.0001: same variable behavior (as expected)... Amp=0, nothing happens (as expected)
//		//				  cyclic behavior must be due to the variable random noise with each initialization
//		// Amp=1e-6
//		// lambda=0.001: consistent 10^-14 accuracy, occasional 10^-15, small accelerations (basically behaves just like lambda=0.0001)
//		// lambda=0.01: consistent uHat 10^-13/-12 after ~300K iterations, occasional 10^-14, noise in and out, small accelerations
//		// Amp=0.01
//		// lambda=0.01: uHat only driven to 10^-10 diff (both relatively noise-free to 10^-13), takes longer to for uHat to decay (makes sense)...
//		//				... friction still =/opposite but 10^-10 magnitude, force eventually driven below machine accuracy (but never constant??)...
//		//				... imparted large - acceleration and occasionaly appears to be adding total mass ~10^-12 scale (maybe a neglibible amount?)...
//		//				... first time behavior appears to be consistent from init-to-init (probably just from larger magnitude scales, though)
//		// ** Never unstable to this point **
//		// lambda=1.0: unstable by 1500 iterations
//		//
//		// 4th run - explicitly take out the gradMuForce terms (even though they are initalized to zero)
//		// lambda=0.01: uHat driven to 10^-11/-12 diff (friction corresponds), large + accel this time (but magnitude was - this time... accel towards zero)
//		// 				virtually the same as the end of 3rd run
//		//
//		// 5th run - put gradMuForce back in
//		// lambda=0.1: unstable by 34,000 iterations
//		//
//		// 6th run - swap order: correct velocities, friction
//		// lambda=0.01: same behavior as 3rd/4th runs
//		// lambda=0.001: same behavior as 1st run
//		friction1[i] = nReduced[i] * (uHat1[i]-uHat2[i]); // friction forces on 1st component
//		friction2[i] = nReduced[i] * (uHat2[i]-uHat1[i]); // friction forces on 2nd component
//		F1[i] = gradMuForce1[i] + lambda*friction1[i];
//		F2[i] = gradMuForce2[i] - lambda*friction2[i];

//		// lambda=0.0001: uHat1 noise 10^-13, uHat2 10^-15... still stable +3M iterations
//		// lambda=0.001: ~500K iterations, both uHat 10^-11/-12 start accelerating + direction
//		friction1[i] = nReduced[i] * (uHat1[i]-uHat2[i]); // friction forces on 1st component
//		friction2[i] = nReduced[i] * (uHat2[i]-uHat1[i]); // friction forces on 2nd component
//		F1[i] = gradMuForce1[i] - lambda*friction1[i];
//		F2[i] = gradMuForce2[i] + lambda*friction2[i];

//		// lambda = 0.0001: ~150K iterations, uHat 10^-6, error+noise starts building until simulation unstable ~300K iterations
//		friction1[i] = nReduced[i] * (uHat1[i]-uHat2[i]); // friction forces on 1st component
//		F1[i] = gradMuForce1[i] - lambda*friction1[i];
//		F2[i] = gradMuForce2[i] + lambda*friction1[i];

//		// lambda = 0.0001: ~150K iterations, uHat 10^-4, error starts building until simulation unstable
//		friction1[i] = nReduced[i] * (uHat1[i]-uHat2[i]); // friction forces on 1st component
//		F1[i] = gradMuForce1[i] + lambda*friction1[i];
//		F2[i] = gradMuForce2[i] - lambda*friction1[i];

//		// lambda = 0.0001: ~150K iterations, uHat 10^-4, error starts building until simulation unstable
//		friction2[i] = nReduced[i] * (uHat2[i]-uHat1[i]); // friction forces on 1st component
//		F1[i] = gradMuForce1[i] - lambda*friction2[i];
//		F2[i] = gradMuForce2[i] + lambda*friction2[i];

//		// lambda = 0.0001: ~150K iterations, uHat 10^-6, error+noise starts building until simulation unstable ~300K iterations
//		friction2[i] = nReduced[i] * (uHat2[i]-uHat1[i]);
//		F1[i] = gradMuForce1[i] + lambda*friction2[i];
//		F2[i] = gradMuForce2[i] - lambda*friction2[i];
	} // end for
} // end function calculateFriction()


/**
 * @brief Function _correctExcessMomentum_ calculates an average force applied to each lattice cell for an iteration.  On the subsequent iteration, that
 * average is subtracted from the actual force applied to the cell.  This adjustment helps to correct for velocity errors that arise due to the discrete
 * nature of the gradients used in force calculations.
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
} // end function correctExcessMomentum()


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
	calculateChemicalPotentials();

	// Forcing derived from chemical potential gradients... a la Gibbs-Duhem (sum of both equals pressure gradient)
	if (useChemicalPotentialForcingMethod != lastMuForcingMethod) { // print out a statement when the forcing method changes
		lastMuForcingMethod = useChemicalPotentialForcingMethod;
		switch (useChemicalPotentialForcingMethod) {
		case 1: // grad non-ideal mu
			printf("\"nid\" gradMu forcing: -nx*grad(MuNidx)\n");
			break;
		case 2: // gradient of mu minus ideal pressure
			printf("\"log\" gradMu forcing: -nx*grad(Mux)-theta*grad(nx)\n");
			break;
		default:
			useChemicalPotentialForcingMethod = 2;
			lastMuForcingMethod = 2;
			printf("Invalid Selection!  Defaulting to \"log\" gradMu forcing method 2: -nx*grad(Mux)-theta*grad(nx)\n");
		} // end switch
	} // end if
	switch (useChemicalPotentialForcingMethod) {
	case 1: // grad non-ideal mu (nid)
		#pragma omp parallel for
		for (int i = 0; i < XDIM; i++) {
			gradMuForce1[i] = -1.*n1[i]*gradient(muNonIdeal1,i);
			gradMuForce2[i] = -1.*n2[i]*gradient(muNonIdeal2,i);
		}
		break;
	case 2: // gradient of mu minus ideal pressure/chemical potential gradient (theta rho grad-log rho = theta grad-rho)
		#pragma omp parallel for
		for (int i = 0; i < XDIM; i++) {
			gradMuForce1[i] = -1. * ( n1[i]*gradient(mu1,i)-theta*gradient(n1,i) );
			gradMuForce2[i] = -1. * ( n2[i]*gradient(mu2,i)-theta*gradient(n2,i) );
		}
		break;
	//default: // do nothing...
	} // end switch

	correctVelocities(); // TODO: should this be before the friction calculation?
	calculateFriction(); // u-hat from previous step used to calculate friction for this step
//	calculatePressureCoupled();
	correctExcessMomentum();

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

} // end function iteration()


