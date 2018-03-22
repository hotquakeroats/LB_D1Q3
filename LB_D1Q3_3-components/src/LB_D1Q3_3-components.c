/*
 * LB_D1Q3_2components.c
 *
 *  Created on: Apr 22, 2016
 *      Author: Kent Ridl
 */


#include "LB_D1Q3_3-components.h"


//
// Define the externs declared in the .h file
//

// Initial densities of each component
double theoreticalDensityA1 = 1.77;
double theoreticalDensityA2 = 0.09;
double theoreticalDensityA3 = 0.09;
double theoreticalDensityA4 = 0.45;
double theoreticalDensityB1 = 0.09;
double theoreticalDensityB2 = 1.77;
double theoreticalDensityB3 = 0.09;
double theoreticalDensityB4 = 0.45;
double theoreticalDensityC1 = 0.09;
double theoreticalDensityC2 = 0.09;
double theoreticalDensityC3 = 1.77;
double theoreticalDensityC4 = 0.45;
int interface1 = (int)(1*PHASE_WIDTH*XDIM); // 50
int interface2 = (int)(2*PHASE_WIDTH*XDIM); // 100
int interface3 = (int)(3*PHASE_WIDTH*XDIM); // 150
int interface4 = (int)(4*PHASE_WIDTH*XDIM); // 200
int interface5 = (int)(5*PHASE_WIDTH*XDIM); // 250
int interface6 = XDIM-1; 					// 299

// Total amount of each component
double volumeTotal = 1.0;
double totalForceA = 0;
double totalForceB = 0;
double totalForceC = 0;
double correctionForceA = 0;
double correctionForceB = 0;
double correctionForceC = 0;

// Evaporation constants
double dn=0; 				// delta n for evaporation
int epos=0;  				// lattice position at which evaporation occurs

// Critical point traits
double tcA = 0.4;
double ncA = 1.0;
double pcA = 1.0;
double tcB = 0.4;
double ncB = 1.0;
double pcB = 1.0;
double tcC = 0.4;
double ncC = 1.0;
double pcC = 1.0;

double pressureFilterNeighbors = 3;
double muFilterNeighbors = 12;
double theoreticalPressure = 0;
double theoreticalMuA = 0;
double theoreticalMuB = 0;

// Equations of motion constants
double T0=0.33333333;  		// initial theta from write-up
double theta = 1./3.;
double nA0 = 1.0;
double nB0 = 1.0;
double nC0 = 1.0;
double Amp=0.01;
double oneOverTau = 1;		// 1/tau from write-up (relaxation constant)
double tau = 1;
double g = 0;          		// gravitational acceleration term
double lambda = 1e-5;     		// friction (F12) coefficient
double gammaP = 1;			// pressure coefficient for forcing rate
double gammaMu = 0.001; //0.001;		// chemical potential coefficient for forcing rate
double kappa = 0.5;

double aA = 0.1;		// VDW constants for each component, interaction
double aB = 0.1;
double aC = 0.1;
double aAB = 0.1;
double aAC = 0.1;
double aBC = 0.1;
double bA = 1./3.;
double bB = 1./3.;
double bC = 1./3.;
double b = 1./3.;		// standard VDW sizing; used in loop control
double vdwInteractionFactor = 0.35; // nu in paper
double interfaceWidth = 0.05;

int useChemicalPotentialForcingMethod = 2;
int lnExplosion = 0;
int useBoundaryConditionsPeriodic = 1; 		// default periodic BCs
int useMuVdwInteraction = 1;
int momentumCorrectionOn = 1;
int phase_iterations = 1000000;
int total_phase_iterations = 0;

// GUI control flags
int ulreq=0;
int ugreq=0;
int pgreq=0;
int tgreq=0;
int muAgreq = 0;
int next=0;
int Pause=1;
int run = 0;
int done=0;
int Repeat=100;
int iterations;
int collectData = 0;


void calculateMassAndVelocities() {
	int i = 0;

	// Sweep across the lattice to conserve mass and determine the pressures at each cell
	for (i = 0; i < XDIM; i++) {

		nA[i] = fA_0[i] + fA_1[i] + fA_2[i]; 				// 1st component mass density for this step
		nB[i] = fB_0[i] + fB_1[i] + fB_2[i]; 				// 2st component mass density for this step
		nC[i] = fC_0[i] + fC_1[i] + fC_2[i];
		n[i] = nA[i] + nB[i] + nC[i];  						// conservation of particles

		if (nA[i] != 0) {
			uA[i] = (fA_1[i]-fA_2[i]) / nA[i];
			
			// Correction to the display of the mean fluid velocity
			// Needed for the forcing methods (pressure method is good with the above)
			uHatA[i] = uA[i] + 0.5/nA[i]*forceA[i];
		}
		else {
			uA[i] = 0;
		}
		
		if (nB[i] != 0) {
			uB[i] = (fB_1[i]-fB_2[i]) / nB[i];

			// Correction to the display of the mean fluid velocity
			// Needed for the forcing methods (pressure method is good with the above)
			uHatB[i] = uB[i] + 0.5/nB[i]*forceB[i];
		}
		else {
			uB[i] = 0;
		}

		if (nC[i] != 0) {
			uC[i] = (fC_1[i]-fC_2[i]) / nC[i];

			// Correction to the display of the mean fluid velocity
			// Needed for the forcing methods (pressure method is good with the above)
			uHatC[i] = uC[i] + 0.5/nC[i]*forceC[i];
		}
		else {
			uC[i] = 0;
		}

		if (n[i] != 0) {
			u[i] = (nA[i]*uA[i] + nB[i]*uB[i] + nC[i]*uC[i]) / n[i];			// bulk average velocity
		}
		else {
			u[i] = 0;
		}
	} // end for
} // end function calculateMassAndVelocities()


void correctExcessMomentum() {
	// Total force per lattice site is correction applied
	correctionForceA = totalForceA / XDIM;
	correctionForceB = totalForceB / XDIM;
	correctionForceC = totalForceC / XDIM;
	for (int i = 0; i < XDIM; i++) {
		forceA[i] -= correctionForceA;
		forceB[i] -= correctionForceB;
		forceC[i] -= correctionForceC;
	}

	// Sum up new total forces for next iteration's correction application
	totalForceA = 0;
	totalForceB = 0;
	totalForceC = 0;
	for (int i = 0; i < XDIM; i++) {
		totalForceA += forceA[i];
		totalForceB += forceB[i];
		totalForceC += forceC[i];
	}
} // end function correctExcessMomentum()


/**
 * This function calculates pressure for the full mixture.
 * Components are coupled together and cannot be separated into meaningful partial pressures.
 */
void calculatePressureCoupled() {
	int i = 0;
	int filterNeighbors = 3;

	// Total pressure
	for (i = 0; i < XDIM; i++) {

		// Bulk pressure
		pressure[i] = nA[i]*theta/(1.-bA*nA[i]-bB*nB[i]-bC*nC[i]) + nB[i]*theta/(1.-bA*nA[i]-bB*nB[i]-bC*nC[i]) + nC[i]*theta/(1.-bA*nA[i]-bB*nB[i]-bC*nC[i])
								- aA*nA[i]*nA[i] - aB*nB[i]*nB[i] - aC*nC[i]*nC[i] - 2*aAB*nA[i]*nB[i] - 2*aAC*nA[i]*nC[i] - 2*aBC*nB[i]*nC[i];

		// Gradient corrections for each single component, including self-interactions
		pressure[i] += -kappa*( nA[i]*laplace(nA,i) + nB[i]*laplace(nB,i) + nC[i]*laplace(nC,i) );
		pressure[i] += -kappa*( 0.5*gradient(nA,i)*gradient(nA,i) + 0.5*gradient(nB,i)*gradient(nB,i) + 0.5*gradient(nC,i)*gradient(nC,i) );
		pressure[i] += kappa*( gradient(nA,i)*gradient(nA,i) + gradient(nB,i)*gradient(nB,i) + gradient(nC,i)*gradient(nC,i) );

		// Gradient corrections for the cross terms (i.e. cross interactions)
		pressure[i] += -kappa*( nA[i]*laplace(nB,i) + nB[i]*laplace(nA,i) + nA[i]*laplace(nC,i) + nC[i]*laplace(nA,i) + nB[i]*laplace(nC,i) + nC[i]*laplace(nB,i) );
		pressure[i] += -kappa*( gradient(nA,i)*gradient(nB,i) + gradient(nA,i)*gradient(nC,i) + gradient(nB,i)*gradient(nC,i) );
		pressure[i] += kappa*( 2.*gradient(nA,i)*gradient(nB,i) + 2.*gradient(nA,i)*gradient(nC,i) + 2.*gradient(nB,i)*gradient(nC,i) );

		pressure[i] *= gammaP;
	}

	// Filtered pressure... some digital noise otherwise
	for (i = 0; i < XDIM; i++) {
		pressureFiltered[i] = filterArrayData(pressure, i, pressureFilterNeighbors);
	}

} // end function calculatePressureCoupled()


void calculateChemicalPotentials() {
	int i = 0;
	double tmp = 0;

	for (i = 0; i < XDIM; i++) {
		tmp = 1.-bA*nA[i]-bB*nB[i]-bC*nC[i];
		volumeExclusion[i] = tmp;

		muA[i] = gammaMu * ( theta*log(nA[i]/(1.-bA*nA[i]-bB*nB[i]-bC*nC[i])) + theta*bA*n[i]/(1.-bA*nA[i]-bB*nB[i]-bC*nC[i])
				- 2*aA*nA[i] - 2*aAB*nB[i] - 2*aAC*nC[i] - kappa*laplace(nA,i)
		- (useMuVdwInteraction ? vdwInteractionFactor : 1)*kappa*laplace(nB,i)
		- (useMuVdwInteraction ? vdwInteractionFactor : 1)*kappa*laplace(nC,i) );
		muB[i] = gammaMu * ( theta*log(nB[i]/(1.-bA*nA[i]-bB*nB[i]-bC*nC[i])) + theta*bB*n[i]/(1.-bA*nA[i]-bB*nB[i]-bC*nC[i])
				- 2*aB*nB[i] - 2*aAB*nA[i] - 2*aBC*nC[i] - kappa*laplace(nB,i)
		- (useMuVdwInteraction ? vdwInteractionFactor : 1)*kappa*laplace(nA,i)
		- (useMuVdwInteraction ? vdwInteractionFactor : 1)*kappa*laplace(nC,i) );
		muC[i] = gammaMu * ( theta*log(nC[i]/(1.-bA*nA[i]-bB*nB[i]-bC*nC[i])) + theta*bC*n[i]/(1.-bA*nA[i]-bB*nB[i]-bC*nC[i])
				- 2*aC*nC[i] - 2*aAC*nA[i] - 2*aBC*nB[i] - kappa*laplace(nC,i)
		- (useMuVdwInteraction ? vdwInteractionFactor : 1)*kappa*laplace(nA,i)
		- (useMuVdwInteraction ? vdwInteractionFactor : 1)*kappa*laplace(nB,i) );

		muNonIdealA[i] = muA[i] - theta*log(nA[i]);
		muNonIdealB[i] = muB[i] - theta*log(nB[i]);
		muNonIdealC[i] = muC[i] - theta*log(nC[i]);
	}

	// Filtered chemical potentials... lots of digital noise otherwise
	for (i = 0; i < XDIM; i++) {
		muAFiltered[i] = filterArrayData(muA, i, muFilterNeighbors);
		muBFiltered[i] = filterArrayData(muB, i, muFilterNeighbors);
		muCFiltered[i] = filterArrayData(muC, i, muFilterNeighbors);
	}

} // end function calculateChemicalPotentials()


void calculateFriction() {
	int i = 0;

	for (i = 0; i < XDIM; i++) {
		frictionA[i] = nA[i]/n[i] * (nB[i]*(uHatA[i]-uHatB[i]) + nC[i]*(uHatA[i]-uHatC[i]));
		frictionB[i] = nB[i]/n[i] * (nA[i]*(uHatB[i]-uHatA[i]) + nC[i]*(uHatB[i]-uHatC[i]));
		frictionC[i] = nC[i]/n[i] * (nA[i]*(uHatC[i]-uHatA[i]) + nB[i]*(uHatC[i]-uHatB[i]));

		forceA[i] += lambda * frictionA[i];
		forceB[i] += lambda * frictionB[i];
		forceC[i] += lambda * frictionC[i];
	} // end for
} // end function calculateFriction()


void collisionForcingNewChemicalPotentialGradient() {
	int i = 0;
	static int lastMuForcingMethod = 0; // default to the first case below

	iterations++;

	calculateMassAndVelocities();
	//	calculateFreeEnergy();
	calculateChemicalPotentials();
	calculatePressureCoupled();

	// Forcing derived from chemical potential gradients... a la Gibbs-Duhem (sum of both equals pressure gradient)
	if (useChemicalPotentialForcingMethod != lastMuForcingMethod) {
		lastMuForcingMethod = useChemicalPotentialForcingMethod;
		switch (useChemicalPotentialForcingMethod) {
		case 1: // grad non-ideal mu
			printf("gradMu forcing: -nx*grad(MuNidx)\n");
			break;
		case 2: // gradient of mu minus ideal pressure
			printf("gradMu forcing: -nx*grad(Mux)-theta*grad(nx)\n");
			break;
		default: // use case 2
			useChemicalPotentialForcingMethod = 2;
			lastMuForcingMethod = 2;
			printf("Invalid Selection!  Using gradMu forcing: -nx*grad(Mux)-theta*grad(nx)\n");
		} // end switch
	} // end if
	switch (useChemicalPotentialForcingMethod) {
	case 1: // grad non-ideal mu
		for (i = 0; i < XDIM; i++) {
			forceA[i] = -1.*nA[i]*gradient(muNonIdealA,i) + nA[i]*g;
			forceB[i] = -1.*nB[i]*gradient(muNonIdealB,i) + nB[i]*g;
			forceC[i] = -1.*nC[i]*gradient(muNonIdealC,i) + nC[i]*g;
		}
		break;
	case 2: // gradient of mu minus ideal pressure/chemical potential gradient (theta rho grad-log rho = theta grad-rho)
		for (i = 0; i < XDIM; i++) {
			forceA[i] = -1. * ( nA[i]*gradient(muA,i)-theta*gradient(nA,i) ) + nA[i]*g;
			forceB[i] = -1. * ( nB[i]*gradient(muB,i)-theta*gradient(nB,i) ) + nB[i]*g;
			forceC[i] = -1. * ( nC[i]*gradient(muC,i)-theta*gradient(nC,i) ) + nC[i]*g;
		}
		break;
	//default: // do nothing...
	} // end switch

	calculateFriction();
	if (momentumCorrectionOn) {
		correctExcessMomentum();
	}

	for (i = 0; i < XDIM; i++) {

		// Correction to the equilibrium distribution that alters the actual forcing
		if (nA[i] !=0) {
			psiA[i] = -oneOverTau * ( (tau-0.25)*forceA[i]*forceA[i]/nA[i] + (1./12.)*laplace(nA,i) );	// subtract psi, so minus sign relative to paper
		}
		else {
			psiA[i] = 0;
		}

		// Calculate particle densities at current lattice spot with forcing included
		fA_0[i] += oneOverTau * ( (nA[i] - nA[i]*theta - nA[i]*uA[i]*uA[i]) - fA_0[i] ) - ( 2.*forceA[i]*uA[i] - psiA[i] );
		fA_1[i] += oneOverTau * ( 0.5*(nA[i]*uA[i]*uA[i]+nA[i]*uA[i]+nA[i]*theta)-fA_1[i]) - ( -forceA[i]*uA[i] - 0.5*forceA[i] + 0.5*psiA[i] );
		fA_2[i] += oneOverTau * ( 0.5*(nA[i]*uA[i]*uA[i]-nA[i]*uA[i]+nA[i]*theta)-fA_2[i]) - ( -forceA[i]*uA[i] + 0.5*forceA[i] + 0.5*psiA[i] );

		// Correction to the equilibrium distribution that alters the actual PGF to pressure is constant in equilibirum
		if (nB[i] != 0) {
			psiB[i] = -oneOverTau * ( (tau-0.25)*forceB[i]*forceB[i]/nB[i] + (1./12.)*laplace(nB,i) );	// subtract psi, so minus sign relative to paper
		}
		else {
			psiB[i] = 0;
		}

		fB_0[i] += oneOverTau * ( (nB[i] - nB[i]*theta - nB[i]*uB[i]*uB[i]) - fB_0[i] ) - ( 2.*forceB[i]*uB[i] - psiB[i] );
		fB_1[i] += oneOverTau * ( 0.5*(nB[i]*uB[i]*uB[i] + nB[i]*uB[i] + nB[i]*theta) - fB_1[i] ) - ( -forceB[i]*uB[i] - 0.5*forceB[i] + 0.5*psiB[i] );
		fB_2[i] += oneOverTau * ( 0.5*(nB[i]*uB[i]*uB[i] - nB[i]*uB[i] + nB[i]*theta) - fB_2[i] ) - ( -forceB[i]*uB[i] + 0.5*forceB[i] + 0.5*psiB[i] );

		// Correction to the equilibrium distribution that alters the actual PGF to pressure is constant in equilibirum
		if (nC[i] != 0) {
			psiC[i] = -oneOverTau * ( (tau-0.25)*forceC[i]*forceC[i]/nC[i] + (1./12.)*laplace(nC,i) );	// subtract psi, so minus sign relative to paper
		}
		else {
			psiC[i] = 0;
		}

		fC_0[i] += oneOverTau * ( (nC[i] - nC[i]*theta - nC[i]*uC[i]*uC[i]) - fC_0[i] ) - ( 2.*forceC[i]*uC[i] - psiC[i] );
		fC_1[i] += oneOverTau * ( 0.5*(nC[i]*uC[i]*uC[i] + nC[i]*uC[i] + nC[i]*theta) - fC_1[i] ) - ( -forceC[i]*uC[i] - 0.5*forceC[i] + 0.5*psiC[i] );
		fC_2[i] += oneOverTau * ( 0.5*(nC[i]*uC[i]*uC[i] - nC[i]*uC[i] + nC[i]*theta) - fC_2[i] ) - ( -forceC[i]*uC[i] + 0.5*forceC[i] + 0.5*psiC[i] );
	} // end for
} // end function collisionForcingNewChemicalPotentialGradient()


void streaming() {
	double tmp;

	/* Original wrap-around end points */
	tmp=fA_1[XDIM-1];                                   // save right end point
	memmove(&fA_1[1],&fA_1[0],(XDIM-1)*sizeof(double)); // shift all cells +1
	fA_1[0]=tmp;                                        // rotate former end to first lattice cell
	tmp=fA_2[0];                                        // save left end point
	memmove(&fA_2[0],&fA_2[1],(XDIM-1)*sizeof(double)); // shift all cells -1
	fA_2[XDIM-1]=tmp;                                   // rotate former first lattice cell to end

	tmp=fB_1[XDIM-1];                                   // save right end point
	memmove(&fB_1[1],&fB_1[0],(XDIM-1)*sizeof(double)); // shift all cells +1
	fB_1[0]=tmp;                                        // rotate former end to first lattice cell
	tmp=fB_2[0];                                        // save left end point
	memmove(&fB_2[0],&fB_2[1],(XDIM-1)*sizeof(double)); // shift all cells -1
	fB_2[XDIM-1]=tmp;                                   // rotate former first lattice cell to end

	tmp=fC_1[XDIM-1];                                   // save right end point
	memmove(&fC_1[1],&fC_1[0],(XDIM-1)*sizeof(double)); // shift all cells +1
	fC_1[0]=tmp;                                        // rotate former end to first lattice cell
	tmp=fC_2[0];                                        // save left end point
	memmove(&fC_2[0],&fC_2[1],(XDIM-1)*sizeof(double)); // shift all cells -1
	fC_2[XDIM-1]=tmp;                                   // rotate former first lattice cell to end

	// Walls at the end points
	// Bounce from lattice origin (0)
	if (!useBoundaryConditionsPeriodic) {
		tmp = fA_1[0];
		fA_1[0] = fA_2[XDIM-1];
		fA_2[XDIM-1]=tmp;

		tmp = fB_1[0];
		fB_1[0] = fB_2[XDIM-1];
		fB_2[XDIM-1]=tmp;

		tmp = fC_1[0];
		fC_1[0] = fC_2[XDIM-1];
		fC_2[XDIM-1]=tmp;
	}
} // end function streaming()


void iteration(){

	// Need to reset the critical and VDW constants each iteration
	// Keeps them all in sync if one is changed during a simulation
	pcA = 3.*tcA/8.; // determine pcA
	ncA = pcA / ((3./8.)*tcA); // fix ncA to be 1
	pcB = (3./8.)*tcB*ncB; // determine pcB
	pcC = (3./8.)*tcC*ncC; // determine pcC

	aA = (27./64.)*(tcA*tcA/pcA);
	aB = (27./64.)*(tcB*tcB/pcB);
	aC = (27./64.)*(tcC*tcC/pcC);
	bA = tcA/(8.*pcA);
	bB = tcB/(8.*pcB);
	bC = tcC/(8.*pcC);
	aAB = sqrt(aA*aB) * vdwInteractionFactor;
	aAC = sqrt(aA*aC) * vdwInteractionFactor;
	aBC = sqrt(aB*aC) * vdwInteractionFactor;

	collision();

	/*Evaporation part */
	//f0[epos]-=(dn-dn*T);
	//f1[epos]-=0.5*(dn*T);
	//f2[epos]-=0.5*(dn*T);

	streaming();

} // end function iteration()

