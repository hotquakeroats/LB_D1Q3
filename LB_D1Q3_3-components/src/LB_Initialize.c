/*
 * LB_Initialize.c
 *
 *  Created on: Jul 7, 2016
 *      Author: clark
 */


#include "LB_D1Q3_3-components.h"


void initializeRandom(int i) {
	// Make sure fluctuations don't give negative particle counts
	if (nA0 < Amp || nB0 < Amp || nC0 < Amp) {
		Amp = 0.1 * (nA0 < nB0 ? (nA0 < nC0 ? nA0 : nC0) : nB0);
	}

	nA[i] = nA0 + Amp*(1.0*rand()/RAND_MAX-0.5);
	nB[i] = nB0 + Amp*(1.0*rand()/RAND_MAX-0.5);
	nC[i] = nC0 + Amp*(1.0*rand()/RAND_MAX-0.5);

} // end function initializeRandom()


void initializeStepProfile(int i) {
	int width = (int)(interfaceWidth * XDIM); // 15
	int half = (int)(0.5 * XDIM * PHASE_WIDTH); // 25
	double transition = 1.0;

//	if (i < interface1-half) {
//		transition = 0.5 + 0.5*tanh((double)((i+XDIM-1-interface6)%XDIM)/(double)width);
//		nA[i] = (1.-transition)*theoreticalDensityA4 + transition*theoreticalDensityA1;
//		nB[i] = (1.-transition)*theoreticalDensityB4 + transition*theoreticalDensityB1;
//		nC[i] = (1.-transition)*theoreticalDensityC4 + transition*theoreticalDensityC1;
//	}
//	else if (i >= interface1-half && i < interface1+half) {
//		transition = 0.5 + 0.5*tanh((double)(i-interface1)/(double)width);
//		nA[i] = (1.-transition)*theoreticalDensityA1 + transition*theoreticalDensityA4;
//		nB[i] = (1.-transition)*theoreticalDensityB1 + transition*theoreticalDensityB4;
//		nC[i] = (1.-transition)*theoreticalDensityC1 + transition*theoreticalDensityC4;
//	}
//	else if (i >= interface2-half && i < interface2+half) {
//		transition = 0.5 + 0.5*tanh((double)(i-interface2)/(double)width);
//		nA[i] = (1.-transition)*theoreticalDensityA4 + transition*theoreticalDensityA2;
//		nB[i] = (1.-transition)*theoreticalDensityB4 + transition*theoreticalDensityB2;
//		nC[i] = (1.-transition)*theoreticalDensityC4 + transition*theoreticalDensityC2;
//	}
//	else if (i >= interface3-half && i < interface3+half) {
//		transition = 0.5 + 0.5*tanh((double)(i-interface3)/(double)width);
//		nA[i] = (1.-transition)*theoreticalDensityA2 + transition*theoreticalDensityA4;
//		nB[i] = (1.-transition)*theoreticalDensityB2 + transition*theoreticalDensityB4;
//		nC[i] = (1.-transition)*theoreticalDensityC2 + transition*theoreticalDensityC4;
//	}
//	else if (i >= interface4-half && i < interface4+half) {
//		transition = 0.5 + 0.5*tanh((double)(i-interface4)/(double)width);
//		nA[i] = (1.-transition)*theoreticalDensityA4 + transition*theoreticalDensityA3;
//		nB[i] = (1.-transition)*theoreticalDensityB4 + transition*theoreticalDensityB3;
//		nC[i] = (1.-transition)*theoreticalDensityC4 + transition*theoreticalDensityC3;
//	}
//	else if (i >= interface5-half && i < interface5+half) {
//		transition = 0.5 + 0.5*tanh((double)(i-interface5)/(double)width);
//		nA[i] = (1.-transition)*theoreticalDensityA3 + transition*theoreticalDensityA4;
//		nB[i] = (1.-transition)*theoreticalDensityB3 + transition*theoreticalDensityB4;
//		nC[i] = (1.-transition)*theoreticalDensityC3 + transition*theoreticalDensityC4;
//	}
//	else { // transition back to vapor
//		transition = 0.5 + 0.5*tanh((double)((i-interface6)%XDIM)/(double)width);
//		nA[i] = (1.-transition)*theoreticalDensityA4 + transition*theoreticalDensityA1;
//		nB[i] = (1.-transition)*theoreticalDensityB4 + transition*theoreticalDensityB1;
//		nC[i] = (1.-transition)*theoreticalDensityC4 + transition*theoreticalDensityC1;
//	}

	if (i < interface1) { // 0-49
		transition = 0.5 + 0.5*tanh((double)(i-interface1+half)/(double)width);
		nA[i] = (1.-transition)*theoreticalDensityA4 + transition*theoreticalDensityA1;
		nB[i] = (1.-transition)*theoreticalDensityB4 + transition*theoreticalDensityB1;
		nC[i] = (1.-transition)*theoreticalDensityC4 + transition*theoreticalDensityC1;
	}
	else if (i >= interface1 && i < interface2) { // 50-99
		transition = 0.5 + 0.5*tanh((double)(i-interface2+half)/(double)width);
		nA[i] = (1.-transition)*theoreticalDensityA1 + transition*theoreticalDensityA4;
		nB[i] = (1.-transition)*theoreticalDensityB1 + transition*theoreticalDensityB4;
		nC[i] = (1.-transition)*theoreticalDensityC1 + transition*theoreticalDensityC4;
	}
	else if (i >= interface2 && i < interface3) { // 100-149
		transition = 0.5 + 0.5*tanh((double)(i-interface3+half)/(double)width);
		nA[i] = (1.-transition)*theoreticalDensityA4 + transition*theoreticalDensityA2;
		nB[i] = (1.-transition)*theoreticalDensityB4 + transition*theoreticalDensityB2;
		nC[i] = (1.-transition)*theoreticalDensityC4 + transition*theoreticalDensityC2;
	}
	else if (i >= interface3 && i < interface4) { // 150-199
		transition = 0.5 + 0.5*tanh((double)(i-interface4+half)/(double)width);
		nA[i] = (1.-transition)*theoreticalDensityA2 + transition*theoreticalDensityA4;
		nB[i] = (1.-transition)*theoreticalDensityB2 + transition*theoreticalDensityB4;
		nC[i] = (1.-transition)*theoreticalDensityC2 + transition*theoreticalDensityC4;
	}
	else if (i >= interface4 && i < interface5) { // 200-249
		transition = 0.5 + 0.5*tanh((double)(i-interface5+half)/(double)width);
		nA[i] = (1.-transition)*theoreticalDensityA4 + transition*theoreticalDensityA3;
		nB[i] = (1.-transition)*theoreticalDensityB4 + transition*theoreticalDensityB3;
		nC[i] = (1.-transition)*theoreticalDensityC4 + transition*theoreticalDensityC3;
	}
	else { // 250-299
		transition = 0.5 + 0.5*tanh((double)((i-interface6+half)%XDIM)/(double)width);
		nA[i] = (1.-transition)*theoreticalDensityA3 + transition*theoreticalDensityA4;
		nB[i] = (1.-transition)*theoreticalDensityB3 + transition*theoreticalDensityB4;
		nC[i] = (1.-transition)*theoreticalDensityC3 + transition*theoreticalDensityC4;
	}

} // end function initializeTheoreticalThreePhases()


void setInitializeRandom() {
	initializeProfile = initializeRandom;
	initialize();
}


void setInitializeStepProfile() {
	initializeProfile = initializeStepProfile;
	initialize();
}


void initialize() {
	int i;

	iterations = 0;
	total_phase_iterations = phase_iterations;
	volumeTotal = 1.0;
	totalForceA = 0;
	totalForceB = 0;
	totalForceC = 0;
	correctionForceA = 0;
	correctionForceB = 0;
	correctionForceC = 0;

	//
	// VDW parameters and critical constants need to be in sync with the LB iteration function
	//

	// tcA, tcB, and ncB vary as the degrees of freedom
	pcA = 3.*tcA/8.; // determine pcA
	ncA = pcA / ((3./8.)*tcA); // fix ncA to be 1
	pcB = (3./8.)*tcB*ncB; // determine pcB
	pcC = (3./8.)*tcC*ncC; // determine pcC

	// Reset the values of the VDW constants for each component
	aA = (27./64.)*(tcA*tcA/pcA);
	aB = (27./64.)*(tcB*tcB/pcB);
	aC = (27./64.)*(tcC*tcC/pcC);
	bA = tcA/(8.*pcA);
	bB = tcB/(8.*pcB);
	bC = tcC/(8.*pcC);
	aAB = sqrt(aA*aB) * vdwInteractionFactor;
	aAC = sqrt(aA*aC) * vdwInteractionFactor;
	aBC = sqrt(aB*aC) * vdwInteractionFactor;

	lnExplosion = 0;

	for (i = 0; i < XDIM; i++) {
		initializeProfile(i);
		n[i] = nA[i] + nB[i] + nC[i];
		uA[i] = 0;
		uB[i] = 0;
		uC[i] = 0;
		uHatA[i] = 0;
		uHatB[i] = 0;
		uHatC[i] = 0;
		u[i] = 0;

		// Initialize 1st component
		fA_0[i] = nA[i] - nA[i]*T0 - nA[i]*uA[i]*uA[i];            		// zero velocity
		fA_1[i] = 0.5 * (nA[i]*uA[i]*uA[i] + nA[i]*uA[i] + nA[i]*T0); 	// +1 velocity (moving right)
		fA_2[i] = 0.5 * (nA[i]*uA[i]*uA[i] - nA[i]*uA[i] + nA[i]*T0); 	// -1 velocity (moving left)
		psiA[i] = 0;
		muA[i] = 0;
		muNonIdealA[i] = 0;

		// Initialize 2nd component
		fB_0[i] = nB[i] - nB[i]*T0 - nB[i]*uB[i]*uB[i];            		// zero velocity
		fB_1[i] = 0.5 * (nB[i]*uB[i]*uB[i] + nB[i]*uB[i] + nB[i]*T0); 	// +1 velocity (moving right)
		fB_2[i] = 0.5 * (nB[i]*uB[i]*uB[i] - nB[i]*uB[i] + nB[i]*T0); 	// -1 velocity (moving left)
		psiB[i] = 0;
		muB[i] = 0;
		muNonIdealB[i] = 0;

		// Initialize 3nd component
		fC_0[i] = nC[i] - nC[i]*T0 - nC[i]*uC[i]*uC[i];            		// zero velocity
		fC_1[i] = 0.5 * (nC[i]*uC[i]*uC[i] + nC[i]*uC[i] + nC[i]*T0); 	// +1 velocity (moving right)
		fC_2[i] = 0.5 * (nC[i]*uC[i]*uC[i] - nC[i]*uC[i] + nC[i]*T0); 	// -1 velocity (moving left)
		psiC[i] = 0;
		muC[i] = 0;
		muNonIdealC[i] = 0;



		// Also scale the theoretical expectations by the appropriate gamma factor
		//		theoreticalMuAArray[i] = theoreticalMuA * gammaMu;
		//		theoreticalMuBArray[i] = theoreticalMuB * gammaMu;
		//		theoreticalPressureArray[i] = theoreticalPressure * gammaP;
		pressure[i] = 0;
		volumeExclusion[i] = 0;
		frictionA[i] = 0;
		frictionB[i] = 0;
		frictionC[i] = 0;
		forceA[i] = 0;
		forceB[i] = 0;
		forceC[i] = 0;
		tg[i] = 0;
	}

} // end function initialize()


// This init() function is only called by main() a single time
void init() {
	collision = collisionForcingNewChemicalPotentialGradient;
	initializeProfile = initializeRandom;
	printf("Initializing with gradMu forcing...\n");
	initialize();
}

