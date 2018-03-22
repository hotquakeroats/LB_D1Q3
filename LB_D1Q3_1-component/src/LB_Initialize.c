/*
 * LB_Analyze.c
 *
 *  Created on: Jul 7, 2016
 *      Author: clark
 */


#include "LB_D1Q3_1-component.h"


void initializeRandom(int i) {

	n1[i] = n0+Amp*(1.0*rand()/RAND_MAX-0.5);

}


void initializeSteps(int i) {

	if (i < wall){
		n1[i] = n1_liquid;
	}
	else {
		n1[i] = n1_gas;
	}

}


void setInitializeRandom() {
	initializeProfile = initializeRandom;

	// Turn the wall off for this profile
	if (wall > 0) {
		wall *= -1;
	}

	initialize();
}


void setInitializeSteps() {
	initializeProfile = initializeSteps;

	// Turn the wall on for this profile
	if (wall < 0) {
		wall *= -1;
	}

	initialize();
}


void initialize() {

	int i;
	double u = 0;

	iterations = 0;
	rho1 = 0;
	excludedVolume1 = 0;

	pc = 3.*tc/8.;	// This is needed to put the original critical parameter formulation on equal footing with the VDW constant formulation
	nc = pc / ((3./8.)*tc);

	// Reset the values of the VDW constants for each component
	a1 = (27./64.)*(tc*tc/pc);
	b1 = tc/(8.*pc);

	for (i = 0; i < XDIM; i++) {

		initializeProfile(i);
		u1[i] = 0;
		uHat1[i] = 0;

		rho1 += n1[i];

		// Initialize 1st component
		f1_0[i] = n1[i] - n1[i]*T0 - n1[i]*u*u;            		// zero velocity
		f1_1[i] = (1./2) * (n1[i]*u*u + n1[i]*u + n1[i]*T0); 	// +1 velocity (moving right)
		f1_2[i] = (1./2) * (n1[i]*u*u - n1[i]*u + n1[i]*T0); 	// -1 velocity (moving left)

		mu1[i] = 0;
		muNonIdeal1[i] = 0;

		A[i] = 0;

		pressure[i] = 0;
		pressure1[i] = 0;
		pressureNonIdeal1[i] = 0;
		pressureCorrected1[i] = 0;
		p[i] = 0;
		pni[i] = 0;
		pf[i] = 0;
		PF[i] = 0;
		ddni[i] = 0;
	}

	// Set the total density and excluded volume constants for the VDW equations
	excludedVolume1 = rho1 * b1;

} // end function initialize()


// This init() function is only called by main() a single time
void init() {

	// Turn the wall off for the initial start-up profile
	if (wall > 0) {
		wall *= -1;
	}

	//collision = collisionForcingNewPressureGradient;
	setCollisionForcingNewChemicalPotentialGradient();
	//initializeProfile = initializeRandom;
	setInitializeRandom();

	initialize();

}

