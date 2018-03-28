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

//	if (i < wall){
//		n1[i] = n1_liquid;
//	}
//	else {
//		n1[i] = n1_gas;
//	}

	int interface = 0.5 * XDIM;
	double transition = 0.0;
	double rhoA1 = theoreticalRhoVapor, rhoA2 = theoreticalRhoLiquid;

	if (i < interface-0.25*XDIM) {
		transition = 0.5 + 0.5*tanh((double)(i%XDIM)/interfaceWidth);
		n1[i] = (1.0-transition)*rhoA2 + transition*rhoA1;
	}
	else if (i >= interface-0.25*XDIM && i <= interface+0.25*XDIM) {
		transition = 0.5 + 0.5*tanh((double)(i-interface)/interfaceWidth);
		n1[i] = (1.0-transition)*rhoA1 + transition*rhoA2;
	}
	else {
		transition = 0.5 + 0.5*tanh((double)((i-XDIM)%XDIM)/interfaceWidth);
		n1[i] = (1.0-transition)*rhoA2 + transition*rhoA1;
	}
}

int getTheoreticalDensities() {
	int readEOF;
	int match = 0;
	double tmpRhoVapor = 0.0, tmpRhoLiquid = 0.0, tmpThetaRatio = 1.0;
	double thetaRatio = theta / tc;
	double matchThreshold = 1e-4;

	FILE *phaseDiagram_densities_twoPhases;
	phaseDiagram_densities_twoPhases = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/single component/phase-data-4/phaseDiagram_rhoVsTemp_theoretical.dat","r");

	if (phaseDiagram_densities_twoPhases) {
		readEOF = fscanf(phaseDiagram_densities_twoPhases, "%17lf %17lf", &tmpRhoVapor, &tmpThetaRatio);
		while (readEOF != EOF) { // find theoretical vapor density
			if (isEqual(tmpThetaRatio, thetaRatio, matchThreshold) && tmpRhoVapor < 1.0) {
				match = 1;
				break;
			}
			readEOF = fscanf(phaseDiagram_densities_twoPhases, "%17lf %17lf", &tmpRhoVapor, &tmpThetaRatio);
		}
		readEOF = fscanf(phaseDiagram_densities_twoPhases, "%17lf %17lf", &tmpRhoLiquid, &tmpThetaRatio);
		while (readEOF != EOF) { // continue to find theoretical liquid density
			if (isEqual(tmpThetaRatio, thetaRatio, matchThreshold) && tmpRhoLiquid > 1.0) {
				match = 2;
				break;
			}
			readEOF = fscanf(phaseDiagram_densities_twoPhases, "%17lf %17lf", &tmpRhoLiquid, &tmpThetaRatio);
		}
		fclose(phaseDiagram_densities_twoPhases);
	}

	if (match == 2) {
		theoreticalRhoVapor = tmpRhoVapor;
		theoreticalRhoLiquid = tmpRhoLiquid;
		printf("t-ratio=%f\trhoV=%f\trhoL=%f\n", thetaRatio, theoreticalRhoVapor, theoreticalRhoLiquid);
	}
	else { // zero out theoretical densities if there are none found
		theoreticalRhoVapor = 0.0;
		theoreticalRhoLiquid = 0.0;
		printf("\nNo matching theoretical densities found... \n\n");
	}

	return match;
} // end function getTheoreticalDensities()


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

	if (initializeProfile == initializeSteps) {
		if (getTheoreticalDensities() && autoKappaGammaMu) {
			interfaceWidth = 1.0 / sqrt(4.0*theoreticalRhoVapor*fabs(tc-theta));
			kappa = 1.0 / (8.0 * theta * theoreticalRhoVapor);
			gammaMu = 1.0 / (6.0 * kappa * theoreticalRhoLiquid);
			gammaMu /= 10.0;
			printf("kappa = %f\tgammaMu = %f\twidth = %f\n", kappa, gammaMu, interfaceWidth);
		}
	}

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

