/*
 * LB_D1Q3_2components.c
 *
 *  Created on: Apr 22, 2016
 *      Author: Kent Ridl
 */


#include "LB_D1Q3_1-component.h"


// Initial densities of each component
double n1_liquid = 1.27;
double n1_gas = 0.6;

// Total amount of each component
double rho1 = 0;
double excludedVolume1 = 0;

int freeEnergyArraySize = 0;

// Evaporation constants
double dn=0; 				// delta n for evaporation
int epos=0;  				// lattice position at which evaporation occurs

// Critical point traits
double tc = .34;
double nc = 1.;
double pc = 1.;

// Equations of motion constants
double T0=0.33333333;  		// initial theta from write-up
double theta = 1./3.;
double n0=1.0;
double Amp=0.01;
double oneOverTau = 1;		// 1/tau from write-up (relaxation constant)
double tau = 1;
double g = 0;          		// gravitational acceleration term
double lambda = 1;     		// friction (F12) coefficient
double gammaP = 1;			// pressure coefficient for forcing rate
double gammaMu = 1;			// chemical potential coefficient for forcing rate
double kappa = 0.1;
double dummy = 1;
double quenchDepth = 0.85;
double pressureMethodCorrection = 0;
double pressureMethodCoefficient = 7./3.;	// u*gradRho + u*gradRho(transpose) + theta*u*gradRho (2 1/3)

double a1 = 0.1;			// VDW constants for each component, interaction
double b1 = 1./3.;

int usePressureCoupled = 0;
int usePressureCriticalParameters = 0;			// default use VDW constants
int useChemicalPotentialCriticalParameters = 0;	// default use VDW constants
int useChemicalPotentialsCoupled = 0;
int useChemicalPotentialNonIdeal = 0;			// default use full mu (not non-ideal part)
int lnExplosion = 0;
int useBoundaryConditionsPeriodic = 1;			// default periodic BCs

// GUI
int ulreq=0;
int ugreq=0;
int pgreq=0;
int tgreq=0;
int mu1greq = 0; // graph requests
int next=0;
int Pause=1;
int done=0;
int Repeat=10;
int iterations;
int collectData = 0;
int wall = XDIM * 0.5;


void calculateMassAndVelocities() {

	int i = 0;

	// Sweep across the lattice to conserve mass and determine the pressures at each cell
	for (i = 0; i < XDIM; i++) {

		n1[i] = f1_0[i] + f1_1[i] + f1_2[i]; 		// 1st component mass density for this step

		if (n1[i] != 0) {
			u1[i] = (f1_1[i]-f1_2[i]) / n1[i];

			// Correction to the display of the mean fluid velocity
			// Needed for the forcing methods (pressure method is good with the above)
			uHat1[i] = u1[i] + 0.5/n1[i]*F1[i];
		}
		else {
			u1[i] = 0;
		}
	}
}


/**
 * This function calculates pressures for each component.
 * It leaves the pressures uncoupled; each components reacts independent from the others.
 */
void calculatePressurePartial() {

	int i = 0;

	for (i = 0; i < XDIM; i++) {
		if (usePressureCriticalParameters) {
			pressure1[i] = n1[i]/(3.-n1[i]) - (9./8.)*tc*n1[i]*n1[i];	// van der Waals gas pressure term
		}
		else {
			pressure1[i] = n1[i]*theta + (b1*n1[i]*n1[i]*theta)/(1.-b1*n1[i]) - a1*n1[i]*n1[i];
		}
	}

	// Corrections for the pressure gradient terms
	for (i = 0; i < XDIM; i++) {	// this correction makes it equal to Alexander's p[] from paper, not pf[] or PF[]
		correctionsPressure1[i] = -kappa*( n1[i]*laplace(n1,i) + (1./2.)*gradient(n1,i)*gradient(n1,i) ) + ( kappa*gradient(n1,i)*gradient(n1,i) );
		pressure1[i] += correctionsPressure1[i];
		pressureNonIdeal1[i] = (gammaP * pressure1[i]) - n1[i]*theta;	// non-ideal part of pressure (vdw - ideal); provides a force

		pressure[i] = gammaP * pressure1[i]; // + pressure2[i]);  // Add partial pressures as if these were ideal gases (which they are not)

		if (usePressureCriticalParameters) {
			pressureCriticalConstants[i] = pressure[i];
		}
		else {
			pressureVDWConstants[i] = pressure[i];
		}
	}
}


void calculatePressureTest() {
	for (int i = 0; i < XDIM; i++) {
		pressureTest[i] = n1[i]*theta + (b1*n1[i]*n1[i]*theta)/(1.-b1*n1[i]) - a1*n1[i]*n1[i];
		//
		// Does something
		//
		pressureTest1[i] = pressureTest[i] - 0.5*kappa * n1[i]*laplace(n1,i); // slightly displaced from next down; better if 1.0*kappa instead
		//		pressureTest2[i] = pressureTest[i] - 0.5*kappa * laplace(n1,i); // slightly displaced from next up?
		pressureTest3[i] = pressureTest[i] - 0.5*kappa * n1[i]*(laplace(n1,i) + laplace(n1,(i+1)%XDIM)); // (a1) mirror of next down; virtually identical to (b1)
		//		pressureTest4[i] = pressureTest[i] - 0.5*kappa * n1[i]*(laplace(n1,i) + laplace(n1,(i-1+XDIM)%XDIM)); // (a2) mirror of next up; virtually identical to (b2)
		//		pressureTest5[i] = pressureTest[i] - 0.5*kappa * n1[(i+1)%XDIM]*(laplace(n1,i) + laplace(n1,(i+1)%XDIM)); // (b1) mirror; same
		//		pressureTest6[i] = pressureTest[i] - 0.5*kappa * n1[(i-1+XDIM)%XDIM]*(laplace(n1,i) + laplace(n1,(i-1+XDIM)%XDIM)); // (b2) mirror; same

		// test to see if g1-g2 = 0... "almost"... robs a bit from 1 peak to add to the other
		pressureTestZero[i]  = pressureTest[i] - 0.5*kappa * (n1[i]-n1[(i+1)%XDIM]) * (laplace(n1,i)+laplace(n1,(i+1)%XDIM));

		//
		// Does nothing
		//
		//
		//		pressureTest[i] += 0.5*kappa * n1[i]*(laplace(n1,i) + laplace(n1,(i+1)%XDIM));
		//		pressureTest[i] += 0.5*kappa * n1[i]*(laplace(n1,i) + laplace(n1,(i-1+XDIM)%XDIM));
		//		pressureTest[i] += 0.5*kappa * n1[(i+1)%XDIM]*(laplace(n1,i) + laplace(n1,(i+1)%XDIM));
		//		pressureTest[i] += 0.5*kappa * n1[(i-1+XDIM)%XDIM]*(laplace(n1,i) + laplace(n1,(i-1+XDIM)%XDIM));
		//		pressureTest[i] += -0.5*kappa * laplace(n1,i)*(n1[i] - n1[(i+1)%XDIM]);
		//		pressureTest[i] += 0.5*kappa * laplace(n1,i)*(n1[i] - n1[(i+1)%XDIM]);
		//		pressureTest[i] += -0.5*kappa * laplace(n1,i)*(n1[i] - n1[(i-1+XDIM)%XDIM]);
		//		pressureTest[i] += 0.5*kappa * laplace(n1,i)*(n1[i] - n1[(i-1+XDIM)%XDIM]);

		//
		// Crashes
		//
		//		pressureTest[i] += -0.5*kappa * laplace((i+1)%XDIM,i)*(n1[i] - n1[(i+1)%XDIM]);
		//		pressureTest[i] += 0.5*kappa * laplace((i+1)%XDIM,i)*(n1[i] - n1[(i+1)%XDIM]);
		//		pressureTest[i] += -0.5*kappa * laplace((i-1+XDIM)%XDIM,i)*(n1[i] - n1[(i-1+XDIM)%XDIM]);
		//		pressureTest[i] += 0.5*kappa * laplace((i-1+XDIM)%XDIM,i)*(n1[i] - n1[(i-1+XDIM)%XDIM]);

	}
}


void calculatePressures() {

	int i;
	int ip = 0;
	int im = 0;
	//	double dn;
	//	double ddn;
	//	double dp;
	//	double ddp;
	double omega = 1./oneOverTau;

	// Functions to calculate pressures, partial and full
	calculatePressurePartial();
	calculatePressureTest();

	for (i = 0; i < XDIM; i++) {
		pressureCCMinusVDW[i] = pressureCriticalConstants[i] - pressureVDWConstants[i];
	}

	// For comparison, these loops calculate pressures in the same manner as Alexander's code
	// Modified to include the gamma "modulation"
	for (i = 0; i < XDIM; i++) {
		ip=(i+1)%XDIM;
		im=(i+XDIM-1)%XDIM;

		dni[i]=gradient(n1,i);
		ddni[i]=laplace(n1,i);

		pni[i]=Pni(n1[i],dni[i],ddni[i]);
		p[i] = P(n1[i],dni[i],ddni[i]);

		p[i] *= gammaP;
		pni[i] *= gammaP;
	}

	for (i = 0; i < XDIM; i++) {
		ip=(i+1)%XDIM;
		im=(i+XDIM-1)%XDIM;

		dpi[i]=gradient(pni,i);
		ddpi[i]=laplace(p,i);

		Forces[i]=-dpi[i];

		pf[i]=p[i]-0.25*Forces[i]*Forces[i]/n1[i]-(1./omega-0.5)*Forces[i]*Forces[i]/n1[i]+0.25*ddpi[i]-1./12.*ddni[i];
		PF[i]=p[i]+0.25*Forces[i]*Forces[i]/n1[i]+0.25*ddpi[i]-1./12.*ddni[i];
	}

}


void calculateChemicalPotentials() {

	int i = 0;
	double tmp = 0;

	if (useChemicalPotentialCriticalParameters) {
		for (i = 0; i < XDIM; i++) {
			if (n1[i] != 0) {
				mu1[i] = gammaMu * (-theta*log(1./n1[i]-1./(3*nc)) + theta/(1-(n1[i]/(3*nc))) - (9./4.)*tc*n1[i]/nc - kappa*laplace(n1,i));
			}
			else {
				mu1[i] = 0; // TODO: Is zeroing this out having an adverse effect?
			}

			muCriticalConstants1[i] = mu1[i];
		}
	}
	else if (!useChemicalPotentialCriticalParameters) {
		// These chemical potentials need density gradients, so calculate them separately
		for (i = 0; i < XDIM; i++) {
			tmp = 1.-b1*n1[i]; //-b2*n2[i];
			if (tmp < 0) {
				lnExplosion = 1;
				printf("LN EXPLOSION!!\t%f\n", tmp);
			}

			mu1[i] = gammaMu * ( theta*log(n1[i]/(1.-b1*n1[i])) + theta*b1*n1[i]/(1.-b1*n1[i]) + theta - 2.*a1*n1[i] - kappa*laplace(n1,i) );
			muVDWConstants1[i] = mu1[i];
		}
	}

	// These chemical potentials need density gradients, so calculate them separately
	for (i = 0; i < XDIM; i++) {
		muNonIdeal1[i] = mu1[i] - theta*log(n1[i]);
		muCCMinusVDW1[i] = muCriticalConstants1[i] - muVDWConstants1[i];
	}

}


void collisionForcingNewPressureGradient() {

	int i = 0;

	iterations++;

	calculateMassAndVelocities();
	calculatePressures();
	calculateChemicalPotentials();

	// Forcing derived from pressure gradients... a la Gibbs-Duhem (pressure gradient - partial pressure - equals chemical potential gradients)
	for (i = 0; i < XDIM; i++) {
		// PGF plus gravity term for each component
		F1[i] = -1.*gradient(pressureNonIdeal1,i) + n1[i]*g;
	}

	// Correction for the effective pressure that is exerted by including the gradient terms
	// This correction does not change the PGF; it only effects what is displayed in a variable separate from the actual pressure
	for (i = 0; i < XDIM; i++) {
		// TODO: the corrected pressure has a bug in the boundary conditions
		if (n1[i] != 0) {
			correctionsPressureDisplay1[i] = -( tau-(1./4.) )*F1[i]*F1[i]/n1[i] + (1./4.)*( laplace(pressure1,i)-theta*laplace(n1,i) );
		}
		else {
			correctionsPressureDisplay1[i] = 0;
		}
		pressureCorrected1[i] = pressure1[i] + correctionsPressureDisplay1[i];	// used to display the effective pressure, not to calculate forces
		pressureGradPMethod[i] = pressure[i];

		muGradPMethod[i] = mu1[i];
	}

	for (i = 0; i < XDIM; i++) {
		F1GradPMethod[i] = F1[i];
		F1GradPGradMuDifference[i] = F1GradPMethod[i] - F1GradMuMethod[i];
	}

	compareGradPRhoGradMu();

	for (i = 0; i < XDIM; i++) {
		// Correction to the equilibrium distribution that alters the actual PGF to pressure is constant in equilibirum
		if (n1[i] !=0) {
			psi1[i] = -oneOverTau * ( (tau-(1./4.))*F1[i]*F1[i]/n1[i] + (1./12.)*laplace(n1,i) );	// subtract psi, so minus sign relative to paper
		}
		else {
			psi1[i] = 0;
		}

		// Calculate particle densities at current lattice spot with forcing included
		f1_0[i] += oneOverTau * ( (n1[i] - n1[i]*theta - n1[i]*u1[i]*u1[i]) - f1_0[i] ) - ( 2.*F1[i]*u1[i] - psi1[i] );
		f1_1[i] += oneOverTau * ( (1./2.)*(n1[i]*u1[i]*u1[i]+n1[i]*u1[i]+n1[i]*theta)-f1_1[i] ) - ( -F1[i]*u1[i] - (1./2.)*F1[i] + (1./2.)*psi1[i] );
		f1_2[i] += oneOverTau * ( (1./2.)*(n1[i]*u1[i]*u1[i]-n1[i]*u1[i]+n1[i]*theta)-f1_2[i] ) - ( -F1[i]*u1[i] + (1./2.)*F1[i] + (1./2.)*psi1[i] );
	}

} // end function collision()


/**
 * This pressure method based on the Holdych-new method cited in Alexander's thermodynamic consistency paper (with the Corr laplace(n) term)
 */
void collisionPressureMethod() {

	int i = 0;

	iterations++;

	calculateMassAndVelocities();
	calculatePressures();
	calculateChemicalPotentials();

	for (i = 0; i < XDIM; i++) {
		pressurePressureMethod[i] = pressure[i];
		muPressureMethod[i] = mu1[i];

		// Calculate the pressure method correction A
		A[i] = pressureNonIdeal1[i] + (tau-0.5)*pressureMethodCoefficient*(n1[i]*u1[i]*gradient(n1,i) / (n1[i] + pressureMethodCorrection*laplace(n1,i)));

		// Calculate particle densities at current lattice spot with pressure method correction included
		f1_0[i] += oneOverTau * ( (n1[i] - n1[i]*theta - n1[i]*u1[i]*u1[i] - A[i]) - f1_0[i] );
		f1_1[i] += oneOverTau * ( (1./2.)*(n1[i]*u1[i]*u1[i] + n1[i]*u1[i] + n1[i]*theta + A[i]) - f1_1[i] );
		f1_2[i] += oneOverTau * ( (1./2.)*(n1[i]*u1[i]*u1[i] - n1[i]*u1[i] + n1[i]*theta + A[i]) - f1_2[i] );
	}

} // end function collision()


void collisionForcingNewChemicalPotentialGradient() {

	int i = 0;

	iterations++;

	calculateMassAndVelocities();
	calculatePressures();
	calculateChemicalPotentials();

	for (i = 0; i < XDIM; i++) {
		pressureGradMuMethod[i] = pressure[i];
		muGradMuMethod[i] = mu1[i];
	}

	// Forcing derived from chemical potential gradients... a la Gibbs-Duhem (sum of both equals pressure gradient)
	if (useChemicalPotentialNonIdeal) { // gradient of non-ideal mu
		for (i = 0; i < XDIM; i++) {
			F1[i] = -1.*n1[i]*gradient(muNonIdeal1,i) + n1[i]*g;
		}
	}
	else {
		for (i = 0; i < XDIM; i++) { // gradient of FULL mu minus ideal pressure!
			F1[i] = -1. * ( n1[i]*gradient(mu1,i)-theta*gradient(n1,i) ) + n1[i]*g;
		}
	}

	for (i = 0; i < XDIM; i++) {
		F1GradMuMethod[i] = F1[i];
		F1GradPGradMuDifference[i] = F1GradPMethod[i] - F1GradMuMethod[i];
	}

	compareGradPRhoGradMu();

	for (i = 0; i < XDIM; i++) {

		// Correction to the equilibrium distribution that alters the actual forcing
		if (n1[i] !=0) {
			psi1[i] = -oneOverTau * ( (tau-(1./4.))*F1[i]*F1[i]/n1[i] + (1./12.)*laplace(n1,i) );	// subtract psi, so minus sign relative to paper
		}
		else {
			psi1[i] = 0;
		}

		// Calculate particle densities at current lattice spot with forcing included
		f1_0[i] += oneOverTau * ( (n1[i] - n1[i]*theta - n1[i]*u1[i]*u1[i]) - f1_0[i] ) - ( 2.*F1[i]*u1[i] - psi1[i] );
		f1_1[i] += oneOverTau * ( (1./2.)*(n1[i]*u1[i]*u1[i]+n1[i]*u1[i]+n1[i]*theta)-f1_1[i]) - ( -F1[i]*u1[i] - (1./2.)*F1[i] + (1./2.)*psi1[i] );
		f1_2[i] += oneOverTau * ( (1./2.)*(n1[i]*u1[i]*u1[i]-n1[i]*u1[i]+n1[i]*theta)-f1_2[i]) - ( -F1[i]*u1[i] + (1./2.)*F1[i] + (1./2.)*psi1[i] );
	}

} // end function collision()

// TODO: dynamic kappa, gammaP/gammaMu values
// TODO: sync up changes from this single component code with the multi-component code
// TODO: look at relative stability among the 3 methods (pressure, gradP, gradMu)
// TODO: even/odd lattice problem comments for thesis
void streaming() {

	double tmp;

	/* Original wrap-around end points */

	tmp=f1_1[XDIM-1];                                   // save right end point
	memmove(&f1_1[1],&f1_1[0],(XDIM-1)*sizeof(double)); // shift all cells +1
	f1_1[0]=tmp;                                        // rotate former end to first lattice cell
	tmp=f1_2[0];                                        // save left end point
	memmove(&f1_2[0],&f1_2[1],(XDIM-1)*sizeof(double)); // shift all cells -1
	f1_2[XDIM-1]=tmp;                                   // rotate former first lattice cell to end

	/* Walls at the end points */
	// Bounce from lattice origin (0)
	if (!useBoundaryConditionsPeriodic) {
		tmp = f1_1[0];
		f1_1[0] = f1_2[XDIM-1];
		f1_2[XDIM-1]=tmp;
	}

	/* Switching cells to simulate wall between lattice spots [wall-1]|[wall] */
	// TODO: this works for wall == 1 (correct) through wall == XDIM (bug)... should be restricted to wall == 1 through XDIM-1
	if (wall > 0) {
		tmp = f1_1[wall];
		f1_1[wall] = f1_2[(wall-1+XDIM)%XDIM];
		f1_2[(wall-1+XDIM)%XDIM] = tmp;
	}

} // end function streaming()


void iteration(){

	// Need to reset the critical and VDW constants each iteration
	// Keeps them all in sync if one is changed during a simulation
	pc = 3.*tc/8.;
	nc = pc / ((3./8.)*tc);
	a1 = (27./64.)*(tc*tc/pc);
	b1 = tc/(8.*pc);

	collision();

	/*Evaporation part */
	//f0[epos]-=(dn-dn*T);
	//f1[epos]-=0.5*(dn*T);
	//f2[epos]-=0.5*(dn*T);

	streaming();

} // end function iteration()

