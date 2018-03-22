/*
 * LB_D1Q3_2-components.h
 *
 *  Created on: Jul 7, 2016
 *      Author: clark
 */

#ifndef LB_D1Q3_2_COMPONENTS_H_
#define LB_D1Q3_2_COMPONENTS_H_

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>
#include <sys/stat.h>

#include "mygraph.h"

#define XDIM 300 // lattice length
#define PHASE_WIDTH (1./6.)

#define abs(a) ( (a < 0) ? -(a) : (a) )
#define isEqual(a, b, threshold) ( (abs((a)-(b)) < (threshold)) ? 1 : 0 )

// Equilibrium distribution of each component velocities
double fA_0[XDIM]; // the static vector
double fA_1[XDIM]; // +1 lattice space
double fA_2[XDIM]; // -1 lattice space
double fB_0[XDIM]; // the static vector
double fB_1[XDIM]; // +1 lattice space
double fB_2[XDIM]; // -1 lattice space
double fC_0[XDIM]; // the static vector
double fC_1[XDIM]; // +1 lattice space
double fC_2[XDIM]; // -1 lattice space

double psiA[XDIM];
double psiB[XDIM];
double psiC[XDIM];

// Number densities of each component
double nA[XDIM]; // 1st component mass/density
double nB[XDIM]; // 2nd component mass/density
double nC[XDIM];
double n[XDIM];  // total mass/density (rho)

//double freeEnergy[XDIM];

// Component and bulk velocities
double uA[XDIM];
double uB[XDIM];
double uC[XDIM];
double uHatA[XDIM];
double uHatB[XDIM];
double uHatC[XDIM];
double u[XDIM];

// Component thermodynamic properties
double muA[XDIM]; // chemical potential
double muNonIdealA[XDIM];
double muAFiltered[XDIM];
double muB[XDIM];
double muNonIdealB[XDIM];
double muBFiltered[XDIM];
double muC[XDIM];
double muNonIdealC[XDIM];
double muCFiltered[XDIM];
double theoreticalMuAArray[XDIM];
double theoreticalMuBArray[XDIM];
double theoreticalMuCArray[XDIM];
double pressure[XDIM];
double pressureFiltered[XDIM];
double theoreticalPressureArray[XDIM];
double volumeExclusion[XDIM];

// Forces on each component; potential and friction
double forceA[XDIM];
double forceB[XDIM];
double forceC[XDIM];
double frictionA[XDIM];
double frictionB[XDIM];
double frictionC[XDIM];

// Graph/GUI variables
double tg[XDIM];


//
// Externs
//

// Initial densities of each component
extern double theoreticalDensityA1;
extern double theoreticalDensityA2;
extern double theoreticalDensityA3;
extern double theoreticalDensityA4;
extern double theoreticalDensityB1;
extern double theoreticalDensityB2;
extern double theoreticalDensityB3;
extern double theoreticalDensityB4;
extern double theoreticalDensityC1;
extern double theoreticalDensityC2;
extern double theoreticalDensityC3;
extern double theoreticalDensityC4;
extern int interface1;
extern int interface2;
extern int interface3;
extern int interface4;
extern int interface5;
extern int interface6;

// Total amount of each component
extern double volumeTotal;
extern double totalForceA;
extern double totalForceB;
extern double totalForceC;
extern double correctionForceA;
extern double correctionForceB;
extern double correctionForceC;

extern double pressureFilterNeighbors;
extern double muFilterNeighbors;
extern double theoreticalPressure;
extern double theoreticalMuA;
extern double theoreticalMuB;

// Evaporation constants
extern double dn; // delta n for evaporation
extern int epos;  // lattice position at which evaporation occurs

// Equations of motion constants
extern double T0;  			// initial theta from write-up
extern double theta;
extern double nA0;
extern double nB0;
extern double nC0;
extern double Amp;
extern double tau;
extern double oneOverTau;   // 1/tau from write-up (relaxation constant)
extern double g;          	// gravitational acceleration term
extern double lambda;     	// friction (F12) coefficient
extern double gammaP; 		// pressure coefficient to control the rate of forcing per time step
extern double gammaMu;		// chemical potential coefficient to control the rate of forcing per time step
extern double kappa;		// coefficient of gradient corrections
extern double aA;
extern double aB;
extern double aC;
extern double aAB;
extern double aAC;
extern double aBC;
extern double bA;
extern double bB;
extern double bC;
extern double b;
extern double vdwInteractionFactor;
extern double interfaceWidth;

// Critical point traits
extern double tcA;
extern double ncA;
extern double pcA;
extern double tcB;
extern double ncB;
extern double pcB;
extern double tcC;
extern double ncC;
extern double pcC;

extern int lnExplosion;
extern int useChemicalPotentialForcingMethod;
extern int useBoundaryConditionsPeriodic;
extern int useMuVdwInteraction;
extern int momentumCorrectionOn;

// GUI
extern int ulreq;
extern int ugreq;
extern int pgreq;
extern int tgreq;
extern int mu1greq; // graph requests
extern int next;
extern int Pause;
extern int run;
extern int done;
extern int Repeat;
extern int iterations;
extern int collectData;
extern int phase_iterations;
extern int total_phase_iterations;


//
// Function Prototypes
//

void (*initializeProfile)(int);
void init();
void initialize();
void initializeRandom(int);
void setInitializeRandom();
void setInitializeStepProfile();

void iteration();
void GUI();
void GetData();

double gradient(double *, int);
double laplace(double *, int);
double filterArrayData(double *, int, int);
void getDensityProfile();
void getPressureProfile();
void getChemicalPotentialProfile();
void logVDWParameters();
void printComponentDensities();
void printComponentMaxMin();
void printChemicalPotentials();
void printVolumeExclusion();

void calculatePressures();

void (*collision)();
void collisionForcingNewChemicalPotentialGradient();
void setCollisionForcingNewChemicalPotentialGradient();


#endif /* LB_D1Q3_2_COMPONENTS_H_ */
