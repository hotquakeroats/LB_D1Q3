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

#include "mygraph.h"


#define XDIM 101 // lattice length
// Even lattice: mu->10^-8 pinched; p->10^-6 no ripple (walls) mu->10^-8 uniform; p->10^-6 ripple (periodic)
// Odd lattice:  mu->10^-6 pinched; p->10^-5 big ripple (walls) mu->10^-9 interfaces; p->10^-6 no ripple (periodic)

#define PB(n) (pc*((n)/3./(1-(n)/(3.*nc))-9./8. *(n)*(n)*tc/nc)) //theta/nc)) //
//#define PB(n) ((n)/(3.-(n))-(9./8.)*tc*(n)*(n))
#define PI(n,dn,ddn) (-pc*kappa*((n)*(ddn)-0.5*(dn)*(dn)))
#define P(n,dn,ddn) (PB(n)+PI(n,dn,ddn))
#define Pni(n,dn,ddn) (P(n,dn,ddn)-n/3.)

#define MAXSTRLEN 255 ///< Maximum length of a file name
#define isEqual(a, b, threshold) ( (fabs((a)-(b)) < (threshold)) ? 1 : 0 ) ///< Equality tester for all the doubles out there.


// Equilibrium distribution of each component velocities
double f1_0[XDIM]; // the static vector
double f1_1[XDIM]; // +1 lattice space
double f1_2[XDIM]; // -1 lattice space

double psi1[XDIM];

// Number densities of each component
double n1[XDIM]; // 1st component mass/density

double *freeEnergyArray;

// Component and bulk velocities
double u1[XDIM];
double uHat1[XDIM];

// Component thermodynamic properties
double mu1[XDIM]; // chemical potential
double muNonIdeal1[XDIM];
double muCriticalConstants1[XDIM];
double muVDWConstants1[XDIM];
double muCCMinusVDW1[XDIM];
double muGradPMethod[XDIM];
double muGradMuMethod[XDIM];
double muPressureMethod[XDIM];

double pressure[XDIM];
double pressure1[XDIM];
double pressureNonIdeal1[XDIM];
double pressureCorrected1[XDIM];
double pressureGradPMethod[XDIM];
double pressureGradMuMethod[XDIM];
double pressurePressureMethod[XDIM];
double pressureCriticalConstants[XDIM];
double pressureVDWConstants[XDIM];
double pressureCCMinusVDW[XDIM];
double correctionsPressure1[XDIM];
double correctionsPressureDisplay1[XDIM];
double pressureTest[XDIM];
double pressureTestZero[XDIM];
double pressureTest1[XDIM];
double pressureTest2[XDIM];
double pressureTest3[XDIM];
double pressureTest4[XDIM];
double pressureTest5[XDIM];
double pressureTest6[XDIM];

double dni[XDIM];
double ddni[XDIM];
double dpi[XDIM];
double ddpi[XDIM];
double p[XDIM];
double pni[XDIM];
double pf[XDIM];
//double pff[XDIM];
double PF[XDIM];

// Friction forces for each component velocity
double F1_0[XDIM];
double F1_1[XDIM];
double F1_2[XDIM];

// Forces on each component; potential and friction
double F1[XDIM];
double F1GradPMethod[XDIM];
double F1GradMuMethod[XDIM];
double F1GradPGradMuDifference[XDIM];
double Forces[XDIM];

double A[XDIM];		// pressure method change to the density distributions

double gradP[XDIM];
double rhoGradMu[XDIM];
double gradPMinusRhoGradMu[XDIM];

// Graph/GUI variables
double u1g[XDIM];
//double pg1[XDIM];
double tg[XDIM];
//double mu1g[xdim];


//
// Externs
//

// Initial densities of each component
extern double n1_liquid;
extern double n1_gas;

// Total amount of each component
extern double rho1;
extern double excludedVolume1;

extern int freeEnergyArraySize;

// Evaporation constants
extern double dn; // delta n for evaporation
extern int epos;  // lattice position at which evaporation occurs

// Equations of motion constants
extern double T0;  			// initial theta from write-up
extern double theta;
extern double n0;
extern double Amp;
extern double tau;
extern double oneOverTau;   // 1/tau from write-up (relaxation constant)
extern double g;          	// gravitational acceleration term
extern double lambda;     	// friction (F12) coefficient
extern double gammaP; 		// pressure coefficient to control the rate of forcing per time step
extern double gammaMu;		// chemical potential coefficient to control the rate of forcing per time step
extern double kappa;		// coefficient of gradient corrections
extern double a1;
extern double b1;
extern double dummy;
extern double quenchDepth;	// depth of tc/theta ratio for calculating phase diagrams
extern double pressureMethodCorrection;		// coefficient of laplace correction in the pressure method
extern double pressureMethodCoefficient;	// coefficient of the Holdych correction in the pressure method

// Critical point traits
extern double tc;
extern double nc;
extern double pc;

extern int usePressureCoupled;
extern int usePressureCriticalParameters;	// chooses between tc,pc,nc and a,b VDW constants
extern int useChemicalPotentialCriticalParameters;
extern int useChemicalPotentialsCoupled;  // chooses between independent or coupled mu's
extern int useChemicalPotentialNonIdeal;
extern int lnExplosion;
extern int useBoundaryConditionsPeriodic;

// GUI
extern int ulreq;
extern int ugreq;
extern int pgreq;
extern int tgreq;
extern int mu1greq; // graph requests
extern int next;
extern int Pause;
extern int done;
extern int Repeat;
extern int iterations;
extern int collectData;
extern int wall;
extern int phase_iterations;


//
// Function Prototypes
//

void (*initializeProfile)(int);
void init();
void initialize();
void initializeSteps(int);
void initializeRandom(int);
void setInitializeRandom();
void setInitializeSteps();

void iteration();

void GUI();
void GetData();

double gradient(double *, int);
double laplace(double *,int);
void compareGradPRhoGradMu();
void getPhaseDiagramDensityVsTemp();
void getDensityProfile();
void getDiffusionSeries();
void getRateOfDiffusion();
void calculatePhaseDiagramRhoVsTempTheoretical();
void calculatePhaseDiagramRhoVsPressureTheoretical();
void calculatePhaseDiagramRhoVsMuTheoretical();
void calculateFreeEnergySingleComponentMinimization();
void calculateFreeEnergyABOrdering();
void calculateFreeEnergyFourPhaseOrdering();
void generateFreeEnergySlices();
//void getForceData();

void calculatePressures();

void (*collision)();
void collisionIdeal();
void collisionForcingNewPressureGradient();
void collisionForcingNewChemicalPotentialGradient();
void collisionForcingAlexander();
void collisionPressureMethod();
void setCollisionForcingNewPressureGradient();
void setCollisionForcingNewChemicalPotentialGradient();
void setCollisionForcingOriginal();
void setCollisionForcingBasic();
void setCollisionForcingAlexander();
void setCollisionForcingKyoto();
void setCollisionPressureMethod();
void setCollisionIdeal();

void minimize();


#endif /* LB_D1Q3_2_COMPONENTS_H_ */
