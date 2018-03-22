/*
 * LB_GUI.c
 *
 *  Created on: Jul 7, 2016
 *      Author: clark
 */


#include "LB_D1Q3_3-components.h"


/**
This routine initializes the Grapical User Interface (GUI). First we
define the fields we want to visualize, and in the second part we
structure the menu where we can manipulate variables during the
simulation.
 */
void GUI() {
	static int xdimi = XDIM;

	/* Initialize the available graphs to plot */
	SetDefaultColor(1); // black
	DefineGraphN_R("n", &n[0], &xdimi, NULL);
	DefineGraphN_R("nA", &nA[0], &xdimi, NULL);
	SetDefaultColor(2); // red
	DefineGraphN_R("nB", &nB[0], &xdimi, NULL);
	SetDefaultColor(3);
	DefineGraphN_R("nC", &nC[0], &xdimi, NULL);
	SetDefaultColor(1);
	DefineGraphN_R("u", &u[0], &xdimi, NULL);
	DefineGraphN_R("uHat-A", &uHatA[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("uHat-B", &uHatB[0], &xdimi, NULL);
	SetDefaultColor(3);
	DefineGraphN_R("uHat-C", &uHatC[0], &xdimi, &ugreq);
	SetDefaultColor(1);
	NewGraph();

	DefineGraphN_R("pressure", &pressure[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("pressure-filtered", &pressureFiltered[0], &xdimi, NULL);
	SetDefaultColor(1);
	//	SetDefaultColor(2); // red
	//	DefineGraphN_R("theoretical pressure", &theoreticalPressureArray[0], &xdimi, NULL);
	//	SetDefaultColor(1); // black
	NewGraph();

	DefineGraphN_R("muA", &muA[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("muB", &muB[0], &xdimi, NULL);
	SetDefaultColor(3);
	DefineGraphN_R("muC", &muC[0], &xdimi, NULL);
	SetDefaultColor(1);
	DefineGraphN_R("muA-filtered", &muAFiltered[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("muB-filtered", &muBFiltered[0], &xdimi, NULL);
	SetDefaultColor(3);
	DefineGraphN_R("muC-filtered", &muCFiltered[0], &xdimi, NULL);
	//	DefineGraphN_R("theoreticalMuA", &theoreticalMuAArray[0], &xdimi, NULL);
	//	SetDefaultColor(2);
	//	DefineGraphN_R("theoreticalMuB", &theoreticalMuBArray[0], &xdimi, NULL);
	//	SetDefaultColor(3);
	//	DefineGraphN_R("theoreticalMuC", &theoreticalMuCArray[0], &xdimi, NULL);
	//	SetDefaultColor(1);
	NewGraph();

	DefineGraphN_R("Force A", &forceA[0], &xdimi, NULL);
	DefineGraphN_R("Force B", &forceB[0], &xdimi, NULL);
	DefineGraphN_R("Force C", &forceC[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("frictionA", &frictionA[0], &xdimi, NULL);
	DefineGraphN_R("frictionB", &frictionB[0], &xdimi, NULL);
	DefineGraphN_R("frictionC", &frictionC[0], &xdimi, NULL);
	SetDefaultColor(1);
	NewGraph();

	/* End of available graphs */

	/* Graphics menu and control */
	StartMenu("D1Q3", 1);

	SetActiveGraph(0);
	DefineGraph(curve2d_, "Graphs - Density, Velocity");
	SetActiveGraph(1);
	DefineGraph(curve2d_, "Graphs - Pressures");
	SetActiveGraph(2);
	DefineGraph(curve2d_, "Graphs - Chemical Potentials");
	SetActiveGraph(3);
	DefineGraph(curve2d_, "Graphs - Forces");
	SetActiveGraph(4);
	DefineGraph(curve2d_, "Graphs - Phase Diagrams");
	DefineBool("Start/Stop", &Pause);
	DefineBool("Next Step", &next);
	DefineBool("Run phase_iterations", &run);

	StartMenu("Initialization", 0);
	DefineDouble("nA0", &nA0);
	DefineDouble("nB0", &nB0);
	DefineDouble("nC0", &nC0);
	DefineDouble("theoreticalDensityA1", &theoreticalDensityA1);
	DefineDouble("theoreticalDensityA2", &theoreticalDensityA2);
	DefineDouble("theoreticalDensityA3", &theoreticalDensityA3);
	DefineDouble("theoreticalDensityA4", &theoreticalDensityA4);
	DefineDouble("theoreticalDensityB1", &theoreticalDensityB1);
	DefineDouble("theoreticalDensityB2", &theoreticalDensityB2);
	DefineDouble("theoreticalDensityB3", &theoreticalDensityB3);
	DefineDouble("theoreticalDensityB4", &theoreticalDensityB4);
	DefineDouble("theoreticalDensityC1", &theoreticalDensityC1);
	DefineDouble("theoreticalDensityC2", &theoreticalDensityC2);
	DefineDouble("theoreticalDensityC3", &theoreticalDensityC3);
	DefineDouble("theoreticalDensityC4", &theoreticalDensityC4);
	DefineFunction("Initialize - Random", &setInitializeRandom);
	DefineFunction("Initialize - Steps", &setInitializeStepProfile);
	EndMenu();

	DefineInt("Total Iterations", &iterations);
	DefineInt("Step Size", &Repeat);

	StartMenu("Generate Data", 0);
	DefineFunction("Print Component Densities", &printComponentDensities);
	DefineFunction("Print max/min densities", &printComponentMaxMin);
	DefineFunction("Print chemical potentials", &printChemicalPotentials);
	DefineFunction("Print volume exclusion", &printVolumeExclusion);
	DefineFunction("Save density profiles to file", &getDensityProfile);
	DefineFunction("Save pressure profile to file", &getPressureProfile);
	DefineFunction("Save chemical potential profiles to file", &getChemicalPotentialProfile);
	DefineFunction("Save simulation VDW parameters to file", &logVDWParameters);
	EndMenu();

	StartMenu("Algorithms", 0);
	DefineFunction("Forcing New - grad mu", &setCollisionForcingNewChemicalPotentialGradient);
	DefineInt("-- Grad-Mu forcing type (1-nid,2-full)", &useChemicalPotentialForcingMethod);
	DefineFunction("(Re)Initialize", &initialize);
	EndMenu();

	StartMenu("Constraints", 0);
	DefineDouble("omega", &oneOverTau);
	DefineDouble("*tc-A", &tcA);
	DefineDouble("nc-A", &ncA);
	DefineDouble("pc-A", &pcA);
	DefineDouble("*tc-B", &tcB);
	DefineDouble("*nc-B", &ncB);
	DefineDouble("pc-B", &pcB);
	DefineDouble("*tc-C", &tcC);
	DefineDouble("*nc-C", &ncC);
	DefineDouble("pc-C", &pcC);
	DefineDouble("Amp", &Amp);
	DefineDouble("T0", &T0);
	DefineDouble("theta", &theta);
	DefineDouble("g", &g);
	DefineDouble("lambda", &lambda);
	DefineDouble("gamma-P", &gammaP);
	DefineDouble("gamma-Mu", &gammaMu);
	DefineDouble("kappa", &kappa);
	DefineDouble("aA", &aA);
	DefineDouble("aB", &aB);
	DefineDouble("aC", &aC);
	DefineDouble("aAB", &aAB);
	DefineDouble("aAC", &aAC);
	DefineDouble("aBC", &aBC);
	DefineDouble("bA", &bA);
	DefineDouble("bB", &bB);
	DefineDouble("bC", &bC);
	DefineDouble("VDW interaction factor", &vdwInteractionFactor);
	DefineBool("-- Use in mu calculations?", &useMuVdwInteraction);
	DefineDouble("Interface width (% of XDIM)", &interfaceWidth);
	DefineFunction("(Re)Initialize", &initialize);
	EndMenu();

	StartMenu("Stuff", 0);
	DefineInt("phase iterations", &phase_iterations);
	DefineBool("Boundary Conditions - Periodic", &useBoundaryConditionsPeriodic);
	DefineBool("Momentum Correction", &momentumCorrectionOn);
	DefineDouble("Pressure filter neighbors", &pressureFilterNeighbors);
	DefineDouble("Mu filter neighbors", &muFilterNeighbors);
	EndMenu();

	DefineBool("EXIT",&done);

	EndMenu();
	/* End of graphics menu */

} // end function GUI()


/**
This routine allows you to calculate quantities to be displayed
in the graphics routine, even if you do not need them during each
iteration.
 */
void GetData(){
	int i;

	if (tgreq) {
		double A = nA0*3*g*XDIM/(1-exp(-3*g*XDIM));
		for (i = 0; i < XDIM; i++) {
			tg[i] = A*exp(-3*g*i);
		}
		tgreq=0;
	}

} // end function GetData()



