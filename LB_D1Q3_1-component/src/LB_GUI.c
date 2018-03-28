/*
 * LB_GUI.c
 *
 *  Created on: Jul 7, 2016
 *      Author: clark
 */


#include "LB_D1Q3_1-component.h"


/**
This routine initializes the Grapical User Interface (GUI). First we
define the fields we want to visualize, and in the second part we
structure the menu where we can manipulate variables during the
simulation.
 */
void GUI() {
	static int xdimi = XDIM;
	static int freeEnergyLocalSize = 1000; //freeEnergyArraySize;

	/* Initialize the available graphs to plot */
	SetDefaultColor(1);
	DefineGraphN_R("n1", &n1[0], &xdimi, NULL);
	//	DefineGraphN_R("u1", &u1g[0], &xdimi, &ugreq);
	DefineGraphN_R("u1", &u1[0], &xdimi, NULL);
	DefineGraphN_R("u-hat1", &uHat1[0], &xdimi, NULL);
	//	DefineGraphN_R("n-theoretical",&tg[0],&xdimi,&tgreq);
	NewGraph();

	SetDefaultColor(2);
	DefineGraphN_R("pressureGradPMethod", &pressureGradPMethod[0], &xdimi, NULL);
	SetDefaultColor(4);
	DefineGraphN_R("pressureGradMuMethod", &pressureGradMuMethod[0], &xdimi, NULL);
	SetDefaultColor(3);
	DefineGraphN_R("pressurePressureMethod", &pressurePressureMethod[0], &xdimi, NULL);
	SetDefaultColor(1);
	DefineGraphN_R("pressure1", &pressure1[0], &xdimi, NULL);
	DefineGraphN_R("pressureTest - no kappa", &pressureTest[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("pressureTest1 - discrete w/ n", &pressureTest1[0], &xdimi, NULL);
	//	SetDefaultColor(3);
	//	DefineGraphN_R("pressureTest2 - discrete w/out n", &pressureTest2[0], &xdimi, NULL);
	SetDefaultColor(3);
	DefineGraphN_R("pressureTest3 - calc n1(l1+l2a)", &pressureTest3[0], &xdimi, NULL);
	//	SetDefaultColor(5);
	//	DefineGraphN_R("pressureTest4 - calc n1(l1+l2b)", &pressureTest4[0], &xdimi, NULL);
	//	SetDefaultColor(2);
	//	DefineGraphN_R("pressureTest5 - calc n2a(l1+l2a)", &pressureTest5[0], &xdimi, NULL);
	//	SetDefaultColor(3);
	//	DefineGraphN_R("pressureTest6 - calc n2b(l1+l2b)", &pressureTest6[0], &xdimi, NULL);
	SetDefaultColor(4);
	DefineGraphN_R("pressureTestZero - e1-e2", &pressureTestZero[0], &xdimi, NULL);
	SetDefaultColor(1);
	DefineGraphN_R("pressureCriticalConstants", &pressureCriticalConstants[0], &xdimi, NULL);
	DefineGraphN_R("pressureVDWConstants", &pressureVDWConstants[0], &xdimi, NULL);
	DefineGraphN_R("pressure CC-VDW", &pressureCCMinusVDW[0], &xdimi, NULL);
	//	DefineGraphN_R("pressureNonIdeal1", &pressureNonIdeal1[0], &xdimi, NULL);
	//	DefineGraphN_R("pressureCorrected1", &pressureCorrected1[0], &xdimi, NULL);
	DefineGraphN_R("p-Alexander", &p[0], &xdimi, NULL);
	DefineGraphN_R("pf-Alexander", &pf[0], &xdimi, NULL);
	//	DefineGraphN_R("pff-Alexander", &pff[0], &xdimi, NULL);
	DefineGraphN_R("PF-Alexander", &PF[0], &xdimi, NULL);
	NewGraph();

	DefineGraphN_R("mu1", &mu1[0], &xdimi, NULL);
	//	DefineGraphN_R("muNonIdeal1", &muNonIdeal1[0], &xdimi, NULL);
	DefineGraphN_R("muCriticalConstants1", &muCriticalConstants1[0], &xdimi, NULL);
	DefineGraphN_R("muVDWConstants1", &muVDWConstants1[0], &xdimi, NULL);
	DefineGraphN_R("mu CC-VDW 1", &muCCMinusVDW1[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("muGradPMethod", &muGradPMethod[0], &xdimi, NULL);
	SetDefaultColor(4);
	DefineGraphN_R("muGradMuMethod", &muGradMuMethod[0], &xdimi, NULL);
	SetDefaultColor(3);
	DefineGraphN_R("muPressureMethod", &muPressureMethod[0], &xdimi, NULL);
	SetDefaultColor(1);
	NewGraph();

	DefineGraphN_R("F1", &F1[0], &xdimi, NULL);
	DefineGraphN_R("F1GradPMethod", &F1GradPMethod[0], &xdimi, NULL);
	DefineGraphN_R("F1GradMuMethod", &F1GradMuMethod[0], &xdimi, NULL);
	DefineGraphN_R("F1 GradP-GradMu", &F1GradPGradMuDifference[0], &xdimi, NULL);
	DefineGraphN_R("GradP", &gradP[0], &xdimi, NULL);
	DefineGraphN_R("GradP - Rho*GradMu", &gradPMinusRhoGradMu[0], &xdimi, NULL);
	DefineGraphN_R("F1_0", &F1_0[0], &xdimi, NULL);
	DefineGraphN_R("F1_1", &F1_1[0], &xdimi, NULL);
	DefineGraphN_R("F1_2", &F1_2[0], &xdimi, NULL);
	NewGraph();

	//	DefineGraphN_R("Free Energy", &freeEnergyArray[0], &freeEnergyLocalSize, NULL);
	//	NewGraph();
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

	StartMenu("Init Variables", 0);
	DefineDouble("n1_liquid", &n1_liquid);
	DefineDouble("n1_gas", &n1_gas);
	DefineFunction("Initialize", &init);
	EndMenu();

	DefineInt("Total Iterations", &iterations);
	DefineInt("Step Size", &Repeat);

	StartMenu("Generate Data", 0);
	DefineFunction("Density Profile", &getDensityProfile);
	DefineFunction("Diffusion Series", &getDiffusionSeries);
	DefineFunction("Rate of Diffusion", &getRateOfDiffusion);
	//	DefineFunction("Free Energy - Single Component Minimization", &calculateFreeEnergySingleComponentMinimization);
	DefineFunction("Phase Diagram (theoretical) - Density vs Temperature (1 component)", &calculatePhaseDiagramRhoVsTempTheoretical);
	DefineFunction("Phase Diagram (theoretical) - Density vs Pressure (1 component)", &calculatePhaseDiagramRhoVsPressureTheoretical);
	DefineFunction("Phase Diagram (theoretical) - Density vs Chemical Potential (1 component)", &calculatePhaseDiagramRhoVsMuTheoretical);
	DefineFunction("Phase Diagram - Density vs Temp", &getPhaseDiagramDensityVsTemp);
	//DefineFunction("Force Data", &getForceData);
	DefineBool("Collect Data", &collectData);
	EndMenu();

	//DefineInt("Phase Iterations", &phase_iterations);

	StartMenu("Algorithms", 0);
	DefineFunction("Forcing New - grad P", &setCollisionForcingNewPressureGradient);
	DefineFunction("Forcing New - grad mu", &setCollisionForcingNewChemicalPotentialGradient);
	DefineFunction("Pressure Method", setCollisionPressureMethod);
	DefineDouble("-- Holdych correction", &pressureMethodCoefficient);
	DefineDouble("-- Laplace term correction", &pressureMethodCorrection);
	DefineFunction("Forcing Basic", &setCollisionForcingBasic);
	DefineFunction("Forcing Alexander", &setCollisionForcingAlexander);
	DefineFunction("Forcing Kyoto", &setCollisionForcingKyoto);
	DefineBool("Pressure - Coupled", &usePressureCoupled);
	DefineBool("-- Critical Constants? (uncoupled only)", &usePressureCriticalParameters);
	DefineBool("Chemical Potential - Critical Constants", &useChemicalPotentialCriticalParameters);
	DefineBool("-- Non-ideal Mu?", &useChemicalPotentialNonIdeal);
	DefineBool("-- ln Explosion", &lnExplosion);
	DefineFunction("Initialize - Random", &setInitializeRandom);
	DefineFunction("Initialize - Steps", &setInitializeSteps);
	DefineFunction("(Re)Initialize", &initialize);
	EndMenu();

	StartMenu("Constraints", 0);
	DefineDouble("omega", &oneOverTau);
	DefineDouble("quench depth", &quenchDepth);
	DefineDouble("tc", &tc);
	DefineDouble("nc", &nc);
	DefineDouble("pc", &pc);
	DefineDouble("Amp", &Amp);
	DefineDouble("n0", &n0);
	DefineDouble("T0", &T0);
	DefineDouble("theta", &theta);
	DefineDouble("g", &g);
	DefineDouble("lambda", &lambda);
	DefineDouble("gamma-P", &gammaP);
	DefineDouble("gamma-Mu", &gammaMu);
	DefineDouble("kappa", &kappa);
	DefineBool("autoKappaGammaMu?", &autoKappaGammaMu);
	DefineDouble("a1", &a1);
	DefineDouble("b1", &b1);
	DefineDouble("pressureMethodCorrection", &pressureMethodCorrection);
	DefineDouble("dummy", &dummy);
	DefineDouble("dn", &dn);
	DefineInt("epos", &epos);
	DefineInt("wall", &wall);
	DefineFunction("(Re)Initialize", &initialize);
	EndMenu();

	StartMenu("Stuff", 0);
	DefineDouble("rho1", &rho1);
	DefineDouble("excludedVolume1", &excludedVolume1);
	DefineBool("Boundary Conditions - Periodic", &useBoundaryConditionsPeriodic);
	DefineFunction("Generate free energy surface data", &generateFreeEnergySlices);
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
		double A = n0*3*g*XDIM/(1-exp(-3*g*XDIM));
		for (i = 0; i < XDIM; i++) {
			tg[i] = A*exp(-3*g*i);
		}
		tgreq=0;
	}

} // end function GetData()



