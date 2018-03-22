/**
 * @file LB_GUI.c
 * @author Alexander Wagner (modified by Kent S. Ridl)
 * @date (last modified: 29 November 2017)
 *
 * The module _LB_GUI.c_ controls the layout and function of the program's graphical user interface (GUI).
 */

#include "LB_friction-test.h"


/**
 * @brief Function _GUI_ defines the layout of the program's graphical user interface.
 *
 * This routine initializes the Grapical User Interface (GUI). First we define the fields we want to visualize, and in the second part we structure the menu
 * where we can manipulate variables during the simulation.
 */
void GUI() {
	static int xdimi = XDIM;

	/* Initialize the available graphs to plot */
	SetDefaultColor(1); // black
	DefineGraphN_R("n", &n[0], &xdimi, NULL);
	DefineGraphN_R("n1", &n1[0], &xdimi, NULL);
	SetDefaultColor(2); // red
	DefineGraphN_R("n2", &n2[0], &xdimi, NULL);
	SetDefaultColor(1);
	DefineGraphN_R("u1", &u1[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("u2", &u2[0], &xdimi, NULL);
	SetDefaultColor(1);
	DefineGraphN_R("u-hat1", &uHat1[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("u-hat2", &uHat2[0], &xdimi, NULL);
	SetDefaultColor(1);
	NewGraph();

	DefineGraphN_R("pressure", &pressure[0], &xdimi, NULL);
	NewGraph();

	DefineGraphN_R("mu1", &mu1[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("mu2", &mu2[0], &xdimi, NULL);
	SetDefaultColor(1);
	NewGraph();

	DefineGraphN_R("F1", &F1[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("F2", &F2[0], &xdimi, NULL);
	SetDefaultColor(1);
	DefineGraphN_R("friction1", &friction1[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("friction2", &friction2[0], &xdimi, NULL);
	SetDefaultColor(1);
	DefineGraphN_R("gradMuForce1", &gradMuForce1[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("gradMuForce2", &gradMuForce2[0], &xdimi, NULL);
	SetDefaultColor(1);
	NewGraph();
	/* End of available graphs */

	/* Graphics menu and control */
	StartMenu("D1Q3", 1);

	DefineBool("Start/Stop", &Pause);
	DefineBool("Next Step", &next);
	DefineInt("-- Step Size", &Repeat);
	DefineBool("Run phase_iterations", &run);
	DefineInt("-- phase iterations", &phase_iterations);
	DefineInt("Total Iterations", &iterations);

	SetActiveGraph(0);
	DefineGraph(curve2d_, "Graphs - Density, Velocity");
	SetActiveGraph(1);
	DefineGraph(curve2d_, "Graphs - Pressures");
	SetActiveGraph(2);
	DefineGraph(curve2d_, "Graphs - Chemical Potentials");
	SetActiveGraph(3);
	DefineGraph(curve2d_, "Graphs - Forces");

	StartMenu("Algorithms", 0);
	DefineFunction("Forcing New - grad mu", &setCollisionForcingNewChemicalPotentialGradient);
	DefineInt("-- Grad-Mu forcing type (1-nid,2-log)", &useChemicalPotentialForcingMethod);
	DefineBool("Pressure derived from partition function?", &usePressurePartitionFunction);
	DefineFunction("(Re)Initialize", &initialize);
	EndMenu();

	StartMenu("Constraints", 0);
	DefineDouble("nA0", &nA0);
	DefineDouble("nB0", &nB0);
	DefineDouble("uA0", &uA0);
	DefineDouble("uB0", &uB0);
	DefineDouble("omega", &oneOverTau);
	DefineDouble("*tc-A", &tcA);
	DefineDouble("nc-A", &ncA);
	DefineDouble("pc-A", &pcA);
	DefineDouble("*tc-B", &tcB);
	DefineDouble("*nc-B", &ncB);
	DefineDouble("pc-B", &pcB);
	DefineDouble("Amp", &Amp);
	DefineDouble("theta", &theta);
	DefineDouble("lambda", &lambda);
	DefineDouble("kappa", &kappa);
	DefineDouble("gamma-Mu", &gammaMu);
	DefineDouble("gamma-P", &gammaP);
	DefineDouble("aA", &aA);
	DefineDouble("aB", &aB);
	DefineDouble("aAB", &aAB);
	DefineDouble("bA", &bA);
	DefineDouble("bB", &bB);
	DefineDouble("VDW interaction factor", &vdwInteractionFactor);
	DefineFunction("Initialize - Random", &setInitializeRandom);
//	DefineFunction("Initialize - 2-phase Steps", &setInitializeTheoreticalTwoPhases);
	DefineFunction("(Re)Initialize", &initialize);
	EndMenu();

	DefineBool("EXIT",&done);
	EndMenu();
	/* End of graphics menu */

} // end function GUI()



