/**
 * @file LB_GUI.c
 * @author Alexander Wagner (modified by Kent S. Ridl)
 * @date (last modified: 29 November 2017)
 *
 * The module _LB_GUI.c_ controls the layout and function of the program's graphical user interface (GUI).
 */

#include "LB_D1Q3_2-components.h"


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
	DefineGraphN_R("u-hat1", &uHat1[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("u-hat2", &uHat2[0], &xdimi, NULL);
	SetDefaultColor(1);
	NewGraph();

	DefineGraphN_R("pressure", &pressure[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("theoretical pressure", &theoreticalPressureArray[0], &xdimi, NULL);
	SetDefaultColor(1);
	NewGraph();

	DefineGraphN_R("mu1", &mu1[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("mu2", &mu2[0], &xdimi, NULL);
	SetDefaultColor(1);
	DefineGraphN_R("theoreticalMuA", &theoreticalMuAArray[0], &xdimi, NULL);
	SetDefaultColor(2);
	DefineGraphN_R("theoreticalMuB", &theoreticalMuBArray[0], &xdimi, NULL);
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

	StartMenu("Initialization", 0);
	DefineDouble("nA0", &nA0);
	DefineDouble("nB0", &nB0);
	DefineDouble("theoreticalDensityA1", &theoreticalDensityA1);
	DefineDouble("theoreticalDensityA2", &theoreticalDensityA2);
	DefineDouble("theoreticalDensityA3", &theoreticalDensityA3);
	DefineDouble("theoreticalDensityB1", &theoreticalDensityB1);
	DefineDouble("theoreticalDensityB2", &theoreticalDensityB2);
	DefineDouble("theoreticalDensityB3", &theoreticalDensityB3);
	DefineDouble("theoreticalVolume1", &theoreticalVolume1);
	DefineDouble("theoreticalVolume2", &theoreticalVolume2);
	DefineDouble("theoreticalVolume3", &theoreticalVolume3);
	DefineBool("-- Use theoretical densities?", &useTheoreticalDensities);
	DefineBool("-- Use theoretical volumina?", &useTheoreticalVolumina);
	DefineFunction("Initialize - Random", &setInitializeRandom);
	DefineFunction("Initialize - 2-phase Steps", &setInitializeTheoreticalTwoPhases);
	DefineBool("-- use theoretical phase 1", &useTheoreticalPhase1);
	DefineBool("-- use theoretical phase 2", &useTheoreticalPhase2);
	DefineBool("-- use theoretical phase 3", &useTheoreticalPhase3);
	DefineFunction("Initialize - 3-phase Steps", &setInitializeTheoreticalThreePhases);
	DefineBool("-- suppress theoretical phase 1", &suppressTheoreticalPhase1);
	DefineBool("-- suppress theoretical phase 2", &suppressTheoreticalPhase2);
	DefineFunction("Initialize - 4-phase Steps", &setInitializeTheoreticalFourPhases);
	DefineBool("-- Isolate 4th Phase?", &isolateFourthPhase);
	EndMenu();

	StartMenu("Generate Data", 0);
	//	DefineFunction("Free Energy Minimization - 2 components, 2 phases", &calculatePhaseDiagramRhoAVsRhoBTwoPhasesTheoretical);
	DefineFunction("Free Energy Minimization - 2 components, 3 phases", &calculatePhaseDiagramRhoAVsRhoBThreePhasesTheoretical);
	DefineFunction("Generate LB Phase Diagram", &generatePhaseDiagramRhoAVsRhoB);
	DefineBool("-- tanh (on) or random (off) profile?", &useStepOrRandomInitialization);
	DefineBool("- LB Phase Diagram - 2-phase Region", &setPhaseDiagramTwoPhaseRegion);
	DefineBool("- LB Phase Diagram - 3-phase Region", &setPhaseDiagramThreePhaseRegion);
	DefineBool("-- Check 2-phase metastable points?", &checkTwoPhaseMetastableLBPoints);
	DefineBool("-- Check 3-phase points?", &checkThreePhaseLBPoints);
	DefineFunction("Pressure/Mu Deviations", &generateMinimizedPressureAndMuDeviations);
	DefineFunction("Print Component Densities", &printComponentDensities);
	DefineFunction("Print max/min densities", &printComponentMaxMin);
	DefineFunction("Save density profiles to file", &getDensityProfile);
	DefineFunction("Save pressure profile to file", &getPressureProfile);
	DefineFunction("Save chemical potential profiles to file", &getChemicalPotentialProfile);
	DefineFunction("Eigenvalue Map", &generateEigenvalueMap);
	EndMenu();

	StartMenu("Algorithms", 0);
	DefineFunction("Forcing New: grad-Mu", &setCollisionForcingNewChemicalPotentialGradient);
	DefineInt("-- grad-Mu forcing type (1-nid,2-log)", &useChemicalPotentialForcingMethod);
	DefineBool("Pressure derived from partition function?", &usePressurePartitionFunction);
	DefineFunction("(Re)Initialize", &initialize);
	EndMenu();

	StartMenu("Constraints", 0);
	DefineDouble("** tc-A **", &tcA);
	DefineDouble("** tc-B **", &tcB);
	DefineDouble("** nc-B **", &ncB);
	DefineDouble("** VDW interaction factor **", &vdwInteractionFactor);
	DefineBool("-- Use in mu calculations?", &useMuVdwInteraction);
	DefineDouble("nc-A", &ncA);
	DefineDouble("pc-A", &pcA);
	DefineDouble("pc-B", &pcB);
	DefineDouble("Amp", &Amp);
	DefineDouble("omega", &oneOverTau);
	DefineDouble("theta", &theta);
	//	DefineDouble("g", &g);
	DefineDouble("lambda", &lambda);
	DefineDouble("interface width", &interfaceWidth);
	DefineBool("-- automatic interface width?", &determineInterfaceWidthAutomatically);
	DefineBool("-- override minimum width (2.0)?", &overrideMinimumInterfaceWidth);
	DefineDouble("kappa", &kappa);
	DefineBool("-- determineKappaAutomatically", &determineKappaAutomatically);
	DefineDouble("kappaFactor", &kappaFactor);
	DefineBool("-- kappa factor determined?", &kappaFactorDetermined);
	DefineDouble("gamma-Mu", &gammaMu);
	DefineDouble("gammaFactor", &gammaFactor);
	DefineDouble("gamma-P", &gammaP);
	DefineDouble("aA", &aA);
	DefineDouble("aB", &aB);
	DefineDouble("aAB", &aAB);
	DefineDouble("bA", &bA);
	DefineDouble("bB", &bB);
	DefineFunction("(Re)Initialize", &initialize);
	EndMenu();

	StartMenu("Stuff", 0);
	DefineDouble("minimizationParticlesStepSize", &minimizationParticlesStepSize);
	DefineDouble("-- step size divisor for minimization", &minimizationStepFactor);
	DefineFunction("Minimize single point", &minimizeSinglePoint);
	DefineDouble("-- A-coordinate", &singlePointRhoA);
	DefineDouble("-- B-coordinate", &singlePointRhoB);
	//	DefineFunction("Generate free energy differences", &generateFreeEnergyData);
	//	DefineFunction("Generate free energy slices/surface", &generateFreeEnergySlices);
	DefineFunction("Explore VDW shield region", &exploreVDWShieldRegion);
	DefineBool("Periodic boundary conditions?", &useBoundaryConditionsPeriodic);
	DefineBool("Apply momentum correction?", &applyMomentumCorrection);
	DefineFunction("Generate Phase Diagram Binodal", &generatePhaseDiagramBinodalLines);
	DefineBool("-- sort A binodal by B (on) or A (off)?", &sortABinodalByIncreasingB);
	DefineBool("-- sort B binodal by A (on) or B (off)?", &sortBBinodalByIncreasingA);
	DefineFunction("Generate Gnuplot Scripts", &gnuplotTwoComponentThreePhase);
	DefineFunction("Generate Publication Scripts", &gnuplotPublicationGraphics);
	DefineFunction("Generate interface width fit data", &generateInterfaceWidthFitValues);
	DefineFunction("Sort metastable minimization results", &sortMetastableMinimizationResults);
	DefineFunction("Print minimization parent of point (below):", &printMinimizationParentPoint);
	DefineDouble("-- child rhoA", &childRhoA);
	DefineDouble("-- child rhoB", &childRhoB);
	DefineFunction("Print LB parent of point (above):", &printLBParentPoint);
	DefineFunction("Sort metastable LB results", &sortMetastableLBResults);
	DefineDouble("--metastableThreshold", &metastableThreshold);
	DefineDouble("--2 vs. 3 phase factor", &metastableThresholdFactor);
	DefineDouble("nAIntegrated", &nAIntegrated);
	DefineDouble("nBIntegrated", &nBIntegrated);
	DefineFunction("van K and Scott parameters", &calculateKandSParameters);
	EndMenu();

	DefineBool("EXIT",&done);
	EndMenu();
	/* End of graphics menu */

} // end function GUI()


///**
// * @brief Function _GetData_ allows you to calculate quantities to be displayed in the graphics routine, even if you do not need them during each iteration.
// */
//void GetData(){
//	if (tgreq) {
//		double A = nA0*3*g*XDIM/(1-exp(-3*g*XDIM));
//		for (int i = 0; i < XDIM; i++) {
//			tg[i] = A*exp(-3*g*i);
//		}
//		tgreq=0;
//	}
//} // end function GetData()



