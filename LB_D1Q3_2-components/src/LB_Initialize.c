/**
 * @file LB_Initialize.c
 * @author Kent S. Ridl
 * @date 7 July 2016 (last modified: 28 November 2017)
 *
 * The module _LB_Initialize.c_ contains the code used to initialize lattice Boltzmann simulations.  This includes fetching the theoretical results from
 * the free energy minimization for a test point, setting the near-equilibrium phases for 2-phase and 3-phase behavior, and optimizing simulation
 * parameters for the anticipated interface width.
 */

#include "LB_D1Q3_2-components.h"

//TODO: add ability to initialize 3-phase point with only 3 interfaces (a-vapor, b-vapor, a-b) instead of 4
//TODO: may need to break kappa symmetry
//TODO: implement a choice to do the full golden search for an optimal value or to simply find one that works and go (to save time)

void (*initializeProfile)(int); ///< function pointer to choose the method to initialize the lattice density profile (random, 2-phase, 3-phase)


/**
 * @brief Function _fitInterfaceWidth_ determines the near-equilibrium interface width of a density profile for a given value of kappa.
 *
 * LB simulations first initialize a tanh density profile by estimating interface width according to single-component theory (see reference); however,
 * this estimate often overbounds the actual interface widths realized in multicomponent LB simulations.  This function modifies the initialized interface
 * width to better fit the chosen value of kappa by first iterating a LB simulation to a near-equilbrium state (50,000 time steps) so that the width
 * can be numerically measured.
 *
 * @see A. Wagner and C. Pooley, Physical Review E 76, 045702 (2007)
 */
void fitInterfaceWidth() {
	double evA = 0.0, evB = 0.0; // dummy variables

	// Run a single simulation at a user-defined kappa; fit a tanh profile width to it if it's stable
	if (!calculateMinimizationPath(nA0, nB0, &evA, &evB)) {
		printf("The point (%f,%f) doesn't have a negative eigenvalue and isn't expected to phase separate.\n"
				"Unable to fit an appropriate interface width.\n\n", nA0, nB0);
	}
	else {
		if (theoreticalPhases == 2) useTwoOrThreePhaseInitialization = 2;
		else if (theoreticalPhases == 3) useTwoOrThreePhaseInitialization = 3;
		kappaFactor = 1.0;

		printf("\nFitting interface width for (%f,%f)...\n", nA0, nB0);

		kappaFactorDetermined = 1; // set this flag here to avoid an infinite initialization loop
		setLBInitializationProfile();
		for (int i = 0; i < phase_iterations; i++) iteration(); // iterate the simulation for the near-eq baseline density profile

		if (!stableSimulation(n1) || !stableSimulation(n2)) {
			printf("Unstable simulation!!... width=%f, kappa=%f, kappaFactor=%f, gammaMu=%f, gammaFactor=%f...\n\n",
					interfaceWidth, kappa, kappaFactor, gammaMu, gammaFactor);
			interfaceWidthForGivenKappa = 10.0 * XDIM;
			kappaFactorDetermined = 0;
		}
		else { // golden section search on the interface widths used for subsequent initializations
			double baseDensityProfileA[XDIM], baseDensityProfileB[XDIM];
			double aSSE = 0.0, bSSE = 0.0, ssePoint1 = 1000.0, ssePoint2 = 1000.0;
			double goldenRatio = 0.5 * (1.0 + sqrt(5));
			double widthIntervalStart = 1.0;
			double widthIntervalEnd = XDIM / 4;
			double widthIntervalPoint1 = widthIntervalEnd - (widthIntervalEnd-widthIntervalStart)/goldenRatio;
			double widthIntervalPoint2 = widthIntervalStart + (widthIntervalEnd-widthIntervalStart)/goldenRatio;

			for (int i = 0; i < XDIM; i++) { // keep this initial profile for the chosen kappa value constant
				baseDensityProfileA[i] = n1[i];
				baseDensityProfileB[i] = n2[i];
			}

			interfaceWidth = widthIntervalPoint1;
			for (int i = 0; i < XDIM; i++) initializeProfile(i); // just want to initialize the tanh profiles; don't need the setInitializationX functions
			aSSE = 0.0;
			bSSE = 0.0;
			for (int i = 0; i < XDIM; i++) {
				aSSE += pow(n1[i]-baseDensityProfileA[i], 2); // measure new interface width initialization against near-eq profile
				bSSE += pow(n2[i]-baseDensityProfileB[i], 2);
			}
			ssePoint1 = (nA0 > 1e-10 ? aSSE : 0.0) + (nB0 > 1e-10 ? bSSE : 0.0); // sum SSE for both A and B density profiles

			interfaceWidth = widthIntervalPoint2;
			for (int i = 0; i < XDIM; i++) initializeProfile(i); // just want to initialize the tanh profiles; don't need the setInitializationX functions
			aSSE = 0.0;
			bSSE = 0.0;
			for (int i = 0; i < XDIM; i++) {
				aSSE += pow(n1[i]-baseDensityProfileA[i], 2); // measure new interface width initialization against near-eq profile
				bSSE += pow(n2[i]-baseDensityProfileB[i], 2);
			}
			ssePoint2 = (nA0 > 1e-10 ? aSSE : 0.0) + (nB0 > 1e-10 ? bSSE : 0.0); // sum SSE for both A and B density profiles

			int j = 1;
			double widthIntervalThreshold = 1e-3;
			while (fabs(widthIntervalEnd-widthIntervalStart) > widthIntervalThreshold) { // golden section search loop

				printf("Interface fit iteration %i: widthIntervalStart=%f, widthIntervalPoint1=%f, widthIntervalPoint2=%f, widthIntervalEnd=%f\n",
						j, widthIntervalStart, widthIntervalPoint1, widthIntervalPoint2, widthIntervalEnd);
				printf("-- ssePoint1=%f ssePoint2=%f\n", ssePoint1, ssePoint2);

				if (ssePoint1 < ssePoint2) {
					widthIntervalEnd = widthIntervalPoint2;
					widthIntervalPoint2 = widthIntervalPoint1;
					widthIntervalPoint1 = widthIntervalEnd - (widthIntervalEnd-widthIntervalStart)/goldenRatio;
					interfaceWidth = widthIntervalPoint1;
					for (int i = 0; i < XDIM; i++) initializeProfile(i); // bypass the full initialization routine; just want to initialize the tanh profiles
					aSSE = 0.0;
					bSSE = 0.0;
					for (int i = 0; i < XDIM; i++) {
						aSSE += pow(n1[i]-baseDensityProfileA[i], 2);
						bSSE += pow(n2[i]-baseDensityProfileB[i], 2);
					}
					ssePoint2 = ssePoint1;
					ssePoint1 = (nA0 > 1e-10 ? aSSE : 0.0) + (nB0 > 1e-10 ? bSSE : 0.0);
				}
				else {
					widthIntervalStart = widthIntervalPoint1;
					widthIntervalPoint1 = widthIntervalPoint2;
					widthIntervalPoint2 = widthIntervalStart + (widthIntervalEnd-widthIntervalStart)/goldenRatio;
					interfaceWidth = widthIntervalPoint2;
					for (int i = 0; i < XDIM; i++) initializeProfile(i); // bypass the full initialization routine; just want to initialize the tanh profiles
					aSSE = 0.0;
					bSSE = 0.0;
					for (int i = 0; i < XDIM; i++) {
						aSSE += pow(n1[i]-baseDensityProfileA[i], 2);
						bSSE += pow(n2[i]-baseDensityProfileB[i], 2);
					}
					ssePoint1 = ssePoint2;
					ssePoint2 = (nA0 > 1e-10 ? aSSE : 0.0) + (nB0 > 1e-10 ? bSSE : 0.0);
				}
				j++; // golden search iteration counter
			} // end while

			interfaceWidthForGivenKappa = 0.5 * (widthIntervalStart + widthIntervalEnd);
			fitInterfaceWidthToKappa = 2; // indicate ok to reinitialize with the determined interface width for the given kappa
			printf("determined interface width=%f for (%f,%f)\n\n", interfaceWidthForGivenKappa, nA0, nB0);
		} // end else (for a stable simulation)
	} // end else (for a nA0,nB0 point that is expected to minimize)
} // end function fitInterfaceWidth()


/**
 * @brief Function _fitGammaFactor_ optimizes the value of gammaMu (p0 in other literature) to minimize the number of time steps needed to reach equilibrium.
 *
 * The parameter gammaMu is a fractional application of the chemical potential forcing term for each time step that is applied to increase the stability of
 * LB simulations.  LB simulations are initialized by default to use a value of gammaMu determined according to single-component theory (see reference).  This
 * is sufficient for most simulations; however, it is prone to failure in multicomponent simulations at high density ratios.  This function performs a golden
 * section search over several LB simulations in an attempt to find a maximum value of gammaMu to stabilize the desired simulation.
 *
 * @note The value of gammaMu determined by single-component theory is divided by a correction factor.  The maximum value of this correction is determined by
 * the global value _endGammaFactor_ and can be modified in the GUI.
 * @note A 3% safety factor is added to the correction value identified by the search algorithm.
 * @note This function will only reduce the value of gammaMu through the correction factor.  If single-component theory provides a value of gammaMu that is
 * also stable for a multicomponent simulation, that value is used even though it may not be a maximum stable value.
 *
 * @see A. Wagner and C. Pooley, Physical Review E 76, 045702 (2007)
 */
int fitGammaFactor() {
	double goldenRatio = 0.5 * (1.0 + sqrt(5));
	double gammaIntervalStart = 1.0;
	double gammaIntervalEnd = endGammaFactor;
	double gammaIntervalPoint1 = gammaIntervalEnd - (gammaIntervalEnd-gammaIntervalStart)/goldenRatio;
	double gammaIntervalPoint2 = gammaIntervalStart + (gammaIntervalEnd-gammaIntervalStart)/goldenRatio;

	printf("Fitting gammaMu value for (%f,%f) at kappa=%f, kappaFactor=%f...\n\n", nA0, nB0, kappa, kappaFactor);

	// Test extreme value first... assume entire interval is not stable if this doesn't work
	gammaFactor = gammaIntervalEnd;
	fitInterfaceWidthToKappa = 1;
	kappaFactorDetermined = 1; // set this flag here to avoid an infinite initialization loop
	setLBInitializationProfile();
	for (int i = 0; i < phase_iterations; i++) iteration();

	if (!stableSimulation(n1) || !stableSimulation(n2)) {
		printf("Extreme gammaMu value isn't stable... assuming entire interval is not stable.\n");
		return 0;
	}
	else { // something in the interval is valid; find the minimum gammaFactor that has a stable simulation (maximizes gammaMu)
		int stablePoint1 = 1, stablePoint2 = 1;

		gammaFactor = gammaIntervalPoint1;
		fitInterfaceWidthToKappa = 1;
		setLBInitializationProfile();
		for (int i = 0; i < phase_iterations; i++) iteration();
		if (!stableSimulation(n1) || !stableSimulation(n2)) stablePoint1 = 0;

		gammaFactor = gammaIntervalPoint2;
		fitInterfaceWidthToKappa = 1;
		setLBInitializationProfile();
		for (int i = 0; i < phase_iterations; i++) iteration();
		if (!stableSimulation(n1) || !stableSimulation(n2)) stablePoint2 = 0;

		int j = 1;
		double gammaIntervalThreshold = 0.1;
		while (fabs(gammaIntervalEnd-gammaIntervalStart) > gammaIntervalThreshold) {

			printf("GammaMu fit iteration %i: gammaIntervalStart=%f, gammaIntervalPoint1=%f, gammaIntervalPoint2=%f, gammaIntervalEnd=%f\n",
					j, gammaIntervalStart, gammaIntervalPoint1, gammaIntervalPoint2, gammaIntervalEnd);
			printf("-- stablePoint1=%i stablePoint2=%i\n", stablePoint1, stablePoint2);

			if (!stablePoint1 && stablePoint2) {
				gammaIntervalStart = gammaIntervalPoint1;
				gammaIntervalPoint1 = gammaIntervalPoint2;
				gammaIntervalPoint2 = gammaIntervalStart + (gammaIntervalEnd-gammaIntervalStart)/goldenRatio;
				gammaFactor = gammaIntervalPoint2;
				fitInterfaceWidthToKappa = 1;
				setLBInitializationProfile();
				for (int i = 0; i < phase_iterations; i++) iteration();
				stablePoint1 = stablePoint2;
				if (!stableSimulation(n1) || !stableSimulation(n2)) stablePoint2 = 0;
			}
			else if (stablePoint1 && stablePoint2) {
				gammaIntervalEnd = gammaIntervalPoint1; // both known stable... skip more
				gammaIntervalPoint1 = gammaIntervalEnd - (gammaIntervalEnd-gammaIntervalStart)/goldenRatio;
				gammaIntervalPoint2 = gammaIntervalStart + (gammaIntervalEnd-gammaIntervalStart)/goldenRatio;

				gammaFactor = gammaIntervalPoint1;
				stablePoint1 = 1;
				fitInterfaceWidthToKappa = 1;
				setLBInitializationProfile();
				for (int i = 0; i < phase_iterations; i++) iteration();
				if (!stableSimulation(n1) || !stableSimulation(n2)) stablePoint1 = 0;

				gammaFactor = gammaIntervalPoint2;
				stablePoint2 = 1;
				fitInterfaceWidthToKappa = 1;
				setLBInitializationProfile();
				for (int i = 0; i < phase_iterations; i++) iteration();
				if (!stableSimulation(n1) || !stableSimulation(n2)) stablePoint2 = 0;
			}
			else { // ((stablePoint1 && !stablePoint2) || (!stablePoint1 && !stablePoint2)) ... restart search over back side of the gammaFactor range
				gammaIntervalStart = gammaIntervalPoint2; // known unstable
				gammaIntervalPoint1 = gammaIntervalEnd - (gammaIntervalEnd-gammaIntervalStart)/goldenRatio;
				gammaIntervalPoint2 = gammaIntervalStart + (gammaIntervalEnd-gammaIntervalStart)/goldenRatio;

				gammaFactor = gammaIntervalPoint1;
				stablePoint1 = 1;
				fitInterfaceWidthToKappa = 1;
				setLBInitializationProfile();
				for (int i = 0; i < phase_iterations; i++) iteration();
				if (!stableSimulation(n1) || !stableSimulation(n2)) stablePoint1 = 0;

				gammaFactor = gammaIntervalPoint2;
				stablePoint2 = 1;
				fitInterfaceWidthToKappa = 1;
				setLBInitializationProfile();
				for (int i = 0; i < phase_iterations; i++) iteration();
				if (!stableSimulation(n1) || !stableSimulation(n2)) stablePoint2 = 0;
			}
			j++; // golden search iteration counter
		}
		gammaFactor = gammaIntervalEnd * 1.03; // should be "guaranteed" stable with an extra 3% correction
		printf("determined gammaFactor=%f, for (%f,%f) with interfaceWidth=%f, kappa=%f, kappaFactor=%f, gammaMu=%f\n\n",
				gammaFactor, nA0, nB0, interfaceWidth, kappa, kappaFactor, gammaMu);
		return 1;
	}
} // end function fitGammaFactor()


/**
 * @brief Function _initializeInterfaceWidth_ sets the interface width used to initialize a lattice Boltzmann simulation.
 *
 * This function provides an interface width to be used in the tanh profile initialization of a LB simulation.  The interface width varies depending on
 * if kappa is automatically set by single-component theory (see reference) or manually set by a user:
 * 1. When kappa is automatically determined by single-component theory, the interface width is also automatically calculated.
 * 2. When kappa is manually set through the GUI, the interface width is set through a 3-step process:
 *    + The interface width is initially set according to single-component theory, and a simulation is iterated to measure the near-equilibrium width.
 *    + The equilibrium width is then used to modify the value of kappa (user-defined kappa / correction kappaFactor) to preserve the proportionality
 *      width = kappaFactor * sqrt(kappa).
 *    + The kappaFactor is saved for subsequent simulations of the (A,B) pair, and the interface width is calculated according to the proportionality.
 *    + When a simulation for a new (A,B) test pair is initialized, this loop begins again.
 *
 * @note This function sets a default minimum interface width of 2 lattice spaces.  The global flag _overrideMinimumInterfaceWidth_ overrides this value.
 *
 * @see A. Wagner and C. Pooley, Physical Review E 76, 045702 (2007)
 */
void initializeInterfaceWidth() {
	double rhoVapor = 0.0, rhoLiquid = 0.0;
	double tc = tcA < tcB ? tcA : tcB;
	double zeroThreshold = 1e-10;
	double tmpGamma = gammaMu;

	// Use the smaller of B1 and A2 as the theoretical smallest component density
	// ... assumes checkVaporPhase() has been run to order the phases properly
	if ( isEqual(nA0, 0.0, zeroThreshold) || (fabs(theoreticalDensityB2-theoreticalDensityB1) > fabs(theoreticalDensityA1-theoreticalDensityA2)) ) {
		rhoVapor = theoreticalDensityB1;
		rhoLiquid = theoreticalDensityB2;
	}
	else {
		rhoVapor = theoreticalDensityA2;
		rhoLiquid = theoreticalDensityA1;
	}

	if (determineKappaAutomatically) {
		interfaceWidth = 1.0 / sqrt(4.0*rhoVapor*fabs(tc-theta)); // added fabs for diagrams where a component is above critical temp
		kappa = 1.0 / (8.0 * theta * rhoVapor);
		kappa /= kappaFactor;
	}
	else { // from the relation width = kappaFactor * sqrt(kappa)
		static double staticKappaFactor = 1.0;

		switch (fitInterfaceWidthToKappa) {
		case 1: // start interface width fit: new nA0,nB0 needs interface width fit to specified kappa
			interfaceWidth = 1.0 / sqrt(4.0*rhoVapor*fabs(tc-theta)); // added fabs for diagrams where a component is above critical temp
			kappaFactor = 1.0;
			interfaceWidthForGivenKappa = XDIM / 20;
			break;
		case 2: // interface width fit complete: during interface width fit, parameters to initialize the LB profile
			interfaceWidth = interfaceWidthForGivenKappa;
			kappaFactor = interfaceWidth / sqrt(kappa);
			staticKappaFactor = kappaFactor;
			break;
		default: // subsequent runs that don't need a new interface width fit:  after RSS run and width achieved, need the kappa factor to determine the width
			kappaFactor = staticKappaFactor;
			interfaceWidth = kappaFactor * sqrt(kappa);
		}
	}
	overrideMinimumInterfaceWidth ? 1 : ((interfaceWidth >= 2.0) ? 1 : (interfaceWidth = 2.0)); // enforce minimum width to help stability at lower kappa values

	gammaMu = 1.0 / (6.0 * kappa * rhoLiquid);
	tmpGamma = gammaMu;
	gammaMu /= gammaFactor;
	printf("%i%i: interfaceWidth=%f, kappa=%f, kappaFactor=%f, gammaMu=%f, gammaFactor=%f (effective gammaMu=%f)\n",
			determineKappaAutomatically, fitInterfaceWidthToKappa, interfaceWidth, kappa, kappaFactor, tmpGamma, gammaFactor, gammaMu);
} // end function initializeInterfaceWidth()


/**
 * @brief Function _initializeRandom_ initializes a fully mixed 2-component simulation (i.e. equal densities across the lattice) with an additional noise term.
 *
 * @note The global value _Amp_ controls the amplitude of the noise term added to the density profile.
 * @note The global values _nA0_ and _nB0_ set the total mass/densities to be initialized.
 *
 * @param [in] i The lattice site to be initialized.
 */
void initializeRandom(int i) {
	if (nA0 < Amp || nB0 < Amp) Amp = 0.1 * (nA0 < nB0 ? nA0 : nB0); // make sure fluctuations don't give negative particle counts

	if (initializeRandomComponents == 1 || initializeRandomComponents == 2) Amp = 0.0; // no noise if only 1 component is random (only done if spoofing zero)
	if (initializeRandomComponents == 1 || initializeRandomComponents == 0) n1[i] = nA0 + Amp*(1.0*rand()/RAND_MAX-0.5);
	if (initializeRandomComponents == 2 || initializeRandomComponents == 0) n2[i] = nB0 + Amp*(1.0*rand()/RAND_MAX-0.5);

	if (initializeRandomComponents == 0) {
		phase1Index = 0;
		phase2Index = 0;
	}
} // end function initializeRandom()


/**
 * @brief Function _initializeTheoreticalTwoPhases_ initializes two density domains for a lattice Boltzmann simulation.
 *
 * This function initializes two density domains for a simulation, one A-rich domain and one B-rich domain, using tanh profiles as an approximation for the
 * domain shape.  The bulk density values are taken from the free energy minimization of the (A,B) pair that is to be simulated.  It's nominal use is for
 * simulations where the expected behavior is phase separation into 2 phases.
 *
 * @note The global values _nA0_ and _nB0_ set the total mass/densities to be initialized.
 *
 * @param [in] i The lattice site to be initialized.
 */
void initializeTheoreticalTwoPhases(int i) {
	int interface = 0.5 * XDIM;
	double transition = 0.0;
	double rhoA1 = theoreticalDensityA1, rhoB1 = theoreticalDensityB1;
	double rhoA2 = theoreticalDensityA2, rhoB2 = theoreticalDensityB2;

	// Chose the theoretical phases to use for the 2-phase initialization
	if (useTheoreticalPhase1 && useTheoreticalPhase2 && useTheoreticalPhase3) {
		useTheoreticalPhase3 = 0;
	}
	else if (useTheoreticalPhase1 && !useTheoreticalPhase2 && useTheoreticalPhase3) {
		rhoA2 = theoreticalDensityA3;
		rhoB2 = theoreticalDensityB3;
	}
	else if (!useTheoreticalPhase1 && useTheoreticalPhase2 && useTheoreticalPhase3) {
		rhoA1 = theoreticalDensityA3;
		rhoB1 = theoreticalDensityB3;
	}

	domain4 = XDIM;
	if (i < interface-0.25*XDIM) {
		transition = 0.5 + 0.5*tanh((double)(i%XDIM)/interfaceWidth);
		n1[i] = (1.0-transition)*rhoA2 + transition*rhoA1;
		n2[i] = (1.0-transition)*rhoB2 + transition*rhoB1;
	}
	else if (i >= interface-0.25*XDIM && i <= interface+0.25*XDIM) {
		transition = 0.5 + 0.5*tanh((double)(i-interface)/interfaceWidth);
		n1[i] = (1.0-transition)*rhoA1 + transition*rhoA2;
		n2[i] = (1.0-transition)*rhoB1 + transition*rhoB2;
	}
	else {
		transition = 0.5 + 0.5*tanh((double)((i-domain4)%XDIM)/interfaceWidth);
		n1[i] = (1.0-transition)*rhoA2 + transition*rhoA1;
		n2[i] = (1.0-transition)*rhoB2 + transition*rhoB1;
	}

	// If one component is practically zero, use the random initializer on only that component
	if (nA0 < 1e-10) {
		initializeRandomComponents = 1;
		initializeRandom(i);
		initializeRandomComponents = 0;
	}
	if (nB0 < 1e-10) {
		initializeRandomComponents = 2;
		initializeRandom(i);
		initializeRandomComponents = 0;
	}

	// Set indices from which to read equilibrium densities in 2-phase simulations
	phase1Index = interface-0.25*XDIM;
	phase2Index = interface+0.25*XDIM;
} // end function initializeTheoreticalTwoPhases()


/**
 * @brief Function _initializeTheoreticalThreePhases_ initializes density domains for three-phase behavior a lattice Boltzmann simulation.
 *
 * This function initializes four density domains for a simulation: one A-rich domain, one B-rich domain, and two vapor domains to separate the A-/B-rich phases.
 * Density profile shapes are approximated by tanh profiles, and the bulk density values for all three phases are taken from the free energy minimization of the
 * (A,B) pair that is to be simulated.  It's nominal use is for simulations where the expected behavior is phase separation into 3 phases.
 *
 * @note The global values _nA0_ and _nB0_ set the total mass/densities to be initialized.
 * @note The global flag _useTheoreticalVolumina_ switches between initializations with equal domain widths or widths according to the free energy minimization.
 * @note The global flag _initializeStepProfile_ replaces the tanh profile approximation with a step profile.
 *
 * @param [in] i The lattice site to be initialized.
 */
void initializeTheoreticalThreePhases(int i) {
	int half = 0.125 * XDIM;
	double rhoA1 = theoreticalDensityA1, rhoA2 = theoreticalDensityA2, rhoA3 = theoreticalDensityA3;
	double rhoB1 = theoreticalDensityB1, rhoB2 = theoreticalDensityB2, rhoB3 = theoreticalDensityB3;

	if (useTheoreticalVolumina) {
		domain1Width = (int)(theoreticalVolume1 * XDIM);
		domain2Width = (int)(0.5 * theoreticalVolume3 * XDIM);
		domain3Width = (int)(theoreticalVolume2 * XDIM);
		domain4Width = domain2Width;
		domain2 = XDIM / 2; // fix center of domain2 (vapor) to the center of the lattice
		domain1 = domain2 - domain2Width/2 - domain1Width/2;
		domain3 = domain2 + domain2Width/2 + domain3Width/2;
		domain4 = (domain3 + domain3Width/2 + domain4Width/2) % XDIM;
	}
	else {
		domain1 = XDIM / 4;
		domain2 = XDIM / 2;
		domain3 = 3 * XDIM / 4;
		domain4 = XDIM;
		domain1Width = domain1;
		domain2Width = domain1Width;
		domain3Width = domain1Width;
		domain4Width = domain2Width;
	}

	if (suppressTheoreticalPhase1 && !suppressTheoreticalPhase2) {
		rhoA1 = rhoA3;
		rhoB1 = rhoB3;
	}
	else if (!suppressTheoreticalPhase1 && suppressTheoreticalPhase2) {
		rhoA2 = rhoA3;
		rhoB2 = rhoB3;
	}
	else {
		suppressTheoreticalPhase1 = 0;
		suppressTheoreticalPhase2 = 0;
	}

	double transition = 0.0;
	if (i >= (domain4 < domain1 ? domain4 : 0) && i < domain1) { // vapor to A-liquid transition
		transition = 0.5 + 0.5*tanh((double)((i-(domain1-domain1Width/2))%XDIM)/interfaceWidth);
		n1[i] = (1.0-transition)*rhoA3 + transition*rhoA1;
		n2[i] = (1.0-transition)*rhoB3 + transition*rhoB1;
	}
	else if (i >= domain1 && i < domain2) { // A-liquid to vapor transition
		transition = 0.5 + 0.5*tanh((double)(i-(domain2-domain2Width/2))/interfaceWidth);
		n1[i] = (1.0-transition)*rhoA1 + transition*rhoA3;
		n2[i] = (1.0-transition)*rhoB1 + transition*rhoB3;
	}
	else if (i >= domain2 && i < domain3) { // vapor to B-liquid transition
		transition = 0.5 + 0.5*tanh((double)(i-(domain3-domain3Width/2))/interfaceWidth);
		n1[i] = (1.0-transition)*rhoA3 + transition*rhoA2;
		n2[i] = (1.0-transition)*rhoB3 + transition*rhoB2;
	}
	else if (i >= domain3 && i < (domain4 > domain3 ? domain4 : domain4 + XDIM)) { // B-liquid to vapor transition
		transition = 0.5 + 0.5*tanh((double)(i-(domain4-domain4Width/2+(domain4 > domain3 ? 0 : XDIM)))/interfaceWidth);
		n1[i] = (1.0-transition)*rhoA2 + transition*rhoA3;
		n2[i] = (1.0-transition)*rhoB2 + transition*rhoB3;
	}
	else { // since domain 4 can be on either side of XDIM, connect it to domain 1 accordingly
		if (domain4 > domain3 && i >= domain4) { // domain4 is "left" of XDIM
			transition = 0.5 + 0.5*tanh((double)((i-(domain1-domain1Width/2+XDIM))%XDIM)/interfaceWidth);
			n1[i] = (1.0-transition)*rhoA3 + transition*rhoA1;
			n2[i] = (1.0-transition)*rhoB3 + transition*rhoB1;
		}
		else if (domain4 < domain1 && i < domain4) { // domain4 is "right" of XDIM
			transition = 0.5 + 0.5*tanh((double)(i-(domain4-domain4Width/2))/interfaceWidth);
			n1[i] = (1.0-transition)*rhoA2 + transition*rhoA3;
			n2[i] = (1.0-transition)*rhoB2 + transition*rhoB3;
		}
	}

	// If one component is practically zero, use the random initializer on only that component
	if (nA0 < 1e-10) {
		initializeRandomComponents = 1;
		initializeRandom(i);
		initializeRandomComponents = 0;
	}
	if (nB0 < 1e-10) {
		initializeRandomComponents = 2;
		initializeRandom(i);
		initializeRandomComponents = 0;
	}

	// Set indicies from which to read equilibrium densities
	phase1Index = domain1-half;
	phase2Index = domain3-half;
} // end function initializeTheoreticalThreePhases()


/**
 * @brief Function _initializeTheoreticalFourPhases_ initializes density domains for a hard-coded four-phase density profile to start a lattice Boltzmann simulation.
 *
 * This function is an attempt to observe 4-phase behavior in a 2-component mixture.  It initializes four density domains for a simulation: one A-rich domain, a 4th intermediate
 * density phase, one B-rich domain, and a vapor domain.  Density profile shapes are approximated by tanh profiles, and the bulk density values for the liquid and vapor phases are
 * taken from the free energy minimization of the mixture defined by tcA = tcB = 0.475, VDW interaction parameter = 0.5636.  The default arrangement of the phase domains across the
 * lattice is vapor, A-liquid, 4th phase, B-liquid.  Toggling the flag _isolateFourthPhase_ changes the domain arrangement to A-liquid, B-liquid, vapor, 4th phase, vapor.
 *
 * @note The densities used to define each phase are not associated with a specific starting (A,B) mixture in the phase diagram.  The triangle points that define the three-phase
 * region for the specified phase diagram were augmented by approximating the density of the 4th phase by phase diagram inspection.
 *
 * @param [in] i The lattice site to be initialized.
 */
void initializeTheoreticalFourPhases(int i) {
	double rhoA1 = 1.896553132497407, rhoA2 = 0.169150765439082, rhoA3 = 0.172977697951899, rhoA4 = 0.71;
	double rhoB1 = 0.169151149061583, rhoB2 = 1.896556831396832, rhoB3 = 0.172978477952882, rhoB4 = 0.71;
	double transition = 0.0;

	nA0 = 4.0;
	nB0 = 4.0;
	theoreticalDensityA1 = 0.0;
	theoreticalDensityB1 = 0.0;
	theoreticalDensityA2 = 0.0;
	theoreticalDensityB2 = 0.0;
	theoreticalDensityA3 = 0.0;
	theoreticalDensityB3 = 0.0;
	theoreticalVolume1 = 0.0;
	theoreticalVolume2 = 0.0;
	theoreticalVolume3 = 0.0;
	tcA = 0.475;
	tcB = 0.475;
	pcB = 1.0;

	interfaceWidth = 2.0;
	kappa = 0.2;
	gammaMu = 0.1;
	vdwInteractionFactor = 0.5636;

	if (isolateFourthPhase) { // alter the default domain arrangement below: A-liquid, B-liquid, vapor, 4th phase, vapor
		domain1 = XDIM / 10;
		domain2 = domain1 + XDIM/5;
		domain3 = domain2 + 3*XDIM/10;
		domain4 = domain3 + XDIM/5;
		int domain5 = XDIM;
		domain1Width = domain1 / 2;
		domain2Width = 2 * domain1;
		domain3Width = 2 * domain1;
		domain4Width = 2 * domain1;
		int domain5Width = domain1Width;
		if (i >= 0 && i < domain1) {
			transition = 0.5 + 0.5*tanh((double)((i-(domain1-domain1Width/2))%XDIM)/interfaceWidth);
			n1[i] = (1.0-transition)*rhoA2 + transition*rhoA1;
			n2[i] = (1.0-transition)*rhoB2 + transition*rhoB1;
		}
		else if (i >= domain1 && i < domain2) {
			transition = 0.5 + 0.5*tanh((double)(i-(domain2-domain2Width/2))/interfaceWidth);
			n1[i] = (1.0-transition)*rhoA1 + transition*rhoA3;
			n2[i] = (1.0-transition)*rhoB1 + transition*rhoB3;
		}
		else if (i >= domain2 && i < domain3) {
			transition = 0.5 + 0.5*tanh((double)(i-(domain3-domain3Width/2))/interfaceWidth);
			n1[i] = (1.0-transition)*rhoA3 + transition*rhoA4;
			n2[i] = (1.0-transition)*rhoB3 + transition*rhoB4;
		}
		else if (i >= domain3 && i < domain4) {
			transition = 0.5 + 0.5*tanh((double)(i-(domain4-domain4Width/2))/interfaceWidth);
			n1[i] = (1.0-transition)*rhoA4 + transition*rhoA3;
			n2[i] = (1.0-transition)*rhoB4 + transition*rhoB3;
		}
		else if (i >= domain4 && i < domain5) {
			transition = 0.5 + 0.5*tanh((double)(i-(domain5-domain5Width/2))/interfaceWidth);
			n1[i] = (1.0-transition)*rhoA3 + transition*rhoA2;
			n2[i] = (1.0-transition)*rhoB3 + transition*rhoB2;
		}
	}
	else { // default arrangement of domains: vapor, A-liquid, 4th phase, B-liquid
		domain1 = XDIM / 4;
		domain2 = XDIM / 2;
		domain3 = 3 * XDIM / 4;
		domain4 = XDIM;
		domain1Width = domain1;
		domain2Width = domain1Width;
		domain3Width = domain1Width;
		domain4Width = domain2Width;
		if (i >= (domain4 < domain1 ? domain4 : 0) && i < domain1) { // vapor to A-liquid transition
			transition = 0.5 + 0.5*tanh((double)((i-(domain1-domain1Width/2))%XDIM)/interfaceWidth);
			n1[i] = (1.0-transition)*rhoA3 + transition*rhoA1;
			n2[i] = (1.0-transition)*rhoB3 + transition*rhoB1;
		}
		else if (i >= domain1 && i < domain2) { // A-liquid to 4th phase transition
			transition = 0.5 + 0.5*tanh((double)(i-(domain2-domain2Width/2))/interfaceWidth);
			n1[i] = (1.0-transition)*rhoA1 + transition*rhoA4;
			n2[i] = (1.0-transition)*rhoB1 + transition*rhoB4;
		}
		else if (i >= domain2 && i < domain3) { // 4th phase to B-liquid transition
			transition = 0.5 + 0.5*tanh((double)(i-(domain3-domain3Width/2))/interfaceWidth);
			n1[i] = (1.0-transition)*rhoA4 + transition*rhoA2;
			n2[i] = (1.0-transition)*rhoB4 + transition*rhoB2;
		}
		else if (i >= domain3 && i < (domain4 > domain3 ? domain4 : domain4 + XDIM)) { // B-liquid to vapor transition
			transition = 0.5 + 0.5*tanh((double)(i-(domain4-domain4Width/2+(domain4 > domain3 ? 0 : XDIM)))/interfaceWidth);
			n1[i] = (1.0-transition)*rhoA2 + transition*rhoA3;
			n2[i] = (1.0-transition)*rhoB2 + transition*rhoB3;
		}
		else { // since domain 4 can be on either side of XDIM, connect it to domain 1 accordingly
			if (domain4 > domain3 && i >= domain4) { // domain4 is "left" of XDIM
				transition = 0.5 + 0.5*tanh((double)((i-(domain1-domain1Width/2+XDIM))%XDIM)/interfaceWidth);
				n1[i] = (1.0-transition)*rhoA3 + transition*rhoA1;
				n2[i] = (1.0-transition)*rhoB3 + transition*rhoB1;
			}
			else if (domain4 < domain1 && i < domain4) { // domain4 is "right" of XDIM
				transition = 0.5 + 0.5*tanh((double)(i-(domain4-domain4Width/2))/interfaceWidth);
				n1[i] = (1.0-transition)*rhoA2 + transition*rhoA3;
				n2[i] = (1.0-transition)*rhoB2 + transition*rhoB3;
			}
		}
	}
} // end function initializeTheoreticalFourPhases()


/**
 * @brief Function _setInitializeRandom_ sets the simulation to initialize a uniform mixture.
 */
void setInitializeRandom() {
	initializeProfile = initializeRandom;
	initialize();
} // end function setInitializeRandom()


/**
 * @brief Function _setInitializeTheoreticalTwoPhases_ sets the simulation to initialize 2-phase behavior.
 */
void setInitializeTheoreticalTwoPhases() {
	initializeProfile = initializeTheoreticalTwoPhases;
	initialize();
} // end function setInitializeTheoreticalTwoPhases()


/**
 * @brief Function _setInitializeTheoreticalThreePhases_ sets the simulation to initialize 3-phase behavior.
 */
void setInitializeTheoreticalThreePhases() {
	initializeProfile = initializeTheoreticalThreePhases;
	initialize();
} // end function setInitializeTheoreticalThreePhases()


/**
 * @brief Function _setInitializeTheoreticalFourPhases_ sets the simulation to initialize a hard-coded 4-phase density profile.
 */
void setInitializeTheoreticalFourPhases() {
	initializeProfile = initializeTheoreticalFourPhases;
	useTheoreticalDensities = 0;
	initialize();
} // end function setInitializeTheoreticalFourPhases()


/**
 * @brief Function _setLBInitializationProfile_ is a the controlling function to choose the simulation initialization method.
 *
 * @note The global value _useTwoPhasesStepInitialization_ controls a switch statement to change the initialization method (1=random, 2=2-phase, 3=3-phase).
 */
void setLBInitializationProfile() {
	if (useStepOrRandomInitialization) {
		switch (useTwoOrThreePhaseInitialization) {
		case 2:
			setInitializeTheoreticalTwoPhases();
			break;
		case 3:
			setInitializeTheoreticalThreePhases();
			break;
		default:
			printf("Defaulting to 2-phase tanh profile...\n\n");
			useTwoOrThreePhaseInitialization = 2;
			setInitializeTheoreticalTwoPhases();
		}
	}
	else setInitializeRandom();
} // end function setTwoPhaseInitialization()


/**
 * @brief Function _getTheoreticalDensities_ retrieves the free energy minimization results for a given (A,B) test pair.
 *
 * @note This function populates all theoretical density and volumina variables with the value 0.0 if no minimization results are found.
 *
 * @return Integer value indicating theoretical values were successfully (1) or unsuccessfully (0) found.
 */
int getTheoreticalDensities() {
	int readEOF;
	int match = 0;
	double tmpParticlesATotal = 0.0, tmpParticlesBTotal = 0.0;
	double tmpRhoA1 = 0.0, tmpRhoB1 = 0.0;
	double tmpRhoA2 = 0.0, tmpRhoB2 = 0.0;
	double tmpRhoA3 = 0.0, tmpRhoB3 = 0.0;
	double tmpVolume1 = 0.0, tmpVolume2 = 0.0, tmpVolume3 = 0.0;
	double matchThreshold = 1e-9;

	//
	// Check 2-phase theoretical densities
	//
	FILE *phaseDiagram_densities_twoPhases = openFile("threePhaseDiagram-densities-twoPhases", ".dat", "r");
	if (phaseDiagram_densities_twoPhases) {
		readEOF = fscanf(phaseDiagram_densities_twoPhases, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf",
				&tmpParticlesATotal, &tmpParticlesBTotal, &tmpRhoA1, &tmpRhoB1, &tmpRhoA2, &tmpRhoB2, &tmpRhoA3, &tmpRhoB3,
				&tmpVolume1, &tmpVolume2, &tmpVolume3);
		while (readEOF != EOF) {
			if (isEqual(tmpParticlesATotal, nA0, matchThreshold) &&	isEqual(tmpParticlesBTotal, nB0, matchThreshold)) {
				match = 1;
				theoreticalPhases = 2;
				break;
			}
			readEOF = fscanf(phaseDiagram_densities_twoPhases, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf",
					&tmpParticlesATotal, &tmpParticlesBTotal, &tmpRhoA1, &tmpRhoB1, &tmpRhoA2, &tmpRhoB2, &tmpRhoA3, &tmpRhoB3,
					&tmpVolume1, &tmpVolume2, &tmpVolume3);
		}
		fclose(phaseDiagram_densities_twoPhases);
	}

	//
	// Check 3-phase theoretical densities
	//
	if (!match) { // no 2-phase matches found
		FILE *phaseDiagram_densities_threePhases = openFile("threePhaseDiagram-densities-threePhases", ".dat", "r");
		if (phaseDiagram_densities_threePhases) {
			readEOF = fscanf(phaseDiagram_densities_threePhases, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf",
					&tmpParticlesATotal, &tmpParticlesBTotal, &tmpRhoA1, &tmpRhoB1, &tmpRhoA2, &tmpRhoB2, &tmpRhoA3, &tmpRhoB3,
					&tmpVolume1, &tmpVolume2, &tmpVolume3);
			while (readEOF != EOF) {
				if (isEqual(tmpParticlesATotal, nA0, matchThreshold) &&	isEqual(tmpParticlesBTotal, nB0, matchThreshold)) {
					match = 1;
					theoreticalPhases = 3;
					break;
				}
				readEOF = fscanf(phaseDiagram_densities_threePhases, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf",
						&tmpParticlesATotal, &tmpParticlesBTotal, &tmpRhoA1, &tmpRhoB1, &tmpRhoA2, &tmpRhoB2, &tmpRhoA3, &tmpRhoB3,
						&tmpVolume1, &tmpVolume2, &tmpVolume3);
			}
			fclose(phaseDiagram_densities_threePhases);
		}
	}

	if (match) {
		theoreticalDensityA1 = tmpRhoA1;
		theoreticalDensityB1 = tmpRhoB1;
		theoreticalDensityA2 = tmpRhoA2;
		theoreticalDensityB2 = tmpRhoB2;
		theoreticalDensityA3 = tmpRhoA3;
		theoreticalDensityB3 = tmpRhoB3;
		theoreticalVolume1 = tmpVolume1;
		theoreticalVolume2 = tmpVolume2;
		theoreticalVolume3 = tmpVolume3;
	}
	else { // zero out theoretical densities if there are none found
		theoreticalDensityA1 = 0.0;
		theoreticalDensityB1 = 0.0;
		theoreticalDensityA2 = 0.0;
		theoreticalDensityB2 = 0.0;
		theoreticalDensityA3 = 0.0;
		theoreticalDensityB3 = 0.0;
		theoreticalVolume1 = 0.0;
		theoreticalVolume2 = 0.0;
		theoreticalVolume3 = 0.0;
		theoreticalPhases = 0;
		printf("\nNo matching theoretical densities found for %f, %f... the step profile initialization may produce an inaccurate LB simulation.\n\n", nA0, nB0);
	}

	return match;
} // end function getTheoreticalDensities()


/**
 * @brief Function _initialize_ is the controlling function to initialize a lattice Boltzmann simulation.
 *
 * This function is the controlling entity for all lattice Boltzmann simulations.  When initializing 2-phase or 3-phase behavior, the default is to use bulk
 * densities from a free energy minimization along with a manually set value for kappa and an interface width that is automatically fit to that kappa.
 *
 * @note The global flag _useTheoreticalDensities_ controls if the minimization results for an (A,B) test pair are to be used to initialize bulk phases.
 * @note The global flag _determineInterfaceWidthAutomatically_ controls if the initialized interface width is automatically or manually set.  If a manual
 * interface width is chosen, the initialization forces kappa to use a manual value as well.
 * @note The global flag _determineKappaAutomatically_ controls if the initialized kappa value is automatically or manually set.
 * @note Once the interface width is automatically fit to a specific manual kappa value, subsequent initializations use the same value of the global parameter
 * _kappaFactor_ to change the interface with when kappa is changed (width = kappaFactor * sqrt(kappa)).  If this results unstable simulations for kappas other
 * than the one used for the interface fit, setting the global flag _kappaFactorDetermined_ to "off" will force a new interface fit for a new kappa value.
 * @note The value for gammaMu is always automatically determined.  However, it may be manually set using a new value for the global parameter _gammaFactor_.
 */
void initialize() {
	static double last_tcA = 0.0, last_tcB = 0.0, last_ncB = 0.0, last_vdw = 0.0;
	double threshold = 1e-6;

	// Detect if this is a new initialization or a re-initialization of the same parameters
	if (!isEqual(last_tcA, tcA, threshold) || !isEqual(last_tcB, tcB, threshold) ||
			!isEqual(last_ncB, ncB, threshold) || !isEqual(last_vdw, vdwInteractionFactor, threshold)) {
		last_tcA = tcA;
		last_tcB = tcB;
		last_ncB = ncB;
		last_vdw = vdwInteractionFactor;
		printf("\n-------------------------------------------------------------------------------------------------------------------------------------------\n");
		printf("- Initializing: tcA=%f tcB=%f ncB=%f vdwInteractionFactor=%f\n", tcA, tcB, ncB, vdwInteractionFactor);
		printf("---------------------------------------------------------------------------------------------------------------------------------------------\n");

		// New initializations need a new 3-phase region defined
		threePhaseRegion[0] = 0;
		threePhaseRegion[1] = 0;
		threePhaseRegion[2] = 0;
		threePhaseRegion[3] = 0;
		threePhaseRegion[4] = 0;
		threePhaseRegion[5] = 0;
		maxArea = 0.0;
		setThreePhaseRegion = 1;
		threePhaseRegionExists = 0;
	}

	// tcA, tcB, and ncB vary as the degrees of freedom
	pcA = 3.*tcA/8.; // determine pcA
	ncA = pcA / ((3./8.)*tcA); // fix ncA to be 1
	pcB = (3./8.)*tcB*ncB; // determine pcB

	// Reset the values of the VDW constants for each component
	aA = (27./64.)*(tcA*tcA/pcA);
	aB = (27./64.)*(tcB*tcB/pcB);
	aAB = sqrt(aA*aB) * vdwInteractionFactor; //0.5*(aA+aB); other random mixing rule
	bA = tcA/(8.*pcA);
	bB = tcB/(8.*pcB);

	// Minimize the initial A,B values to determine theoretical values for pressure and chemical potential
	if (useTheoreticalDensities) {
		if (getTheoreticalDensities()) {

			// Define a density pair for the A-rich and B-rich phases (not using vapor phase)
			double rhoA1 = theoreticalDensityA1;
			double rhoB1 = theoreticalDensityB1;
			double rhoA2 = theoreticalDensityA2;
			double rhoB2 = theoreticalDensityB2;

			// Calculate the bulk values only for the theoretical expectations for pressure and chemical potential
			theoreticalPressure = (rhoA1+rhoB1)*theta*(1 + (bA*rhoA1+bB*rhoB1)/(1.-bA*rhoA1-bB*rhoB1)) - aA*rhoA1*rhoA1 - 2*aAB*rhoA1*rhoB1 - aB*rhoB1*rhoB1; // just use phase 1 here
			if (rhoA1 != 0.0) theoreticalMuA = theta*log(rhoA1/(1.-bA*rhoA1-bB*rhoB1)) + theta*bA*(rhoA1+rhoB1)/(1.-bA*rhoA1-bB*rhoB1) - 2*aA*rhoA1 - 2*aAB*rhoB1; // A-rich for muA
			else theoreticalMuA = 0.0;
			if (rhoB2 != 0.0) theoreticalMuB = theta*log(rhoB2/(1.-bA*rhoA2-bB*rhoB2)) + theta*bB*(rhoA2+rhoB2)/(1.-bA*rhoA2-bB*rhoB2) - 2*aB*rhoB2 - 2*aAB*rhoA2; // B-rich for muB
			else theoreticalMuB = 0.0;

			if (determineInterfaceWidthAutomatically && initializeProfile != initializeRandom) {
				static double lastNA0 = 0.0, lastNB0 = 0.0;

				if (!determineKappaAutomatically && (!isEqual(nA0, lastNA0, threshold) || !isEqual(nB0, lastNB0, threshold) || !kappaFactorDetermined)) {
					lastNA0 = nA0;
					lastNB0 = nB0;
					gammaFactor = 1.0;
					fitInterfaceWidthToKappa = 1; // need to get a kappaCoefficient for the new point (switch in initializeInterfaceWidth())
					goodInterfaceFit = 0; // assume the fit will work
					fitInterfaceWidth();
					if (interfaceWidthForGivenKappa > XDIM) { // invalid interface width fit set to 10.0*XDIM in function fitInterfaceWidth()
						if (fitGammaFactor()) { // if a new gamma factor is found to stabilize LB for an interface fit...
							fitInterfaceWidth();
							if (interfaceWidthForGivenKappa < XDIM) goodInterfaceFit = 1;
						}
						else {
							printf("\nFailed to find a stable LB simulation for (%f,%f) at kappa=%f, gammaMu(p0) from %f to %f.\n",
									nA0, nB0, kappa, gammaMu, gammaMu/gammaFactor);
							printf("Falling back to the automatic single-component theory values.\n\n");
							kappaFactor = 1.0;
							gammaFactor = 1.0;
							determineKappaAutomatically = 1;
							kappaFactorDetermined = 0;
						}
					}
					else goodInterfaceFit = 1;
				}
				initializeInterfaceWidth();
				if (goodInterfaceFit) fitInterfaceWidthToKappa = 0; // new kappaCoefficient; can be used from here on (switch in initializeInterfaceWidth())
				else fitInterfaceWidthToKappa = 1; // if can't fit interface but still want a manual kappa...
			}
			else {
				determineKappaAutomatically = 0;
				fitInterfaceWidthToKappa = 9; // manual interface width
				printf("%i%i: interfaceWidth=%f, kappa=%f, gammaMu=%f\n", determineKappaAutomatically, fitInterfaceWidthToKappa, interfaceWidth, kappa, gammaMu);
			}
		}
		else { // no theoretical densities were found
			printf("Initializing a random profile since there are no theoretical densities\n\n");
			initializeProfile = initializeRandom;
			theoreticalPressure = 0.0;
			theoreticalMuA = 0.0;
			theoreticalMuB = 0.0;
		}
	} // end if (useTheoreticalDensities)

	iterations = 0;
	tmp_phase_iterations = phase_iterations;
	setLineAPoint = 1;
	setLineBPoint = 1;

	// Loop to initialize lattice values
	nAIntegrated = 0.0;
	nBIntegrated = 0.0;
	for (int i = 0; i < XDIM; i++) {
		initializeProfile(i);
		n[i] = n1[i] + n2[i]; // div by XDIM to scale the total particle count appropriately per lattice site
		nReduced[i] = n[i];
		nAIntegrated += n1[i];
		nBIntegrated += n2[i];
		u1[i] = 0.0;
		u2[i] = 0.0;
		uHat1[i] = 0.0;
		uHat2[i] = 0.0;
		u[i] = 0.0;

		// Initialize 1st component
		f1_0[i] = n1[i] - n1[i]*theta - n1[i]*u1[i]*u1[i];            		// zero velocity
		f1_1[i] = 0.5 * (n1[i]*u1[i]*u1[i] + n1[i]*u1[i] + n1[i]*theta); 	// +1 velocity (moving right)
		f1_2[i] = 0.5 * (n1[i]*u1[i]*u1[i] - n1[i]*u1[i] + n1[i]*theta); 	// -1 velocity (moving left)
		psi1[i] = 0.0;
		mu1[i] = 0.0;
		muNonIdeal1[i] = 0.0;

		// Initialize 2nd component
		f2_0[i] = n2[i] - n2[i]*theta - n2[i]*u2[i]*u2[i];            		// zero velocity
		f2_1[i] = 0.5 * (n2[i]*u2[i]*u2[i] + n2[i]*u2[i] + n2[i]*theta); 	// +1 velocity (moving right)
		f2_2[i] = 0.5 * (n2[i]*u2[i]*u2[i] - n2[i]*u2[i] + n2[i]*theta); 	// -1 velocity (moving left)
		psi2[i] = 0.0;
		mu2[i] = 0.0;
		muNonIdeal2[i] = 0.0;

		// Also scale the theoretical expectations by the appropriate gamma factor
		theoreticalMuAArray[i] = theoreticalMuA * gammaMu;
		theoreticalMuBArray[i] = theoreticalMuB * gammaMu;
		theoreticalPressureArray[i] = theoreticalPressure * gammaP;

		pressure[i] = 0.0;
		gradMuForce1[i] = 0.0;
		gradMuForce2[i] = 0.0;
		friction1[i] = 0.0;
		friction2[i] = 0.0;
		F1[i] = 0.0;
		F2[i] = 0.0;
	} // end for

	nAIntegrated /= XDIM;
	nBIntegrated /= XDIM;
} // end function initialize()


/**
 * @brief Function _init_ is a wrapper function for the very first program/simulation initialization.  It sets the default forcing method to the "log"
 * implementation of chemical potential forcing and the random profile for the density initialization.
 *
 * @note This function is only called once by the main() function on program start-up.
 */
void init() {
	setCollisionForcingNewChemicalPotentialGradient();
	initializeProfile = initializeRandom;
	initialize();
} // end function init()

