/*
 * LB_Analyze.c
 *
 *  Created on: Jul 7, 2016
 *      Author: clark
 */


#include "LB_D1Q3_1-component.h"


int phase_iterations = 250000;


double gradient(double *var, int i) {

	double result = 0;
	int ipp = i + 2;
	int ip = i + 1;
	int im = i - 1;
	int imm = i - 2;

	// Gradient for a lattice with periodic BCs
	// Instead of a halo, re-set the out-of-bounds index to the opposite end of the lattice
	if (useBoundaryConditionsPeriodic) {
		if (ip == wall || ip == XDIM) {
			ip = 0;
		}
		else if (im == -1 || im == wall-1) {
			im = XDIM - 1;
		}
		result = 1./2. * (var[ip]-var[im]);
	}
	// Gradient for a lattice with walls at the boundaries
	// Adjust for existence of the walls with forward/backward differencing
	else {
		if (ip == wall || ip == XDIM) {			// backward difference
			//result = var[i] - var[im];
			result = 1./2. * (var[imm] - 4.*var[im] + 3.*var[i]);
		}
		else if (im == -1 || im == wall-1) {	// forward difference
			//result = var[ip] - var[i];
			result = 1./2. * (-3.*var[i] + 4.*var[ip] - var[ipp]);
		}
		else {									// central difference... normal case
			result = 1./2. * (var[ip]-var[im]);
		}
	}

	return result;

} // end function gradient()


double laplace(double *var, int i) {

	double result = 0;
	int ippp = i + 3;
	int ipp = i + 2;
	int ip = i + 1;
	int im = i - 1;
	int imm = i - 2;
	int immm = i - 3;

	// Gradient for a lattice with periodic BCs
	// Instead of a halo, re-set the out-of-bounds index to the opposite end of the lattice
	if (useBoundaryConditionsPeriodic) {
		if (ip == wall || ip == XDIM) {
			ip = 0;
		}
		else if (im == -1 || im == wall-1) {
			im = XDIM - 1;
		}
		result = var[ip] - 2.*var[i] + var[im];
	}
	// Gradient for a lattice with walls at the boundaries
	// Adjust for existence of the walls with forward/backward differencing
	else {
		if (ip == wall || ip == XDIM) {					// backward difference
			//result = var[i] - 2.*var[im] + var[imm];
			result = -1.*var[immm] + 4.*var[imm] -5.*var[im] + 2.*var[i];
		}
		else if (im == -1 || im == wall-1) {			// forward difference
			//result = var[ipp] - 2.*var[ip] + var[i];
			result = 2.*var[i] - 5.*var[ip] + 4.*var[ipp] - var[ippp];
		}
		else {											// central difference... normal case
			result = var[ip] - 2.*var[i] + var[im];
		}
	}
	// TODO: neutralize BCs if needed so no wetting
 	return result;

} // end function laplace()


/**
 * This function returns the index of the array location containing the minimum value of that array.
 */
int findIndexOfMinimumValue(double *array, int arrayLength) {

    int i = 0;
    int indexOfMinimum = 0;

    for (i = 1; i < arrayLength; i++) {
        if (array[i] < array[indexOfMinimum]) {
            indexOfMinimum = i;
        }
    }

    return indexOfMinimum;

} // end function findIndexOfMinimumValue()


/**
 * This function returns the index of the array location containing the minimum value of that array.
 * Only handles square 2-D arrays.
 */
void findIndexOfMinimumValue2DArray(double **array, int arraySquareSize, int *minimumRow, int *minimumColumn) {

    int i = 0;
    int j = 0;
    int currentMinimumRow = 0;
    int currentMinimumColumn = 0;

    for (i = 0; i < arraySquareSize; i++) {
    	for (j = 0; j < arraySquareSize; j++) {
    		if (array[i][j] < array[currentMinimumRow][currentMinimumColumn]) {
    			currentMinimumRow = i;
    			currentMinimumColumn = j;
    		}
    	}
    }

    *minimumRow = currentMinimumRow;
    *minimumColumn = currentMinimumColumn;

} // end function findIndexOfMinimumValue2DArray()


/**
 * This function calculates the difference between gradP and rhoGradMu forcing.
 * Data is read directly from global arrays.
 */
void compareGradPRhoGradMu() {

	int i = 0;

	for (i = 0; i < XDIM; i++) {
		gradP[i] = gradient(pressure, i);
		rhoGradMu[i] = n1[i]*gradient(mu1,i);
	}

	for (i = 0; i < XDIM; i++) {
		gradPMinusRhoGradMu[i] = gradP[i] - rhoGradMu[i];
	}

} // end function compareGradPRhoGradMu()


/**
 * Calculate the points needed to graphically display a free energy curve for a given temperature.
 * Also calculates the points to graphically display the phase-separated densities that minimize the overall free energy.
 * This function only considers single-component phases without gradient terms in the free energy.
 * The resolution of the density step size limits the application of this function to the creation of coarse visuals only.
 */
//void calculateFreeEnergySingleComponentMinimization() {
//
//	int i = 0;
//	int rhoTail = 0;
//	int minimumCurvatureIndex = 0;
//	int trackMinimumCurvature = 1;
//	double volumeTotal = 1.;
//	double volumeLiquid = 0.5;
//	double volumeVapor = volumeTotal - volumeLiquid;  // sets constant volume fractions of .5
//	double rhoTotal = 0;
//	double rhoMin = 0.01;
//	double rhoMax = 5.01;
//	double rhoStepSize = 0.001;
//	double rho = 0;
//	double minRho1 = 0;
//	double minFreeEnergy1 = 0;
//	double minRho2 = 0;
//	double minFreeEnergy2 = 0;
//	double minimumCurvature = 0;  // start this at zero because I'll need a negative curvature to find the right spot to start
//	double freeEnergy = 0;
//	double freeEnergyTmp = 0;
//	double freeEnergyTotal = 0;
//	double volumeExclusion = 0;
//
//	freeEnergyArraySize = (rhoMax - rhoMin) / rhoStepSize;
//	freeEnergyArray = malloc(freeEnergyArraySize * sizeof(double));
//
//	double rhoArray[freeEnergyArraySize];
//	double freeEnergySlopeArray[freeEnergyArraySize];
//	double freeEnergyCurvatureArray[freeEnergyArraySize];
//
//	FILE *free_energy_results;
//	FILE *free_energy_minimum;
//	FILE *free_energy_slopes;
//	FILE *free_energy_curvature;
//	free_energy_results = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/free-energy-curve.dat","w");
//	free_energy_minimum = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/free-energy-minimum-points.dat","w");
//	free_energy_slopes = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/free-energy-slopes.dat","w");
//	free_energy_curvature = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/free-energy-curvature.dat","w");
//
//	// Develop the free energy curve that is to be minimized
//	for (i = 0; i < freeEnergyArraySize; i++) {
//		rho += rhoStepSize;
//		volumeExclusion = 1 - b1*rho;
//		if (volumeExclusion > 0) {
//			freeEnergy = rho*theta*log(rho/volumeExclusion) - a1*rho*rho; // - 0.5*kappa*gradient(rhoArray,i)*gradient(rhoArray,i);
//			if (freeEnergy < 0) {  // cut the curve off once it crosses into positive territory
//				rhoArray[i] = rho;
//				freeEnergyArray[i] = freeEnergy;
//				fprintf(free_energy_results, "%f %f\n", rhoArray[i], freeEnergyArray[i]);
//				rhoTail = i;
//			}
//			else {
//				break;
//			}
//		}
//		else {
//			break;
//		}
//	}
//
//	// Slopes and gradients along the free energy curve to help visualize the curve's details better
//	for (i = 0; i < freeEnergyArraySize; i++) {
//		freeEnergySlopeArray[i] = gradient(freeEnergyArray, i);
//		freeEnergyCurvatureArray[i] = laplace(freeEnergyArray, i);
//		if (trackMinimumCurvature &&
//				(i > 0) &&
//				(freeEnergyCurvatureArray[i] > 0) &&
//				(freeEnergyCurvatureArray[i-1] < 0)) { // turn off flag to track minimum curvature after passing from - to + for the first time (avoids the discontinuity)
//			trackMinimumCurvature = 0;
//		}
//		if (trackMinimumCurvature &&
//				(freeEnergyCurvatureArray[i] < minimumCurvature)) {
//			minimumCurvature = freeEnergyCurvatureArray[i];
//			minimumCurvatureIndex = i;
//		}
//		fprintf(free_energy_slopes, "%i\t%.15f\n", i, freeEnergySlopeArray[i]);
//		fprintf(free_energy_curvature, "%i\t%.15f\n", i, freeEnergyCurvatureArray[i]);
//		if (i == rhoTail) {
//			break;
//		}
//	}
//
//	// Determine the phase-separated densities that minimize the free energy for this curve
//	for (i = 0; i < minimumCurvatureIndex + 1; i++) {
//		rhoTotal = volumeVapor*rhoArray[minimumCurvatureIndex-i] + volumeLiquid*rhoArray[minimumCurvatureIndex+i];
//		freeEnergyTmp = volumeVapor*freeEnergyArray[minimumCurvatureIndex-i] + volumeLiquid*freeEnergyArray[minimumCurvatureIndex+i];
//		if (freeEnergyTmp < freeEnergyTotal) {
//			freeEnergyTotal = freeEnergyTmp;
//			minRho1 = rhoArray[minimumCurvatureIndex-i];
//			minFreeEnergy1 = freeEnergyArray[minimumCurvatureIndex-i];
//			minRho2 = rhoArray[minimumCurvatureIndex+i];
//			minFreeEnergy2 = freeEnergyArray[minimumCurvatureIndex+i];
//		}
//		printf("N=%f\tF=%f\n", rhoTotal, freeEnergyTmp);
//	}
//
//	fprintf(free_energy_minimum,"%.15f %.15f\n%.15f %.15f\n", minRho1, minFreeEnergy1, minRho2, minFreeEnergy2);
//	printf("free energy array size = %i\n", freeEnergyArraySize);
//	printf("min free energy = %.15f\n", freeEnergyTotal);
//	printf("rho1 = %f\trho2 = %.15f\n", minRho1, minRho2);
//	printf("rho tail = %i\n", rhoTail);
//	printf("min curvature index = %i\tmin curvature=%.15f\n", minimumCurvatureIndex, minimumCurvature);
//
//	fclose(free_energy_results);
//	fclose(free_energy_minimum);
//	fclose(free_energy_slopes);
//	fclose(free_energy_curvature);
//	printf("complete\n\n");
//
//} // end function calculateFreeEnergySingleComponentMinimization()


/**
 * This function calculates the total free energy for a given particles/volume trial pair and a given temperature.
 * The arrays to store trial info are also passed in along with the current index to keep all data in sync.
 */
void calculateFreeEnergyMinimumTrial(double particlesTrial,
									 double volumeTrial,
									 double theta,
									 double *freeEnergyTrialArray,
									 double *rho1TrialArray,
									 double *rho2TrialArray,
									 int trialNumber) {

	double particles1 = particlesTrial;
	double particles2 = 0;
	double particlesTotal = nc;
	double volume1 = volumeTrial;
	double volume2 = 0;
	double volumeTotal = nc;
	double rho1 = 0;
	double rho2 = 0;
	double invalidDensity = 0;
	double freeEnergy1 = 0;
	double freeEnergy2 = 0;
	double invalidFreeEnergy = 1000.* nc;

	if (particles1 <= particlesTotal && volume1 <= volumeTotal) { // do nothing if conservation statements violated
		if (particles1 > 0 && volume1 > 0) {
			rho1 = particles1 / volume1;
		}
		else {
			rho1 = invalidDensity;
		}

		if (rho1 <= nc) { // only look at trial points that are on the left of the critical density to avoid the free energy discontinuity and beyond
			rho1TrialArray[trialNumber] = rho1;

			particles2 = particlesTotal - particles1;	// particle number for second density determined by mass conservation
			volume2 = volumeTotal - volume1;			// volume for second density determined by mass conservation

			if (particles2 > 0 && volume2 > 0) {
				rho2 = particles2 / volume2;
			}
			else {
				rho2 = invalidDensity;
			}
			rho2TrialArray[trialNumber] = rho2;

			if (rho1 > 0 && (rho1*b1 < 1)) {
				freeEnergy1 = rho1*theta*log(rho1/(1.-b1*rho1)) - a1*rho1*rho1;
			}
			else {
				freeEnergy1 = invalidFreeEnergy;
			}
			if (rho2 > 0 && (rho2*b1 < 1)) {
				freeEnergy2 = rho2*theta*log(rho2/(1.-b1*rho2)) - a1*rho2*rho2;
			}
			else {
				freeEnergy2 = invalidFreeEnergy;
			}

			freeEnergyTrialArray[trialNumber] = volume1*freeEnergy1 + volume2*freeEnergy2;
		}
		else {
			rho1TrialArray[trialNumber] = invalidDensity;
			rho2TrialArray[trialNumber] = invalidDensity;
			freeEnergyTrialArray[trialNumber] = invalidFreeEnergy;
		}
	}
	else {
		rho1TrialArray[trialNumber] = invalidDensity;
		rho2TrialArray[trialNumber] = invalidDensity;
		freeEnergyTrialArray[trialNumber] = invalidFreeEnergy;
	}

} // end function calculateFreeEnergyMinimumTrial()


/**
 * This function calculates the total free energy for a given particles/volume trial pair and a given temperature.
 * The arrays to store trial info are also passed in along with the current indices to keep all data in sync.
 * The arrays passed in are 2-dimensional (the other version handles 1-D arrays).
 */
void calculateFreeEnergyMinimumTrial2DArray(double particlesTrial,
									 double volumeTrial,
									 double theta,
									 double **freeEnergyTrialArray,
									 double **rho1TrialArray,
									 double **rho2TrialArray,
									 int rowNumber,
									 int columnNumber) {

	double particles1 = particlesTrial;
	double particles2 = 0;
	double particlesTotal = nc;
	double volume1 = volumeTrial;
	double volume2 = 0;
	double volumeTotal = nc;
	double rho1 = 0;
	double rho2 = 0;
	double invalidDensity = 0;
	double freeEnergy1 = 0;
	double freeEnergy2 = 0;
	double invalidFreeEnergy = 1000.* nc;

	if (particles1 <= particlesTotal && volume1 <= volumeTotal) { // do nothing if conservation statements violated
		if (particles1 > 0 && volume1 > 0) {
			rho1 = particles1 / volume1;
		}
		else {
			rho1 = invalidDensity;
		}

		if (rho1 <= nc) { // only look at trial points that are on the left of the critical density to avoid the free energy discontinuity and beyond
			rho1TrialArray[rowNumber][columnNumber] = rho1;

			particles2 = particlesTotal - particles1;	// particle number for second density determined by mass conservation
			volume2 = volumeTotal - volume1;			// volume for second density determined by mass conservation

			if (particles2 > 0 && volume2 > 0) {
				rho2 = particles2 / volume2;
			}
			else {
				rho2 = invalidDensity;
			}
			rho2TrialArray[rowNumber][columnNumber] = rho2;

			if (rho1 > 0 && (rho1*b1 < 1)) {
				freeEnergy1 = rho1*theta*log(rho1/(1.-b1*rho1)) - a1*rho1*rho1;
			}
			else {
				freeEnergy1 = invalidFreeEnergy;
			}
			if (rho2 > 0 && (rho2*b1 < 1)) {
				freeEnergy2 = rho2*theta*log(rho2/(1.-b1*rho2)) - a1*rho2*rho2;
			}
			else {
				freeEnergy2 = invalidFreeEnergy;
			}

			freeEnergyTrialArray[rowNumber][columnNumber] = volume1*freeEnergy1 + volume2*freeEnergy2;
		}
		else {
			rho1TrialArray[rowNumber][columnNumber] = invalidDensity;
			rho2TrialArray[rowNumber][columnNumber] = invalidDensity;
			freeEnergyTrialArray[rowNumber][columnNumber] = invalidFreeEnergy;
		}
	}
	else {
		rho1TrialArray[rowNumber][columnNumber] = invalidDensity;
		rho2TrialArray[rowNumber][columnNumber] = invalidDensity;
		freeEnergyTrialArray[rowNumber][columnNumber] = invalidFreeEnergy;
	}

} // end function calculateFreeEnergyMinimumTrial2DArray()


/**
 * This function only considers single-component phases without gradient terms in the free energy.
 * This optimizes the version in comments below. (a bit)
 */
void calculatePhaseDiagramRhoVsTempTheoretical() {

	int i = 0;
	int numberOfDensityTrials = 8; // 9 total, but 0 assumed to be the minimum and excluded from the arrays
	int trialNumber = 0;
	int freeEnergyMinimumIndex = 0;
	int numberOfThetaPoints = 0;
	double volumeTotal = 1.;
	double volumeTrial = 0;
	double particlesTotal = 0;
	double particlesTrial = 0;
	double minRho1 = 0;
	double minRho2 = 0;
	double minFreeEnergy = 0;
	double freeEnergyThreshold = 1e-6;
	double thetaTmp = tc;
	double particlesStepSize = 0;
	double thetaStepSize = 0.00001;

	// Arrays to hold the information for each trial density
	double freeEnergyTrialArray[numberOfDensityTrials];
	double rho1TrialArray[numberOfDensityTrials];
	double rho2TrialArray[numberOfDensityTrials];
	double particles1TrialArray[numberOfDensityTrials];
	double volume1TrialArray[numberOfDensityTrials];
	for (i = 0; i < numberOfDensityTrials; i++) { // initialize all to zero
		freeEnergyTrialArray[i] = 0;
		rho1TrialArray[i] = 0;
		rho2TrialArray[i] = 0;
		particles1TrialArray[i] = 0;
		volume1TrialArray[i] = 0;
	}

	// Arrays and dimension to hold the density/temperature pairs
	// Stored for later read-out to customize the order of the output file for plotting
	numberOfThetaPoints = thetaTmp / thetaStepSize; // integer division of doubles may miss a point
	printf("number of theta points = %i\n", numberOfThetaPoints);
	double minRhoArray1[numberOfThetaPoints];
	double minRhoArray2[numberOfThetaPoints];
	double thetaRatioArray[numberOfThetaPoints];
	for (i = 0; i < numberOfThetaPoints; i++) {
		minRhoArray1[i] = 0;
		minRhoArray2[i] = 0;
		thetaRatioArray[i] = 0;
	}
	i = 0; // reset the index to use for the resulting minimization data below


	FILE *phaseDiagram_rhoVsTemp_theoretical;
	FILE *phaseDiagram_densities_theoretical;
	phaseDiagram_rhoVsTemp_theoretical = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/phaseDiagram_rhoVsTemp_theoretical.dat","w");
	phaseDiagram_densities_theoretical = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/phaseDiagram_densities_theoretical.dat","w");

	// To allow a more complete phase diagram in less iterations, this loop varies the actual temp (lowers to zero from the critical temp) instead of the critical temp
	for (thetaTmp = tc; thetaTmp > 0; thetaTmp -= thetaStepSize) {
		freeEnergyMinimumIndex = 0;
		particlesStepSize = 0.25 * nc;														// initial step size of 0.25
		particlesTotal = nc * volumeTotal;													// use the critical density particle number as the total conserved; volume total initialized to 1
		minRho1 = nc;																		// assume the minimum is the critical density as an initial condition (minRho2 already is zero)
		minFreeEnergy = minRho1*thetaTmp*log(minRho1/(1.-b1*minRho1)) - a1*minRho1*minRho1; // initialize the min free energy value to be that of the critical density
		particlesTrial = 0.5 * particlesTotal;												// initialize the first trial values
		volumeTrial = 0.5 * volumeTotal;

		// Each particle number, volume trial pair is numbered.
		// Equal step sizes are used for each of the 9 possible trial pairs.
		// The minimum free energy of the 9 trials is selected, and the step size is halved for the next loop iteration.
		//   - Note: this can be collapsed to a couple embedded for loops, but with a bit of index trickery that isn't clear... sticking with this explicit version
		while (particlesStepSize > freeEnergyThreshold) { //freeEnergyDifference > freeEnergyThreshold) {

			// 0. The center point of the 9 trials
			// This is the assumed minimum, either by initialization or because it was returned as such from a previous loop iteration.
			// Nothing needs to happen with this trial.
			//			particles1 = particlesTrial;
			//			volume1 = volumeTrial;

			// 1. Move just the particle number in + direction
			// Explicitly set this trial number every time to ensure indices correspond properly
			trialNumber = 0;
			particles1TrialArray[trialNumber] = particlesTrial + particlesStepSize;
			volume1TrialArray[trialNumber] = volumeTrial;
			calculateFreeEnergyMinimumTrial(particles1TrialArray[trialNumber],
											volume1TrialArray[trialNumber],
											thetaTmp,
											freeEnergyTrialArray,
											rho1TrialArray,
											rho2TrialArray,
											trialNumber);

			// 2. Move just the particle number in the - direction
			trialNumber += 1;
			particles1TrialArray[trialNumber] = particlesTrial - particlesStepSize;
			volume1TrialArray[trialNumber] = volumeTrial;
			calculateFreeEnergyMinimumTrial(particles1TrialArray[trialNumber],
											volume1TrialArray[trialNumber],
											thetaTmp,
											freeEnergyTrialArray,
											rho1TrialArray,
											rho2TrialArray,
											trialNumber);

			// 3. Move just the volume number in the + direction
			trialNumber += 1;
			particles1TrialArray[trialNumber] = particlesTrial;
			volume1TrialArray[trialNumber] = volumeTrial + particlesStepSize;
			calculateFreeEnergyMinimumTrial(particles1TrialArray[trialNumber],
											volume1TrialArray[trialNumber],
											thetaTmp,
											freeEnergyTrialArray,
											rho1TrialArray,
											rho2TrialArray,
											trialNumber);

			// 4. Move just the volume number in the - direction
			trialNumber += 1;
			particles1TrialArray[trialNumber] = particlesTrial;
			volume1TrialArray[trialNumber] = volumeTrial - particlesStepSize;
			calculateFreeEnergyMinimumTrial(particles1TrialArray[trialNumber],
											volume1TrialArray[trialNumber],
											thetaTmp,
											freeEnergyTrialArray,
											rho1TrialArray,
											rho2TrialArray,
											trialNumber);

			// 5. Move both particles, volume in the + direction
			trialNumber += 1;
			particles1TrialArray[trialNumber] = particlesTrial + particlesStepSize;
			volume1TrialArray[trialNumber] = volumeTrial + particlesStepSize;
			calculateFreeEnergyMinimumTrial(particles1TrialArray[trialNumber],
											volume1TrialArray[trialNumber],
											thetaTmp,
											freeEnergyTrialArray,
											rho1TrialArray,
											rho2TrialArray,
											trialNumber);

			// 6. Move both the particles, volume in the - direction
			trialNumber += 1;
			particles1TrialArray[trialNumber] = particlesTrial - particlesStepSize;
			volume1TrialArray[trialNumber] = volumeTrial - particlesStepSize;
			calculateFreeEnergyMinimumTrial(particles1TrialArray[trialNumber],
											volume1TrialArray[trialNumber],
											thetaTmp,
											freeEnergyTrialArray,
											rho1TrialArray,
											rho2TrialArray,
											trialNumber);

			// 7. Move the particles +, volume -
			trialNumber += 1;
			particles1TrialArray[trialNumber] = particlesTrial + particlesStepSize;
			volume1TrialArray[trialNumber] = volumeTrial - particlesStepSize;
			calculateFreeEnergyMinimumTrial(particles1TrialArray[trialNumber],
											volume1TrialArray[trialNumber],
											thetaTmp,
											freeEnergyTrialArray,
											rho1TrialArray,
											rho2TrialArray,
											trialNumber);

			// 8. Move the particles -, volume +
			trialNumber += 1;
			particles1TrialArray[trialNumber] = particlesTrial - particlesStepSize;
			volume1TrialArray[trialNumber] = volumeTrial + particlesStepSize;
			calculateFreeEnergyMinimumTrial(particles1TrialArray[trialNumber],
											volume1TrialArray[trialNumber],
											thetaTmp,
											freeEnergyTrialArray,
											rho1TrialArray,
											rho2TrialArray,
											trialNumber);

			// Find the minimum of the 8 trials (0 is already the assumed minimum)
			freeEnergyMinimumIndex = findIndexOfMinimumValue(freeEnergyTrialArray, numberOfDensityTrials);
			if (freeEnergyTrialArray[freeEnergyMinimumIndex] < minFreeEnergy) { // set the new minimum values of free energy and densities
				minFreeEnergy = freeEnergyTrialArray[freeEnergyMinimumIndex];
				minRho1 = rho1TrialArray[freeEnergyMinimumIndex];
				minRho2 = rho2TrialArray[freeEnergyMinimumIndex];

				// The next trial values or particle number and volume correspond to the newly found minRho1 value
				particlesTrial = particles1TrialArray[freeEnergyMinimumIndex];
				volumeTrial = volume1TrialArray[freeEnergyMinimumIndex];
			}
			else { // set up for the next loop iteration if necessary by halving the step size
				particlesStepSize *= 0.5;
			}
		}

		//		fprintf(phaseDiagram_rhoVsTemp_theoretical,"%.15f %.15f\n", minRho1, thetaTmp/tc);
		minRhoArray1[i] = minRho1;
		minRhoArray2[i] = minRho2;
		thetaRatioArray[i] = thetaTmp / tc;
		i += 1;
		fprintf(phaseDiagram_densities_theoretical,"%.15f\n%.15f\n", minRho1, minRho2); // log densities to create other phase diagrams
		printf("next tc...%f\n", thetaTmp);
	}

	// Read out density and temperature data in a format to let gnuplot plot the theory curve with a line
	for (i = 0; i < numberOfThetaPoints; i++) { // read minRho1 backwards, theta ratio backwards
		fprintf(phaseDiagram_rhoVsTemp_theoretical,"%.15f %.15f\n", minRhoArray1[numberOfThetaPoints-1-i], thetaRatioArray[numberOfThetaPoints-1-i]);
	}
	for (i = 0; i < numberOfThetaPoints; i++) { // read minRho2 forwards, theta ratio forward
		fprintf(phaseDiagram_rhoVsTemp_theoretical,"%.15f %.15f\n", minRhoArray2[i], thetaRatioArray[i]);
	}

	fclose(phaseDiagram_rhoVsTemp_theoretical);
	fclose(phaseDiagram_densities_theoretical);
	printf("complete\n\n");

} // end function calculatePhaseDiagramRhoVsTempTheoretical()


/**
 * This function calculates the theoretical density vs. pressure phase diagram curve for a given temperature (theta).
 * It is meant to depend on the densities produced by calculatePhaseDiagramRhoVsTempTheoretical() for the density vs. temp phase diagram.
 */
void calculatePhaseDiagramRhoVsPressureTheoretical() {

	int readEOF = 0;
	double tmpTheta = theta;
	double tmpRho = 0;
	double tmpPressure = 0;

	FILE *phaseDiagram_rhoVsPressure_theoretical; // write new pressure data here
	FILE *phaseDiagram_densities_theoretical; // read the densities from this file
	phaseDiagram_rhoVsPressure_theoretical = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/phaseDiagram_rhoVsPressure_theoretical.dat","w");
	phaseDiagram_densities_theoretical = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/phaseDiagram_densities_theoretical.dat","r");

	readEOF = fscanf(phaseDiagram_densities_theoretical, "%lf", &tmpRho);
	while (readEOF != EOF) {
		printf("next density...%.15f\n", tmpRho);

		tmpPressure = tmpRho*tmpTheta + (b1*tmpRho*tmpRho*tmpTheta)/(1.-b1*tmpRho) - a1*tmpRho*tmpRho;
		fprintf(phaseDiagram_rhoVsPressure_theoretical,"%.15f %.15f\n", tmpRho, tmpPressure);
		readEOF = fscanf(phaseDiagram_densities_theoretical, "%lf", &tmpRho);
	}

	fclose(phaseDiagram_rhoVsPressure_theoretical);
	fclose(phaseDiagram_densities_theoretical);
	printf("complete\n\n");

} // end function calculatePhaseDiagramRhoVsPressureTheoretical()


/**
 * This function calculates the theoretical density vs. chemical potential phase diagram curve for a given temperature (theta).
 * It is meant to depend on the densities produced by calculatePhaseDiagramRhoVsTempTheoretical() for the density vs. temp phase diagram.
 */
void calculatePhaseDiagramRhoVsMuTheoretical() {

	int readEOF = 0;
	double tmpTheta = theta;
	double tmpRho = 0;
	double tmpMu = 0;

	FILE *phaseDiagram_rhoVsMu_theoretical; // write new pressure data here
	FILE *phaseDiagram_densities_theoretical; // read the densities from this file
	phaseDiagram_rhoVsMu_theoretical = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/phaseDiagram_rhoVsMu_theoretical.dat","w");
	phaseDiagram_densities_theoretical = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/phaseDiagram_densities_theoretical.dat","r");

	readEOF = fscanf(phaseDiagram_densities_theoretical, "%lf", &tmpRho);
	while (readEOF != EOF) {
		printf("next density...%.15f\n", tmpRho);

		tmpMu = tmpTheta*log(tmpRho/(1.-b1*tmpRho)) + tmpTheta*b1*tmpRho/(1.-b1*tmpRho) + tmpTheta - 2.*a1*tmpRho;
		fprintf(phaseDiagram_rhoVsMu_theoretical,"%.15f %.15f\n", tmpRho, tmpMu);
		readEOF = fscanf(phaseDiagram_densities_theoretical, "%lf", &tmpRho);
	}

	fclose(phaseDiagram_rhoVsMu_theoretical);
	fclose(phaseDiagram_densities_theoretical);
	printf("complete\n\n");

} // end function calculatePhaseDiagramRhoVsMuTheoretical()


/**
 * This function measures the performance of the LB code by generating a density vs temp phase diagram.
 */
void getPhaseDiagramDensityVsTemp() {
	int i = 0;
	int ip, im;

	FILE *phase_results;
	phase_results = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/single component/phase-data-4/phase-data.dat","w");

	// Loop through a range of critical temperatures (tc)
	for (tc = theta+0.0005; theta/tc > quenchDepth; tc += 0.005) { //0.0005
		printf("\ntemp ratio = %f\ttc = %f\t", theta/tc, tc);

		setInitializeSteps();

		// At each tc, iterate enough times for the phase to settle
		for (i=0; i < phase_iterations; i++) {
			iteration();
		}

		// After phases settle, identify the min/max values across the lattice
		double max = n1[0];
		double min = n1[0];
		for (i = 0; i < XDIM; i++) {
			if (n1[i] > max) {
				max = n1[i];
			}
			if (n1[i] < min) {
				min = n1[i];
			}

			ip = i+1;
			if (ip == XDIM) {
				ip = XDIM - 1;
			}
			im = i-1;
			if (im == -1) {
				im = 0;
			}

		}

		fprintf(phase_results,"%.15f %.15f\n%.15f %.15f\n", min, theta/tc, max, theta/tc);
		Events(1);
		DrawGraphs();
	}

	fclose(phase_results);

} // end function getPhaseDiagramDensityVsTemp()


void getDensityProfile() {
	int i = 0;

	FILE *density_profile;
	density_profile = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/density_profile.dat","w");

	for (i = 0; i < XDIM; i++) {
		fprintf(density_profile, "%f\n", n1[i]);
	}

	fclose(density_profile);
	printf("density profile done\n");
}


void getDiffusionSeries() {
	int i = 0, j = 0;
	int iterations = 30000;
	int time_steps = 10;
	char name[500];

	FILE *diffusion_series;

	// Loop through a range of time steps
	for (i = 0; i < time_steps+10; i++) {

		if (i == 10) {
			wall *= -1;
		}

		// At each time step, iterate enough times to propagate the diffusion
		for (j = 0; j < iterations; j++) {
			iteration();
		}

		if (i >= 9) {
			sprintf(name,"/home/clark/school/Lattice Boltzmann/data/diffusion-series_%01i.dat",i-10);
			diffusion_series = fopen(name,"w");

			// Log the n1 lattice values
			for (j = 0; j < XDIM; j++) {
				fprintf(diffusion_series, "%i %e\n",j, n1[j]);
			}

			fprintf(diffusion_series, "\n");
			fclose(diffusion_series);
		}
		Events(1);
		DrawGraphs();
	}

} // end function getDiffusionSeries()


void getRateOfDiffusion() {
	char name[500];
	char append_time[100];
	time_t rawTime;
	struct tm *time_struct;
	double midpoint;
	int i, front, old_wall;

	midpoint = n1_liquid - (n1_liquid-n1_gas)*.5; // half of initial liquid-gas densities
	old_wall = wall;
	collectData = 1; // start data collection

	FILE *rate_of_diffusion;
	time(&rawTime);
	time_struct = localtime(&rawTime);
	strftime(append_time, sizeof(append_time), "%Y%m%d%H%M%S", time_struct);
	sprintf(name,"/home/clark/school/Lattice Boltzmann/data/rate-of-diffusion_%s.dat", append_time);
	rate_of_diffusion = fopen(name,"w");

	printf("midpoint = %f\n", midpoint);
	front = wall; // the diffusion front starts at the wall
	wall *= -1; // remove the wall immediately to begin the diffusion

	// Collect diffusion data until commanded to stop
	while (collectData) {
		iteration();

		// Log the lattice spot of the center of the diffusion front
		//  Only need to look left of original wall location
		for (i = 0; i <= old_wall; i++) {

			// Densities decrease from left to right
			//  - the first density less than middle value is the transition
			//  - log the time/total iterations and previous location (when/where the front is)
			if (n1[i+1] < midpoint && i < front) {
				fprintf(rate_of_diffusion, "%i %i\n", iterations, i);
				front = i; // track the new front location
				break;
			}
		} // end for

		Events(1);
		DrawGraphs();
	} // end while

	fclose(rate_of_diffusion);

} // end function getRateOfDiffusion()


FILE * openFile(char *fileName, char *extension, char *attributes) {
	int currentLength = 0;
	char file_name[MAXSTRLEN];
	char *dataDirectory = "/home/clark/school/Lattice Boltzmann/Maxwell Construction/single component/phase-data-4";
	currentLength = snprintf(file_name, MAXSTRLEN, "%s/", dataDirectory);
	currentLength += snprintf(file_name+currentLength, MAXSTRLEN-currentLength, "%s", fileName);
	currentLength += snprintf(file_name+currentLength, MAXSTRLEN-currentLength, "%s", extension);

	printf("Attempting to open %s\n", file_name);

	FILE * file_handle;
	file_handle = fopen(file_name, attributes);
	return file_handle;
} // end function openFile()


double calculateFreeEnergyDensity(double rho, double theta) {
	double a = a1, b = b1;
	double invalidFreeEnergy = 1000. * nc;
	double freeEnergy = invalidFreeEnergy;
	double zeroThreshold = 1e-14;

	if (!isEqual(rho, 0.0, zeroThreshold)) {
		freeEnergy = rho*theta*log(rho/(1.0-rho*b)) - a*rho*rho - theta*rho;
	}
	else {
		freeEnergy = invalidFreeEnergy;
	}

	return freeEnergy;
}


void generateFreeEnergySlices() {
	int readEOF = 0;
	double rho = 0.0;
	double A1 = 0.0, A2 = 0.0;
	double sStep = 0.01;
	double tc = 0.4, theta = 0.4, thetaRatio = 1.0;
	double freeEnergyDensity = 0.0;

	printf("Calculating free energy tie line slices to generate a surface...\n");

	FILE *twoPhase_densities = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/single component/phaseDiagram_rho-rhoVsTemp.dat", "r");
	FILE *free_energy_curves = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/single component/free-energy-curves.dat", "w");

	readEOF = fscanf(twoPhase_densities, "%17lf %17lf %17lf", &A1, &A2, &thetaRatio);
	while (readEOF != EOF) {
		for (double s = -100.0*sStep; s < (1.0 + 100.0*sStep); s += sStep) { // add extra steps on either side of the actual minimum densities to visualize curvature better
			rho = A1 + s*(A2-A1);
			theta = thetaRatio * tc;
			freeEnergyDensity = calculateFreeEnergyDensity(rho, theta);
			if (!isnan(freeEnergyDensity)) {
				fprintf(free_energy_curves, "%f %.15f %.15f\n", s, thetaRatio, freeEnergyDensity);
			}
		}
		readEOF = fscanf(twoPhase_densities, "%17lf %17lf %17lf", &A1, &A2, &thetaRatio);
	}

	fclose(free_energy_curves);
	fclose(twoPhase_densities);
	printf("Surface data created!\n\n");
} // end function generateFreeEnergySlices()
