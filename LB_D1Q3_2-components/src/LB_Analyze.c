/**
 * @file LB_Analyze.c
 * @author Kent S. Ridl
 * @date 7 July 2016 (last modified: 26 November 2017)
 *
 * The module _LB_Analyze.c_ contains the code used to set theoretical expectations for the lattice Boltzmann simulations and to process simulation results.
 * It includes everything from the free energy minimizer, gradient stencils, logging functions, etc.
 */

#include "LB_D1Q3_2-components.h"


/**
 * @brief Function _gradient_ returns a spatial gradient at a specified location.
 *
 * This function calculates a spatial gradient for both periodic and solid boundaries of a simulation domain.
 * It implements 2nd order forward, backward, and central differences.
 * It does not implement a halo to account for gradients at edges of the domain.
 *
 * The global flag _useBoundaryConditionsPeriodic_ toggles between the periodic and solid boundary calculations.
 * - 1: uses periodic boundary conditions (central difference only)
 * - 0: uses solid walls at domain boundaries (forward, backward, central differences allowed)
 *
 * @param [in] *var An input array (pointer) representing values across space.
 * @param [in] i The array index of the location for which the spatial gradient is to be calculated.
 *
 * @return The calculated spatial gradient.
 */
double gradient(double *var, int i) {
	double result = 0;
	int ipp = i+2, ip = i+1; /// use indices for 2 cells before and after cell i
	int im = i-1, imm = i-2;

	if (useBoundaryConditionsPeriodic) { // gradient for a lattice with periodic BCs
		if (ip == XDIM) ip = 0;
		else if (im == -1) im = XDIM - 1;
		result = 1./2. * (var[ip]-var[im]);
	}
	else { // gradient for a lattice with walls at the boundaries; adjust for existence of the walls with forward/backward differencing
		if (ip == XDIM)  result = 1./2. * (var[imm] - 4.*var[im] + 3.*var[i]); // backward difference
		else if (im == -1) result = 1./2. * (-3.*var[i] + 4.*var[ip] - var[ipp]); // forward difference
		else result = 1./2. * (var[ip]-var[im]); // central difference... normal case
	}
	return result;
} // end function gradient()


/**
 * @brief Function _laplace_ returns a spatial laplacian at a specified location.
 *
 * This function calculates a laplacian for both periodic and solid boundaries of a simulation domain.
 * It implements 2nd order forward, backward, and central differences.
 * It does not implement a halo to account for the laplacian at edges of the domain.
 *
 * The global flag _useBoundaryConditionsPeriodic_ toggles between the periodic and solid boundary calculations.
 * - 1: uses periodic boundary conditions (central difference only)
 * - 0: uses solid walls at domain boundaries (forward, backward, central differences allowed)
 *
 * @param [in] *var An input array (pointer) representing values across space.
 * @param [in] i The array index of the location for which the laplacian is to be calculated.
 *
 * @return The calculated laplacian.
 */
double laplace(double *var, int i) {
	double result = 0;
	int ippp = i+3, ipp = i+2, ip = i+1; /// use indices for 3 cells before and after cell i
	int im = i-1, imm = i-2, immm = i-3;

	if (useBoundaryConditionsPeriodic) { // gradient for lattice with periodic BCs; instead of a halo, reset the out-of-bounds index to opposite end the lattice
		if (ip == XDIM) ip = 0;
		else if (im == -1) im = XDIM - 1;
		result = var[ip] - 2.*var[i] + var[im];
	}
	else { // gradient for a lattice with walls at the boundaries; adjust for existence of the walls with forward/backward differencing
		if (ip == XDIM) result = -1.*var[immm] + 4.*var[imm] -5.*var[im] + 2.*var[i]; // backward difference
		else if (im == -1) result = 2.*var[i] - 5.*var[ip] + 4.*var[ipp] - var[ippp]; // forward difference
		else result = var[ip] - 2.*var[i] + var[im]; // central difference... normal case
	}
	return result;
} // end function laplace()


/**
 * @brief Function _stableSimulation_ searches an array for NaN values.
 *
 * This function searches an array to see if any NaN values exist.  The existence of a NaN value is equivalent to an instability.
 * It is typically used to evaluate the results of a lattice Boltzmann simulation as an error-checking mechanism.
 *
 * @note The array passed in must be of size XDIM.
 *
 * @param [in] *var An input array (pointer) representing values across space.
 *
 * @return Integer value indicating a stable (1) or unstable (0) simulation.
 */
int stableSimulation(double *var) {
	for (int i = 0; i < XDIM; i++) {
		if (isnan(var[i])) {
			printf("Instability detected: index %i\n", i);
			return 0;
		}
	}
	return 1;
} // end function stableSimulation()


/**
 * @brief Function _findIndexOfMinimumValue3DArray_ finds the minimum value of a given 3-dimensional array.
 *
 * This function returns the three indices of the array location containing the minimum value of that array.
 * It will only handles square 3-D arrays.
 *
 * @note This function is only used when a free energy minimization allows for only 2 phases in a phase diagram.
 *
 * @param [in] array The input 3-D array.
 * @param [out] *minIIndex Pointer to the variable to hold the I-index for the result.
 * @param [out] *minJIndex Pointer to the variable to hold the J-index for the result.
 * @param [out] *minKIndex Pointer to the variable to hold the K-index for the result.
 */
void findIndexOfMinimumValue3DArray(double array[3][3][3], int *minIIndex, int *minJIndex, int *minKIndex) {
    int currentMinI = 1, currentMinJ = 1, currentMinK = 1; /// assume the center of the 3-D array is the starting minimum

    for (int i = 0; i < 3; i++) {
    	for (int j = 0; j < 3; j++) {
    		for (int k = 0; k < 3; k++) {
    			if (array[i][j][k] < array[currentMinI][currentMinJ][currentMinK]) {
    				currentMinI = i;
    				currentMinJ = j;
    				currentMinK = k;
    			}
    		}
    	}
    }

    *minIIndex = currentMinI;
    *minJIndex = currentMinJ;
    *minKIndex = currentMinK;
} // end function findIndexOfMinimumValue3DArray()


/**
 * @brief Function _findIndexOfMinimumValue6DArray_ finds the minimum value of a given 6-dimensional array.
 *
 * This function returns the three indices of the array location containing the minimum value of that array.
 * It will only handles square 6-D arrays.
 *
 * @note This function is only used when a free energy minimization allows for only 3 phases in a phase diagram.
 *
 * @param [in] array The input 6-D array.
 * @param [out] *minIIndex Pointer to the variable to hold the I-index for the result.
 * @param [out] *minJIndex Pointer to the variable to hold the J-index for the result.
 * @param [out] *minKIndex Pointer to the variable to hold the K-index for the result.
 * @param [out] *minLIndex Pointer to the variable to hold the L-index for the result.
 * @param [out] *minMIndex Pointer to the variable to hold the M-index for the result.
 * @param [out] *minNIndex Pointer to the variable to hold the N-index for the result.
 */
void findIndexOfMinimumValue6DArray(struct FreeEnergy array[3][3][3][3][3][3],
		int *minIIndex, int *minJIndex, int *minKIndex, int *minLIndex, int *minMIndex, int *minNIndex) {
    int currentMinI = 1, currentMinJ = 1, currentMinK = 1, currentMinL = 1, currentMinM = 1, currentMinN = 1; /// assume center of 6-D array is starting minimum

    // TODO: may need threshold to keep junk data from influencing the comparison
    for (int i = 0; i < 3; i++) {
    	for (int j = 0; j < 3; j++) {
    		for (int k = 0; k < 3; k++) {
    		    for (int l = 0; l < 3; l++) {
    		    	for (int m = 0; m < 3; m++) {
    		    		for (int n = 0; n < 3; n++) {
    		    			if (array[i][j][k][l][m][n].freeEnergyTotal <
    		    					array[currentMinI][currentMinJ][currentMinK][currentMinL][currentMinM][currentMinN].freeEnergyTotal) {
    		    				currentMinI = i;
    		    				currentMinJ = j;
    		    				currentMinK = k;
    		    				currentMinL = l;
    		    				currentMinM = m;
    		    				currentMinN = n;
    		    			}
    		    		}
    		    	}
    		    }
    		}
    	}
    }

    *minIIndex = currentMinI;
    *minJIndex = currentMinJ;
    *minKIndex = currentMinK;
    *minLIndex = currentMinL;
    *minMIndex = currentMinM;
    *minNIndex = currentMinN;
} // end function findIndexOfMinimumValue6DArray()


/**
 * @brief Function _findIndicesOfMinMaxValues1DArray_ finds the minimum and maximum values of a given 1-dimensional array.
 *
 * This function returns index of the array location containing the minimum value of that array.
 * This function also returns index of the array location containing the maximum value of that array.
 * It will only handles square 1-D arrays.
 *
 * @note This function supports identification of liquid and vapor phases in the 2-phase regions of a phase diagram.
 *
 * @param [in] array The input 1-D array.
 * @param [out] *minIndex Pointer to the variable to hold the index for the array minimum value.
 * @param [out] *maxIndex Pointer to the variable to hold the index for the array maximum value.
 */
void findIndicesOfMinMaxValues1DArray(double *array, int *minIndex, int *maxIndex) {
	int tmpMinIndex = 0, tmpMaxIndex = 0;

	for (int i = 0; i < XDIM; i++) {
		if (array[i] > array[tmpMaxIndex]) tmpMaxIndex = i;
		else if (array[i] < array[tmpMinIndex]) tmpMinIndex = i;
	}

	*minIndex = tmpMinIndex;
	*maxIndex = tmpMaxIndex;
} // end function findIndicesOfMinMaxValues1DArray()


/**
 * @brief Function _arrayInsertBinodalPoint_ inserts a value into an array at a specified location.
 *
 * The insertion of a new value into an array can happen at any specified location, not just the end.
 * All array elements at and after the insertion location are shifted right by 1 spot,
 * and the new element is inserted at the desired location.
 * The tail of the array is updated to reflect a successful insertion.
 *
 * @note As of now, this function will only operate on the data type _BinodalPoint_.
 * @note This function currently assumes the input is error checked and the array will not overflow.
 *
 * @param [in,out] *array The array (pointer) in which to insert the value.
 * @param [in] insertValue The value to insert into the array.
 * @param [in] insertIndex The index of the array location where the new value is to be inserted.
 * @param [in,out] *arrayTail A pointer to the index of the tail of the array.
 *
 * @todo: Handle failed array insertions.
 * @todo: Update to handle any data type, then refactor function name to just _arrayInsert_
 */
void arrayInsertBinodalPoint(struct BinodalPoint *array, struct BinodalPoint insertValue, int insertIndex, int *arrayTail) {
	for (int i = *arrayTail; i >= insertIndex; i--) array[i+1] = array[i];
	array[insertIndex] = insertValue;
	*arrayTail += 1;
	//TODO: handle failed array insertions
	//TODO: refactor to handle any data type
} // end function arrayInsertBinodalPoint()


/**
 * @brief Function _printRunTime_ formats benchmark times and prints them to the console.
 *
 * This function is used to output benchmark data collected in select functions.
 * The number of processor cycles is scaled by the global constant CLOCKS_PER_SEC.
 * Both the processor time and wall clock time are output to the console in HH:MM:SS format.
 *
 * @param [in] processorRunTime The number of processor clock cycles.
 * @param [in] clockRunTime The wall clock time (seconds).
 */
void printRunTime(clock_t processorRunTime, double clockRunTime) {
	int hours = 0, minutes = 0, seconds = 0;

	processorRunTime /= CLOCKS_PER_SEC;

	hours = processorRunTime / 3600;
	processorRunTime = processorRunTime % 3600;
	minutes = processorRunTime / 60;
	processorRunTime = processorRunTime % 60;
	seconds = processorRunTime;

	printf("\nProcessor time... %.2d:%.2d:%.2d\n", hours, minutes, seconds);

	hours = clockRunTime / 3600;
	clockRunTime = (int) clockRunTime % 3600;
	minutes = clockRunTime / 60;
	clockRunTime = (int) clockRunTime % 60;
	seconds = clockRunTime;

	printf("Clock time... %.2d:%.2d:%.2d\n", hours, minutes, seconds);
} // end function printRunTime()


/**
 * @brief Function _printMinimizationResults_ prints the free energy minimization results for a given test point to the console.
 *
 * @param [in] particlesA The number of A particles in the given test point.
 * @param [in] particlesB the number of B particles in the given test point.
 */
void printMinimizationResults(double particlesA, double particlesB) {
	int readEOF = 0;
	int parentFound = 0;
	double threshold = 1e-6;
	double particlesATotal = 0.0, particlesBTotal = 0.0;
	double rhoA1 = 0.0, rhoA2 = 0.0, rhoA3 = 0.0;
	double rhoB1 = 0.0, rhoB2 = 0.0, rhoB3 = 0.0;

	FILE *two_phase_coords = openFile("threePhaseDiagram-densities-twoPhases", ".dat", "r");
	readEOF = fscanf(two_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
			&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
	while (readEOF != EOF) {
		if ( isEqual(particlesA, particlesATotal, threshold) && isEqual(particlesB, particlesBTotal, threshold) ) {
			printf("\n2-phase point found...\n");
			printf("particlesATotal=%.15f\tparticlesBTotal=%.15f\n", particlesATotal, particlesBTotal);
			printf("%.15f %.15f %.15f %.15f %.15f %.15f\n\n", rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);
			parentFound = 1;
			break;
		}
		else {
			readEOF = fscanf(two_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
					&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
		}
	}
	fclose(two_phase_coords);

	if (!parentFound) {
		FILE *three_phase_coords = openFile("threePhaseDiagram-densities-threePhases", ".dat", "r");
		readEOF = fscanf(three_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
				&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
		while (readEOF != EOF) {
			if ( isEqual(particlesA, particlesATotal, threshold) && isEqual(particlesB, particlesBTotal, threshold) ) {
				printf("\n3-phase point found...\n");
				printf("particlesATotal=%.15f\tparticlesBTotal=%.15f\n", particlesATotal, particlesBTotal);
				printf("%.15f %.15f %.15f %.15f %.15f %.15f\n\n", rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);
				break;
			}
			else {
				readEOF = fscanf(three_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
						&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
			}
		}
		fclose(three_phase_coords);
	}
} // end function printMinimizationResults()


/**
 * @brief Function _calculateMinimizationPath_ assesses the theoretical stability of a given test point.
 *
 * This function performs a stability analysis on a given point in the phase diagram to see if it is expected to phase separate or not.
 * A negative eigenvalue of the Hessian matrix of the free energy indicates the test point is unstable and prone to phase separation.
 * The directions of any instabilities are given by the eigenvectors associated with negative eigenvalues of the Hessian matrix.
 *
 * @note This function is specific to the free energy of a 2-component mixture of van Der Waals (VDW) fluids (i.e. 2x2 Hessian matrix, the determinant of
 * which is the product of the 2 real eigenvalues).
 * @note The eigenvectors associated with the eigenvalues in this implementation always have a B-coordinate equal to 1.0.
 * @note The output eigenvector is normalized by the vector magnitude.
 * @note Testing with the 2-component VDW free energy showed only one negative eigenvalue in each test case, which is why only one eigenvector is output.
 *
 * @param [in] rhoA The A-coordinate of the test point.
 * @param [in] rhoB The B-coordinate of the test point.
 * @param [out] *eigenvectorA Pointer to the A-coordinate of the instability direction.
 * @param [out] *eigenvectorB Pointer to the B-coordinate of the instability direction.
 *
 * @return Integer values of 1 and 0 respectively correspond to an unstable or stable test point.
 *
 * @todo This function may be further generalized to return more than one eigenvector with variable B-coordinates.
 */
int calculateMinimizationPath(double rhoA, double rhoB, double *eigenvectorA, double *eigenvectorB) {
	double D = 1.0 - bA*rhoA - bB*rhoB; // denominator common to all terms below
	double a = theta*(1.0-bB*rhoB) / (rhoA*D) + theta*bA*(1.0-bB*rhoB+bA*rhoB) / (D*D) - 2.0*aA;
	double b = theta*(1.0-bA*rhoA) / (rhoB*D) + theta*bB*(1.0-bA*rhoA+bB*rhoA) / (D*D) - 2.0*aB;
	double c = theta*(bA + bB - bA*bA*rhoA - bB*bB*rhoB) / (D*D) - 2.0*aAB;
	double determinant = a*b - c*c;
	double eigenvalue1 = 0.5 * (a + b + sqrt((a-b)*(a-b) + 4*c*c));
	double eigenvector1 = 2.0*c / (b - a + sqrt((a-b)*(a-b) + 4*c*c));
	double eigenvector2 = 2.0*c / (b - a - sqrt((a-b)*(a-b) + 4*c*c));
	double magnitude = 0.0;
	double zeroThreshold = 1e-14;

	/// Eigenvalue2 always seems to be the negative one, so default to that eigenvector for the instability direction
	if (eigenvalue1 < 0) {
		magnitude = sqrt(eigenvector1*eigenvector1 + 1.0);
		*eigenvectorA = eigenvector1 / magnitude;
	}
	else {
		magnitude = sqrt(eigenvector2*eigenvector2 + 1.0);
		*eigenvectorA = eigenvector2 / magnitude;
	}
	*eigenvectorB = 1.0 / magnitude;

	// Zero density protection... movement along the A- or B-axis expected
	if (isEqual(rhoA, 0.0, zeroThreshold)) {
		*eigenvectorA = 0.0;
		*eigenvectorB = 1.0;
	}
	if (isEqual(rhoB, 0.0, zeroThreshold)) {
		*eigenvectorA = 1.0;
		*eigenvectorB = 0.0;
	}

	/// The sign of the determinant of a 2x2 Hessian has the signs of the 2 eigenvalues encoded in it
	/// Checking eigenvalue1 catches when the determinant may be positive but both eigenvalues are negative (have never seen this, though)
	return determinant < 0 ? 1 : (eigenvalue1 < 0 ? 1 : 0);
} // end function calculateMinimizationPath


/**
 * @brief Function _generateEigenvalueMap_ generates stability analysis data for a free energy.
 *
 * This function calculates eigenvalues, eigenvectors, a determinant, and elements of the Hessian matrix used in a stability analysis.
 *
 * @note This function is specific to the free energy of a 2-component mixture of van Der Waals (VDW) fluids (i.e. 2x2 Hessian matrix).
 * @note Only the A-coordinates of the eigenvectors are output; the B-coordinates are implicitly equal to 1.0.
 */
void generateEigenvalueMap() {
	double rhoA = 0.0, rhoB = 0.0;
	double aComponentMaximum = 1. / bA, bComponentMaximum = 1. / bB;
	double a = 0, b = 0, c = 0, D = 0; // D is denominator common to all terms below
	double determinant = 0, eigenvalue1 = 0, eigenvalue2 = 0, eigenvector1 = 0, eigenvector2 = 0;

	int spinodalArrayASize = (int) round(aComponentMaximum / minimizationParticlesStepSize);
	int spinodalArrayBSize = (int) round(bComponentMaximum / minimizationParticlesStepSize);

	FILE *eigenstuff_data = openFile("eigenstuff", ".dat", "w");
	for (int j = spinodalArrayBSize-1; j >= 0; j--) {
		for (int i = 0; i < spinodalArrayASize; i++) {

			// Determine total particles based on the indices
			rhoA = (double) i * minimizationParticlesStepSize;
			rhoB = (double) j * minimizationParticlesStepSize;

			// Loops are over a square grid; filter out the invalid region of the phase diagram (other side of VDW singularity line)
			if ((0-spinodalArrayASize)*j  - spinodalArrayBSize*(i-spinodalArrayASize) <= 0) continue;

			D = 1.0 - bA*rhoA - bB*rhoB;
			a = theta*(1.0-bB*rhoB) / (rhoA*D) + theta*bA*(1.0-bB*rhoB+bA*rhoB) / (D*D) - 2.0*aA;
			b = theta*(1.0-bA*rhoA) / (rhoB*D) + theta*bB*(1.0-bA*rhoA+bB*rhoA) / (D*D) - 2.0*aB;
			c = theta*(bA + bB - bA*bA*rhoA - bB*bB*rhoB) / (D*D) - 2.0*aAB;
			determinant = a*b - c*c;
			eigenvalue1 = 0.5*(a + b + sqrt((a-b)*(a-b) + 4*c*c));
			eigenvalue2 = 0.5*(a + b - sqrt((a-b)*(a-b) + 4*c*c));
			eigenvector1 = 2*c/(b-a+sqrt((a-b)*(a-b)+4*c*c));
			eigenvector2 = 2*c/(b-a-sqrt((a-b)*(a-b)+4*c*c));
			fprintf(eigenstuff_data, "%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n",
					rhoA, rhoB, determinant, eigenvalue1, eigenvalue2, eigenvector1, eigenvector2, a, b, c);
		}
	}
	fclose(eigenstuff_data);

	// Need to sed-ify the eigenstuff to replace inf and nan values with ? for the gnuplot scripts
	char sed_eigenstuff_data[255];
	sprintf(sed_eigenstuff_data,"sed -i 's/-\\?\\(nan\\|inf\\)/?/ig' \"%s/eigenstuff_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\"", dataDirectory, tcA, tcB, aA, aAB, aB);
	system(sed_eigenstuff_data);
	gnuplotEigenstuffData();
	printf("Eigenstuff generated!\n\n");
} // end function generateEigenvalueMap()


/**
 * @brief Function _checkVaporPhase_ enforces ordering to the phases in the various log files and initialization functions.
 *
 * This function accepts pointers to a set of densities and reorders them to result in the following:
 * - Phase 1 is the A-rich phase.
 * - Phase 2 is the B-rich phase.
 * - Phase 3 is the vapor phase.
 *
 * In the case where the resulting A-rich phase also happens to have the highest concentration of B particles and there are more total B particles
 * than total A particles, phases 1 and 2 will be flipped to ensure phase 2 has the highest concentration of B particles.  In other isolated cases
 * where the A-rich and B-rich phases are practically equal, phase 3 will be swapped with either phase 1 or 2 depending on if there are more total
 * A or B particles, respectively.
 *
 * @note The ordering imposed by this function is a prerequisite for good behavior in almost every function that operates based on reading in data.
 *
 * @param [in,out] *rhoA1 Pointer to the A component density in phase 1.
 * @param [in,out] *rhoB1 Pointer to the B component density in phase 1.
 * @param [in,out] *rhoA2 Pointer to the A component density in phase 2.
 * @param [in,out] *rhoB2 Pointer to the B component density in phase 2.
 * @param [in,out] *rhoA3 Pointer to the A component density in phase 3.
 * @param [in,out] *rhoB3 Pointer to the B component density in phase 3.
 * @param [in] particlesATotal Total number of A particles associated with a set of densities.
 * @param [in] particlesBTotal Total number of B particles associated with a set of densities.
 */
void checkVaporPhase(double *rhoA1, double *rhoB1, double *rhoA2, double *rhoB2, double *rhoA3, double *rhoB3, double particlesATotal, double particlesBTotal) {
	double tmpDouble = 0.0;
	double threshold = 1e-4;

	if (isnan(*rhoA1) && isnan(*rhoB1)) { // swap 1 & 3 if phase 1 are NaNs... rest of logic should order 1 and 2 correctly then
		tmpDouble = *rhoA3;
		*rhoA3 = *rhoA1;
		*rhoA1 = tmpDouble;

		tmpDouble = *rhoB3;
		*rhoB3 = *rhoB1;
		*rhoB1 = tmpDouble;
	}
	if (*rhoA3 > *rhoA2) { // bubble the highest concentration of A particles up to phase 1
		tmpDouble = *rhoA2;
		*rhoA2 = *rhoA3;
		*rhoA3 = tmpDouble;

		tmpDouble = *rhoB2;
		*rhoB2 = *rhoB3;
		*rhoB3 = tmpDouble;
	}
	if (*rhoA2 > *rhoA1) {
		tmpDouble = *rhoA2;
		*rhoA2 = *rhoA1;
		*rhoA1 = tmpDouble;

		tmpDouble = *rhoB2;
		*rhoB2 = *rhoB1;
		*rhoB1 = tmpDouble;
	}
	if (*rhoB3 > *rhoB2) { // bubble the higher B concentration down to up 2
		tmpDouble = *rhoA2;
		*rhoA2 = *rhoA3;
		*rhoA3 = tmpDouble;

		tmpDouble = *rhoB2;
		*rhoB2 = *rhoB3;
		*rhoB3 = tmpDouble;
	}

	// In isolated cases where the largest A and largest B values are both in phase 1, double check
	// to make sure the largest B is in phase 2, but only when there are more B particles present
	if ( (particlesBTotal > particlesATotal) && (*rhoB1 > *rhoB2) ) {
		tmpDouble = *rhoA2;
		*rhoA2 = *rhoA1;
		*rhoA1 = tmpDouble;

		tmpDouble = *rhoB2;
		*rhoB2 = *rhoB1;
		*rhoB1 = tmpDouble;
	}

	// In isolated cases where the A and B values are essentially equal in both in phase 1 and phase 2,
	// swap phases depending on total number of particles to ensure a valid vapor is available for initialization
	if (isEqual(*rhoA1, *rhoA2, threshold) && isEqual(*rhoB1, *rhoB2, threshold)) {
		if (particlesATotal >= particlesBTotal) {
			tmpDouble = *rhoA2;
			*rhoA2 = *rhoA3;
			*rhoA3 = tmpDouble;

			tmpDouble = *rhoB2;
			*rhoB2 = *rhoB3;
			*rhoB3 = tmpDouble;
		}
		else {
			tmpDouble = *rhoA1;
			*rhoA1 = *rhoA3;
			*rhoA3 = tmpDouble;

			tmpDouble = *rhoB1;
			*rhoB1 = *rhoB3;
			*rhoB3 = tmpDouble;
		}
	}
} // end function checkVaporPhase()


/**
 * @brief Function _calculateFreeEnergyMinimumTwoPhaseTrial_ calculates the free energy of a 2-component mixture allowing for 2 phases.
 *
 * This function calculates the free energy of a mixture of 2 components, assuming only 2 phases are possible.
 * It conserves both particle counts and volume by using the input parameters to calculate the particles and volume for phase 2.
 * Particle counts and volumes are checked for validity to prevent the free energy logarithm from exploding.
 * - When a valid free energy cannot be calculated for a phase, an invalid free energy value is used for that phase.
 * - When only a single component exists, the the Helmholtz free energy is calculated as a van der Waals fluid for only that component.
 * - When A and B components are present, the Helmholtz free energy is calculated as a mixture of van der Waals fluids with appropriate component interactions.
 *
 * @note The value (1000 x component A critical density) is assumed to be a sufficiently large invalid free energy so as to not affect the minimization.
 * @note The free energy calculated in this function is a total free energy, not a free energy density.
 * @note This function was used as an intermediate step only while developing the 3-phase algorithm.  It has not been maintained and is kept for posterity.
 *
 * @param [in] theta The temperature of the mixture.
 * @param [in] particlesA1 The number of A particles in phase 1.
 * @param [in] particlesB1 The number of B particles in phase 1.
 * @param [in] volume1 The volume occupied by phase 1.
 * @param [in] particlesATotal The total number of A particles in the mixture.
 * @param [in] particlesBTotal The total number of B particles in the mixture.
 * @param [in] volumeTotal The total volume occupied by the mixture.
 *
 * @return The total free energy (i.e. the sum of the free energy of phase 1 and free energy of phase 2).
 */
double calculateFreeEnergyMinimumTwoPhaseTrial(double theta, double particlesA1, double particlesB1, double volume1,
												double particlesATotal, double particlesBTotal, double volumeTotal){
	double volume2 = volumeTotal - volume1;
	double particlesA2 = particlesATotal - particlesA1;
	double particlesB2 = particlesBTotal - particlesB1;
	double excludedVolume1 = 0, excludedVolume2 = 0;
	double freeEnergy1 = 0, freeEnergy2 = 0, invalidFreeEnergy = 1000.* ncA;

	excludedVolume1 = particlesA1*bA + particlesB1*bB;
	if ((excludedVolume1 >= volume1) || (particlesA1 < 0) || (particlesB1 < 0) || (volume1 < 0)) freeEnergy1 = invalidFreeEnergy;
	else if ((particlesA1 == 0) && (particlesB1 == 0)) freeEnergy1 = 0;
	else if (particlesA1 == 0) freeEnergy1 = particlesB1*theta*log(particlesB1/(volume1-excludedVolume1)) - aB*particlesB1*particlesB1/volume1;
	else if (particlesB1 == 0) freeEnergy1 = particlesA1*theta*log(particlesA1/(volume1-excludedVolume1)) - aA*particlesA1*particlesA1/volume1;
	else {
		freeEnergy1 = particlesA1*theta*log(particlesA1/(volume1-excludedVolume1)) + particlesB1*theta*log(particlesB1/(volume1-excludedVolume1))
		- aA*particlesA1*particlesA1/volume1 - 2*aAB*particlesB1*particlesA1/volume1 - aB*particlesB1*particlesB1/volume1;
	}

	excludedVolume2 = particlesA2*bA + particlesB2*bB;
	if ((excludedVolume2 >= volume2) || (particlesA2 < 0) || (particlesB2 < 0) || (volume2 < 0)) freeEnergy2 = invalidFreeEnergy;
	else if ((particlesA2 == 0) && (particlesB2 == 0)) freeEnergy2 = 0;
	else if (particlesA2 == 0) freeEnergy2 = particlesB2*theta*log(particlesB2/(volume2-excludedVolume2)) - aB*particlesB2*particlesB2/volume2;
	else if (particlesB2 == 0) freeEnergy2 = particlesA2*theta*log(particlesA2/(volume2-excludedVolume2)) - aA*particlesA2*particlesA2/volume2;
	else {
		freeEnergy2 = particlesA2*theta*log(particlesA2/(volume2-excludedVolume2)) + particlesB2*theta*log(particlesB2/(volume2-excludedVolume2))
		- aA*particlesA2*particlesA2/volume2 - 2*aAB*particlesB2*particlesA2/volume2 - aB*particlesB2*particlesB2/volume2;
	}

	return freeEnergy1 + freeEnergy2;
} // end function calculateFreeEnergyMinimumTwoPhaseTrial()


/**
 * @brief Function _calculateFreeEnergyMinimumThreePhaseTrial_ calculates the free energy of a 2-component mixture allowing for 3 phases.
 *
 * This function calculates the free energy of a mixture of 2 components, assuming 3 phases are possible.
 * It conserves both particle counts and volume by using the input parameters to calculate the particles and volume for phase 3.
 * Particle counts and volumes are checked for validity to prevent the free energy logarithm from exploding.
 * - When only a single component exists, the the Helmholtz free energy is calculated as a van der Waals fluid for only that component.
 * - When A and B components are present, the Helmholtz free energy is calculated as a mixture of van der Waals fluids with appropriate component interactions.
 *
 * @note The value (1000 x component A critical density) is assumed to be a sufficiently large invalid free energy so as to not affect the minimization.
 * @note The free energy calculated is the total free energy summed over the specific volumina in each allowed phase (not free energy densities).
 *
 * @param [in] theta The temperature of the mixture.
 * @param [in] particlesA1 The number of A particles in phase 1.
 * @param [in] particlesB1 The number of B particles in phase 1.
 * @param [in] volume1 The volume occupied by phase 1.
 * @param [in] particlesA2 The number of A particles in phase 2.
 * @param [in] particlesB2 The number of B particles in phase 2.
 * @param [in] volume2 The volume occupied by phase 2.
 * @param [in] particlesATotal The total number of A particles in the mixture.
 * @param [in] particlesBTotal The total number of B particles in the mixture.
 * @param [in] volumeTotal The total volume occupied by the mixture.
 * @param [out] *freeEnergy1 Pointer to store the free energy of phase 1.
 * @param [out] *freeEnergy2 Pointer to store the free energy of phase 2.
 * @param [out] *freeEnergy3 Pointer to store the free energy of phase 3.
 *
 * @return The total free energy (i.e. the sum of *freeEnergy1, *freeEnergy2, and *freeEnergy3).
 */
double calculateFreeEnergyMinimumThreePhaseTrial(double theta,
												 double particlesA1, double particlesB1, double volume1,
												 double particlesA2, double particlesB2, double volume2,
												 double particlesATotal, double particlesBTotal, double volumeTotal,
												 double *freeEnergy1, double *freeEnergy2, double *freeEnergy3){
	double volume3 = volumeTotal - volume1 - volume2;
	double particlesA3 = particlesATotal - particlesA1 - particlesA2;
	double particlesB3 = particlesBTotal - particlesB1 - particlesB2;
	double excludedVolume1 = 0, excludedVolume2 = 0, excludedVolume3 = 0;
	double invalidFreeEnergy = 1000.* ncA;
	double threshold = 1e-14;

	*freeEnergy1 = invalidFreeEnergy;
	excludedVolume1 = particlesA1*bA + particlesB1*bB;
	if (isEqual(particlesA1, 0.0, threshold) && isEqual(particlesB1, 0.0, threshold) && isEqual(volume1, 0.0, threshold)) *freeEnergy1 = 0;
	else if (!isEqual(volume1, 0.0, threshold) && volume1 > 0 && !isEqual(volume1, excludedVolume1, threshold) && volume1 > excludedVolume1) {
		if (isEqual(particlesA1, 0.0, threshold) && !isEqual(particlesB1, 0.0, threshold) && particlesB1 > 0) {
			*freeEnergy1 = particlesB1*theta*log(particlesB1/(volume1-excludedVolume1)) - aB*particlesB1*particlesB1/volume1 - theta*particlesB1;
		}
		else if (isEqual(particlesB1, 0.0, threshold) && !isEqual(particlesA1, 0.0, threshold) && particlesA1 > 0) {
			*freeEnergy1 = particlesA1*theta*log(particlesA1/(volume1-excludedVolume1)) - aA*particlesA1*particlesA1/volume1 - theta*particlesA1;
		}
		else if (!isEqual(particlesA1, 0.0, threshold) && particlesA1 > 0 && !isEqual(particlesB1, 0.0, threshold) && particlesB1 > 0) {
			*freeEnergy1 = particlesA1*theta*log(particlesA1/(volume1-excludedVolume1)) + particlesB1*theta*log(particlesB1/(volume1-excludedVolume1))
			- aA*particlesA1*particlesA1/volume1 - 2*aAB*particlesB1*particlesA1/volume1 - aB*particlesB1*particlesB1/volume1 - theta*(particlesA1+particlesB1);
		}
	}

	*freeEnergy2 = invalidFreeEnergy;
	excludedVolume2 = particlesA2*bA + particlesB2*bB;
	if (isEqual(particlesA2, 0.0, threshold) && isEqual(particlesB2, 0.0, threshold) && isEqual(volume2, 0.0, threshold)) *freeEnergy2 = 0;
	else if (!isEqual(volume2, 0.0, threshold) && volume2 > 0 && !isEqual(volume2, excludedVolume2, threshold) && volume2 > excludedVolume2) {
		if (isEqual(particlesA2, 0.0, threshold) && !isEqual(particlesB2, 0.0, threshold) && particlesB2 > 0) {
			*freeEnergy2 = particlesB2*theta*log(particlesB2/(volume2-excludedVolume2)) - aB*particlesB2*particlesB2/volume2 - theta*particlesB2;
		}
		else if (isEqual(particlesB2, 0.0, threshold) && !isEqual(particlesA2, 0.0, threshold) && particlesA2 > 0) {
			*freeEnergy2 = particlesA2*theta*log(particlesA2/(volume2-excludedVolume2)) - aA*particlesA2*particlesA2/volume2 - theta*particlesA2;
		}
		else if (!isEqual(particlesA2, 0.0, threshold) && particlesA2 > 0 && !isEqual(particlesB2, 0.0, threshold) && particlesB2 > 0) {
			*freeEnergy2 = particlesA2*theta*log(particlesA2/(volume2-excludedVolume2)) + particlesB2*theta*log(particlesB2/(volume2-excludedVolume2))
			- aA*particlesA2*particlesA2/volume2 - 2*aAB*particlesB2*particlesA2/volume2 - aB*particlesB2*particlesB2/volume2 - theta*(particlesA2+particlesB2);
		}
	}

	*freeEnergy3 = invalidFreeEnergy;
	excludedVolume3 = particlesA3*bA + particlesB3*bB;
	if (isEqual(particlesA3, 0.0, threshold) && isEqual(particlesB3, 0.0, threshold) && isEqual(volume3, 0.0, threshold)) *freeEnergy3 = 0;
	else if (!isEqual(volume3, 0.0, threshold) && volume3 > 0 && !isEqual(volume3, excludedVolume3, threshold) && volume3 > excludedVolume3) {
		if (isEqual(particlesA3, 0.0, threshold) && !isEqual(particlesB3, 0.0, threshold) && particlesB3 > 0) {
			*freeEnergy3 = particlesB3*theta*log(particlesB3/(volume3-excludedVolume3)) - aB*particlesB3*particlesB3/volume3 - theta*particlesB3;
		}
		else if (isEqual(particlesB3, 0.0, threshold) && !isEqual(particlesA3, 0.0, threshold) && particlesA3 > 0) {
			*freeEnergy3 = particlesA3*theta*log(particlesA3/(volume3-excludedVolume3)) - aA*particlesA3*particlesA3/volume3 - theta*particlesA3;
		}
		else if (!isEqual(particlesA3, 0.0, threshold) && particlesA3 > 0 && !isEqual(particlesB3, 0.0, threshold) && particlesB3 > 0) {
			*freeEnergy3 = particlesA3*theta*log(particlesA3/(volume3-excludedVolume3)) + particlesB3*theta*log(particlesB3/(volume3-excludedVolume3))
			- aA*particlesA3*particlesA3/volume3 - 2*aAB*particlesB3*particlesA3/volume3 - aB*particlesB3*particlesB3/volume3 - theta*(particlesA3+particlesB3);
		}
	}

	return *freeEnergy1 + *freeEnergy2 + *freeEnergy3;
} // end function calculateFreeEnergyMinimumThreePhaseTrial()


/**
 * @brief Function _minimizeFreeEnergyTwoComponentsThreePhases_ determines the configuration that minimizes the free energy at a given test point.
 *
 * This function will determine the configuration of particles and volumes among 3 phases that minimizes the free energy of a given test point
 * from a theoretical phase diagram.  An initial allocation of particles and volumes among the 3 phases must be provided to start the minimization.
 * The minimization is governed by a loop that iterates until reaching a specified convergence threshold which dictates when the free energy
 * minimization is "good enough."  Within each loop iteration, a matrix is constructed that holds the trial values for each physical degree of
 * freedom that are used to calculate the candidate free energies.  Each physical value is modified by a defined step above and below a reference
 * value.  The configuration of particles and volumes that provides the minimal free energy value is returned.
 * - If the free energy resulting from a minimization is invalid, the loop is broken; however, the invalid values will be returned.
 * - If the loop iteration results in a new minimal free energy, the loop continues to minimize with the new configuration; the step size remains unchanged.
 * - If the loop iteration results in an unchanged minimal free energy, the step size is halved for the next iteration.
 * - The loop terminates when the trial step size for physical quantities is less than a given threshold.
 *
 * @note The trial matrix for 2 components, 3 phases is a 6-dimensional matrix with 720 entries to test per iteration of the minimization loop.
 * @note The current threshold to determine when a free energy has been minimized is when the trial step size is less than 1e-12.
 *
 * @param [in] temperature The temperature of the mixture.
 * @param [in] particlesATotal The total number of A particles in the mixture.
 * @param [in] particlesBTotal The total number of B particles in the mixture.
 * @param [in] volumeTotal The total volume occupied by the mixture.
 * @param [in,out] *particlesA1 Pointer to the number of A particles in phase 1.
 * @param [in,out] *particlesB1 Pointer to the number of B particles in phase 1.
 * @param [in,out] *volume1 Pointer to the volume occupied by phase 1.
 * @param [in,out] *particlesA2 Pointer to the number of A particles in phase 2.
 * @param [in,out] *particlesB2 Pointer to the number of B particles in phase 2.
 * @param [in,out] *volume2 Pointer to the volume occupied by phase 2.
 * @param [out] *freeEnergyReturn Pointer to store the total minimized free energy.
 * @param [out] *freeEnergy1Return Pointer to store the free energy of phase 1.
 * @param [out] *freeEnergy2Return Pointer to store the free energy of phase 2.
 * @param [out] *freeEnergy3Return Pointer to store the free energy of phase 3.
 *
 * @return A boolean flag indicating if the free energy minimization produced a phase change.
 */
int minimizeFreeEnergyTwoComponentsThreePhases(double temperature,
											   double particlesATotal, double particlesBTotal, double volumeTotal,
											   double *particlesA1, double *particlesB1, double *volume1,
											   double *particlesA2, double *particlesB2, double *volume2,
											   double *freeEnergyReturn, double *freeEnergy1Return, double *freeEnergy2Return, double *freeEnergy3Return) {
	int phaseChangeOccurred = 0;
	int freeEnergyMinimumIndexA1 = 1, freeEnergyMinimumIndexB1 = 1, freeEnergyMinimumIndexV1 = 1;
	int freeEnergyMinimumIndexA2 = 1, freeEnergyMinimumIndexB2 = 1, freeEnergyMinimumIndexV2 = 1; // initial min is center of 6-D trial matrix (1,1,1,1,1,1)
	double freeEnergy1 = 0.0, freeEnergy2 = 0.0, freeEnergy3 = 0.0;
	double invalidFreeEnergy = 500.* ncA;
	double freeEnergyThreshold = 1e-12; /// free energy changes smaller than this value terminate the minimization loop
	double freeEnergyComparisonThreshold = 1e-12; /// asymmetric phase diagrams had invalid phase separations along singularity line with 1e-13

	// This is the modification to each physical quantity (degree of freedom) for the trial free energy calculations
	// dummy defaults to 0.5 (half of step size used to calculate the theoretical phase diagram)
	double originalStepSize = minimizationStepFactor * minimizationParticlesStepSize;
	double trialStepSize = originalStepSize;
	double stepSizeThreshold = 1e-6;

	// 6-D array to hold trial free energies from which to choose the minimum value (particlesA1, particlesB1, volume1, particlesA2, particlesB2, volume2)
	// Each physical quantity - degree of freedom - has 3 trial values (reference value - step size, reference value, reference value + step size)
	struct FreeEnergy freeEnergyTrialArray[3][3][3][3][3][3];

	#ifdef DEBUG_MINIMIZATION_ON
	double debugValue = 0.57;
	double * debugVariable = &particlesATotal;
	FILE *debugData;
	if ( (*debugVariable < debugValue+freeEnergyThreshold) && (*debugVariable > debugValue-freeEnergyThreshold)) {
		char debugData_name[255];
		sprintf(debugData_name, "/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/debug/"
				"threePhaseDiagram-debugData_%f_%f_%f_%f_%f_%f_%f_%f_%f_%f_%f.dat",
				particlesATotal, particlesBTotal, *particlesA1, *particlesB1, *volume1, *particlesA2, *particlesB2, *volume2,
				(particlesATotal-*particlesA1-*particlesA2), (particlesBTotal-*particlesB1-*particlesB2), (volumeTotal-*volume1-*volume2));
		debugData = fopen(debugData_name, "w");
	}
	#endif

	// Initialize the minimum free energy to be that of the current test point on the phase diagram
	double freeEnergyMinimum = calculateFreeEnergyMinimumThreePhaseTrial(temperature, *particlesA1, *particlesB1, *volume1, *particlesA2, *particlesB2, *volume2,
			particlesATotal, particlesBTotal, volumeTotal, &freeEnergy1, &freeEnergy2, &freeEnergy3);

	#ifdef DEBUG_MINIMIZATION_ON
	if ((*debugVariable < debugValue+freeEnergyThreshold) && (*debugVariable > debugValue-freeEnergyThreshold)) {
		fprintf(debugData, "Starting free energy minimum = %.18f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n",
				freeEnergyMinimum, particlesATotal, particlesBTotal, volumeTotal, *particlesA1, *particlesB1, *volume1, *particlesA2, *particlesB2, *volume2,
				(particlesATotal - *particlesA1 - *particlesA2), (particlesBTotal - *particlesB1 - *particlesB2), (volumeTotal - *volume1 - *volume2));
		fprintf(debugData, "Total particles = %.15f\n\n", (particlesATotal+particlesBTotal));
	}
	#endif

	//
	// Minimize the free energy at a given test point on the phase diagram
	//

	while (trialStepSize > freeEnergyThreshold) {
		#pragma omp parallel for collapse(6) firstprivate(freeEnergy1, freeEnergy2, freeEnergy3) // OpenMP parallelizes this 6-level for loop
		for (int i = 0; i < 3; i++) { // particlesA1 trials
			for (int j = 0; j < 3; j++) { // particlesB1 trials
				for (int k = 0; k < 3; k++){ // volume1 trials
					for (int l = 0; l < 3; l++) { // particlesA2 trials
						for (int m = 0; m < 3; m++) { // particlesB2 trials
							for (int n = 0; n < 3; n++){ // volume2 trials
								freeEnergyTrialArray[i][j][k][l][m][n].freeEnergyTotal = calculateFreeEnergyMinimumThreePhaseTrial(temperature,
										*particlesA1+trialStepSize*(i-1), *particlesB1+trialStepSize*(j-1), *volume1+trialStepSize*(k-1),
										*particlesA2+trialStepSize*(l-1), *particlesB2+trialStepSize*(m-1), *volume2+trialStepSize*(n-1),
										particlesATotal, particlesBTotal, volumeTotal, &freeEnergy1, &freeEnergy2, &freeEnergy3);
								freeEnergyTrialArray[i][j][k][l][m][n].freeEnergy1 = freeEnergy1;
								freeEnergyTrialArray[i][j][k][l][m][n].freeEnergy2 = freeEnergy2;
								freeEnergyTrialArray[i][j][k][l][m][n].freeEnergy3 = freeEnergy3;

								#ifdef DEBUG_MINIMIZATION_ON
								if ((*debugVariable < debugValue+freeEnergyThreshold) && (*debugVariable > debugValue-freeEnergyThreshold)) {
									fprintf(debugData, "%.18f %.15f %.15f %.15f %f %f %.15f %.15f %.15f %f %f %.15f %.15f %.15f %f %f %i %i %i %i %i %i %.15f\n",
											freeEnergyTrialArray[i][j][k][l][m][n].freeEnergyTotal,
											*particlesA1+trialStepSize*(i-1), *particlesB1+trialStepSize*(j-1), *volume1+trialStepSize*(k-1),
											(*particlesA1+trialStepSize*(i-1)) / (*volume1+trialStepSize*(k-1)),
											(*particlesB1+trialStepSize*(j-1)) / (*volume1+trialStepSize*(k-1)),
											*particlesA2+trialStepSize*(l-1), *particlesB2+trialStepSize*(m-1), *volume2+trialStepSize*(n-1),
											(*particlesA2+trialStepSize*(l-1)) / (*volume2+trialStepSize*(n-1)),
											(*particlesB2+trialStepSize*(m-1)) / (*volume2+trialStepSize*(n-1)),
											(particlesATotal-(*particlesA1+trialStepSize*(i-1))-(*particlesA2+trialStepSize*(l-1))),
											(particlesBTotal-(*particlesB1+trialStepSize*(j-1))-(*particlesB2+trialStepSize*(m-1))),
											(volumeTotal-(*volume1+trialStepSize*(k-1))-(*volume2+trialStepSize*(n-1))),
											(particlesATotal-(*particlesA1+trialStepSize*(i-1))-(*particlesA2+trialStepSize*(l-1))) /
											(volumeTotal-(*volume1+trialStepSize*(k-1))-(*volume2+trialStepSize*(n-1))),
											(particlesBTotal-(*particlesB1+trialStepSize*(j-1))-(*particlesB2+trialStepSize*(m-1))) /
											(volumeTotal-(*volume1+trialStepSize*(k-1))-(*volume2+trialStepSize*(n-1))),
											i, j, k, l, m, n, trialStepSize);
								}
								#endif

							} // end for n
						} // end for m
					} // end for l
				} // end for k
			} // end for j
		} // end for i

		#ifdef DEBUG_MINIMIZATION_ON
		if ((*debugVariable < debugValue+freeEnergyThreshold) && (*debugVariable > debugValue-freeEnergyThreshold)) {
			fprintf(debugData, "\n");
		}
		#endif

		findIndexOfMinimumValue6DArray(freeEnergyTrialArray, &freeEnergyMinimumIndexA1, &freeEnergyMinimumIndexB1, &freeEnergyMinimumIndexV1,
									   &freeEnergyMinimumIndexA2, &freeEnergyMinimumIndexB2, &freeEnergyMinimumIndexV2);

		#ifdef DEBUG_MINIMIZATION_ON
		if ((*debugVariable < debugValue+freeEnergyThreshold) && (*debugVariable > debugValue-freeEnergyThreshold)) {
			fprintf(debugData, "Current free energy minimum = %.18f\tCompared to... %.18f\t%i %i %i %i %i %i\n\n",
					freeEnergyMinimum,
					freeEnergyTrialArray[freeEnergyMinimumIndexA1][freeEnergyMinimumIndexB1][freeEnergyMinimumIndexV1]
					                    [freeEnergyMinimumIndexA2][freeEnergyMinimumIndexB2][freeEnergyMinimumIndexV2].freeEnergyTotal,
					freeEnergyMinimumIndexA1, freeEnergyMinimumIndexB1, freeEnergyMinimumIndexV1,
					freeEnergyMinimumIndexA2, freeEnergyMinimumIndexB2, freeEnergyMinimumIndexV2);
		}
		#endif

		if (freeEnergyTrialArray[freeEnergyMinimumIndexA1][freeEnergyMinimumIndexB1][freeEnergyMinimumIndexV1]
		                        [freeEnergyMinimumIndexA2][freeEnergyMinimumIndexB2][freeEnergyMinimumIndexV2].freeEnergyTotal > invalidFreeEnergy) {
			trialStepSize = freeEnergyThreshold / 10; // dope the step size when the initialized free energy is not valid; escape the minimization loop now
		}
		else if ( (freeEnergyTrialArray[freeEnergyMinimumIndexA1][freeEnergyMinimumIndexB1][freeEnergyMinimumIndexV1]
		                               [freeEnergyMinimumIndexA2][freeEnergyMinimumIndexB2][freeEnergyMinimumIndexV2].freeEnergyTotal +
		                               freeEnergyComparisonThreshold) < freeEnergyMinimum ) {

			// If a new minimum free energy is found, indicate a phase change has occurred and save the minimum information
			phaseChangeOccurred = 1;

			#ifdef DEBUG_MINIMIZATION_ON
			if ((*debugVariable < debugValue+freeEnergyThreshold) && (*debugVariable > debugValue-freeEnergyThreshold)) {
				fprintf(debugData, "Phase change occurred... replacing %.18f with %.18f\t%i %i %i %i %i %i\n\n",
						freeEnergyMinimum,
						freeEnergyTrialArray[freeEnergyMinimumIndexA1][freeEnergyMinimumIndexB1][freeEnergyMinimumIndexV1]
						                    [freeEnergyMinimumIndexA2][freeEnergyMinimumIndexB2][freeEnergyMinimumIndexV2].freeEnergyTotal,
						freeEnergyMinimumIndexA1, freeEnergyMinimumIndexB1, freeEnergyMinimumIndexV1,
						freeEnergyMinimumIndexA2, freeEnergyMinimumIndexB2, freeEnergyMinimumIndexV2);
			}
			#endif

			*particlesA1 += (freeEnergyMinimumIndexA1-1) * trialStepSize;
			*particlesB1 += (freeEnergyMinimumIndexB1-1) * trialStepSize;
			*volume1 += (freeEnergyMinimumIndexV1-1) * trialStepSize;

			*particlesA2 += (freeEnergyMinimumIndexA2-1) * trialStepSize;
			*particlesB2 += (freeEnergyMinimumIndexB2-1) * trialStepSize;
			*volume2 += (freeEnergyMinimumIndexV2-1) * trialStepSize;

			freeEnergyMinimum = freeEnergyTrialArray[freeEnergyMinimumIndexA1][freeEnergyMinimumIndexB1][freeEnergyMinimumIndexV1]
			                                        [freeEnergyMinimumIndexA2][freeEnergyMinimumIndexB2][freeEnergyMinimumIndexV2].freeEnergyTotal;

			// Expand the step size when a good minimization is found to help speed up the process
			if (!isEqual(trialStepSize, originalStepSize, stepSizeThreshold)) trialStepSize *= 2.0;
		}
		else trialStepSize /= 2; // If there is not a new free energy minimum (previous step is still the min), reduce the step size for the next iteration
	} // end while loop; free energy minimization complete

	#ifdef DEBUG_MINIMIZATION_ON
	if ((*debugVariable < debugValue+freeEnergyThreshold) && (*debugVariable > debugValue-freeEnergyThreshold)) fclose(debugData);
	#endif

	//
	// Return the free energy minimum and if a phase change occurred
	//

	*freeEnergyReturn = freeEnergyMinimum;
	*freeEnergy1Return = freeEnergyTrialArray[freeEnergyMinimumIndexA1][freeEnergyMinimumIndexB1][freeEnergyMinimumIndexV1]
	                                         [freeEnergyMinimumIndexA2][freeEnergyMinimumIndexB2][freeEnergyMinimumIndexV2].freeEnergy1;
	*freeEnergy2Return = freeEnergyTrialArray[freeEnergyMinimumIndexA1][freeEnergyMinimumIndexB1][freeEnergyMinimumIndexV1]
	                                         [freeEnergyMinimumIndexA2][freeEnergyMinimumIndexB2][freeEnergyMinimumIndexV2].freeEnergy2;
	*freeEnergy3Return = freeEnergyTrialArray[freeEnergyMinimumIndexA1][freeEnergyMinimumIndexB1][freeEnergyMinimumIndexV1]
	                                         [freeEnergyMinimumIndexA2][freeEnergyMinimumIndexB2][freeEnergyMinimumIndexV2].freeEnergy3;
	return phaseChangeOccurred;
} // end function minimizeFreeEnergyTwoComponentsThreePhases()


/**
 * @brief Function _detectTwoOrThreePhaseRegion_ determines if a particle, volume configuration for a free energy is in a 2- or 3-phase region and
 * logs the configuration to an output file.
 *
 * This function accepts a particle and volume configuration associated with a free energy, calculates the densities of the A and B components,
 * and compares the densities among all 3 phases to determine if a phase separation occurred in the 2-phase or 3-phase region of a theoretical
 * phase diagram.  If the configuration phase separated to a region where it is still unstable, no information is logged and the return value
 * indicates which phase is still unstable (-1,-2,-3).  If the configuration separated to stable phases, the number of phases in the configuration (2 or 3)
 * is prepended to the particle, volume free energy configuration and logged to an output file.
 *
 * @param [in] phaseDetectionThreshold Threshold to use for determining if one phase's density is equal to another.
 * @param [out] densities_file Pointer to the file to which the free energy configuration is written.
 * @param [in,out] particlesATotal The total number of A particles in the mixture.
 * @param [in,out] particlesBTotal The total number of B particles in the mixture.
 * @param [in,out] volumeTotal The total volume occupied by the mixture.
 * @param [in,out] particlesA1 The number of A particles in phase 1.
 * @param [in,out] particlesB1 The number of B particles in phase 1.
 * @param [in,out] volume1 The volume occupied by phase 1.
 * @param [in,out] particlesA2 The number of A particles in phase 2.
 * @param [in,out] particlesB2 The number of B particles in phase 2.
 * @param [in,out] volume2 The volume occupied by phase 2.
 * @param [out] freeEnergy The total minimized free energy.
 * @param [out] freeEnergy1 The free energy of phase 1.
 * @param [out] freeEnergy2 The free energy of phase 2.
 * @param [out] freeEnergy3 The free energy of phase 3.
 *
 * @return Integer indicating if the configuration phase separated to an unstable region (-1,-2,-3) or is in a stable 2-phase (2) or 3-phase (3) region.
 * (invalid returns 0)
 */
int detectTwoOrThreePhaseRegion(double phaseDetectionThreshold, FILE *densities_file,
								double particlesATotal, double particlesBTotal, double volumeTotal,
								double particlesA1, double particlesB1, double volume1,
								double particlesA2, double particlesB2, double volume2,
								double freeEnergy, double freeEnergy1, double freeEnergy2, double freeEnergy3) {
	int twoOrThreePhases = 0;
	double zeroThreshold = 1e-4; // particles/volumina below this value from the free energy minimization are considered to be zero
	double phaseSeparationThreshold = 5.0 * minimizationParticlesStepSize; // used to filter out false positives for 3-phase behavior
	double ev1 = 0.0, ev2 = 0.0; // dummy variables for determining stability

	// Calculate the densities from the minimization
	double particlesA3 = particlesATotal-particlesA1-particlesA2, particlesB3 = particlesBTotal-particlesB1-particlesB2, volume3 = volumeTotal-volume1-volume2;
	double rhoA1 = particlesA1 / volume1, rhoA2 = particlesA2 / volume2, rhoA3 = particlesA3 / volume3;
	double rhoB1 = particlesB1 / volume1, rhoB2 = particlesB2 / volume2, rhoB3 = particlesB3 / volume3;

	if ( (isEqual(particlesA1, 0.0, zeroThreshold) && isEqual(particlesB1, 0.0, zeroThreshold) && isEqual(volume1, 0.0, zeroThreshold)) ||
			(isEqual(particlesA2, 0.0, zeroThreshold) && isEqual(particlesB2, 0.0, zeroThreshold) && isEqual(volume2, 0.0, zeroThreshold)) ||
			(isEqual(particlesA3, 0.0, zeroThreshold) && isEqual(particlesB3, 0.0, zeroThreshold) && isEqual(volume3, 0.0, zeroThreshold)) ) {
		twoOrThreePhases = 2;
	}
	else if (calculateMinimizationPath(rhoA1, rhoB1, &ev1, &ev2)) twoOrThreePhases = -1; // flag ph1 for more minimizing; try to force the true 3-phase behavior
	else if (calculateMinimizationPath(rhoA2, rhoB2, &ev1, &ev2)) twoOrThreePhases = -2; // flag ph2 for more minimizing; try to force the true 3-phase behavior
	else if (calculateMinimizationPath(rhoA3, rhoB3, &ev1, &ev2)) twoOrThreePhases = -3; // flag ph3 for more minimizing; try to force the true 3-phase behavior
	else if ( (isEqual(rhoA1, rhoA2, phaseDetectionThreshold) && isEqual(rhoB1, rhoB2, phaseDetectionThreshold)) ||	// two-phase region 1=2
			(isEqual(rhoA1, rhoA3, phaseDetectionThreshold) && isEqual(rhoB1, rhoB3, phaseDetectionThreshold)) ||	// two-phase region 1=3
			(isEqual(rhoA2, rhoA3, phaseDetectionThreshold) && isEqual(rhoB2, rhoB3, phaseDetectionThreshold)) ) {	// two-phase region 2=3
		twoOrThreePhases = 2;
	}
	else { // three-phase region 1!=2, 2!=3, 3!=1
		double phase1Separation = sqrt(pow(fabs(particlesATotal-rhoA1),2) + pow(fabs(particlesBTotal-rhoB1),2));
		double phase2Separation = sqrt(pow(fabs(particlesATotal-rhoA2),2) + pow(fabs(particlesBTotal-rhoB2),2));
		double phase3Separation = sqrt(pow(fabs(particlesATotal-rhoA3),2) + pow(fabs(particlesBTotal-rhoB3),2));

		// phase1 ? (phase2||phase3 ? 1 : 0) : (phase2&&phase3 ? 1 : 0)
		if ( (phase1Separation < phaseSeparationThreshold) ?
				((phase2Separation < phaseSeparationThreshold || phase3Separation < phaseSeparationThreshold) ? 1 : 0) :
				((phase2Separation < phaseSeparationThreshold && phase3Separation < phaseSeparationThreshold) ? 1 : 0) ) {
			twoOrThreePhases = 2;
		}
		else twoOrThreePhases = 3;
	}

	if (twoOrThreePhases == 2 || twoOrThreePhases == 3) {
		fprintf(densities_file, "%i %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %f %f %f %f\n",
				twoOrThreePhases,
				particlesATotal, particlesBTotal, volumeTotal,
				particlesA1, particlesB1, volume1, (particlesA1/volume1), (particlesB1/volume1),
				particlesA2, particlesB2, volume2, (particlesA2/volume2), (particlesB2/volume2),
				(particlesATotal-particlesA1-particlesA2), (particlesBTotal-particlesB1-particlesB2), (volumeTotal-volume1-volume2),
				(particlesATotal-particlesA1-particlesA2) / (volumeTotal-volume1-volume2),
				(particlesBTotal-particlesB1-particlesB2) / (volumeTotal-volume1-volume2),
				freeEnergy, freeEnergy1, freeEnergy2, freeEnergy3);
	}

	return twoOrThreePhases;
} // end function detectTwoOrThreePhases()


/**
 * @brief Function _insideThreePhaseRegion_ determines if a phase diagram test point is inside of a defined 3-phase region.  The 3-phase region may include points that exhibit
 * metastable 2-phase behavior.
 *
 * It is assumed that the 3-phase region of a phase diagram is given by a triangular shape.  The coordinates of the 3-phase region are
 * given by the global array _threePhaseRegion[]_.  The global flag _setThreePhaseRegion_ indicates if the triangle coordinates have already
 * been set of they need to be read from a file.
 *
 * This function takes a point and a triangle and determines if the point is bounded by the triangle.
 * The calculation is in terms of barycentric coordinates.
 * - particlesATotal = x
 * - particlesBTotal = y
 * - threePhaseRegion coordinate indices: 0 = x1, 1 = y1, 2 = x2, 3 = y2, 4 = x3, 5 = y3
 *
 * @param [in] particlesATotal Total number of A particles in the given phase diagram test point.
 * @param [in] particlesBTotal Total number of B particles in the given phase diagram test point.
 *
 * @return Boolean value indicating if the point is inside (1) or outside (0) the triangle.
 */
int insideThreePhaseRegion(double particlesATotal, double particlesBTotal) {
	if (setThreePhaseRegion) readThreePhaseRegionCoordinates(); // make sure there are coordinates for the 3-phase region

	double denominator = ( (threePhaseRegion[3]-threePhaseRegion[5])*(threePhaseRegion[0]-threePhaseRegion[4]) +
			        (threePhaseRegion[4]-threePhaseRegion[2])*(threePhaseRegion[1]-threePhaseRegion[5]) );
	double a = ( (threePhaseRegion[3]-threePhaseRegion[5])*(particlesATotal-threePhaseRegion[4]) +
		  (threePhaseRegion[4]-threePhaseRegion[2])*(particlesBTotal-threePhaseRegion[5]) ) / denominator;
	double b = ( (threePhaseRegion[5]-threePhaseRegion[1])*(particlesATotal-threePhaseRegion[4]) +
		  (threePhaseRegion[0]-threePhaseRegion[4])*(particlesBTotal-threePhaseRegion[5]) ) / denominator;
	double c = 1 - a - b;

	return 0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1;
} // end function insideThreePhaseRegion()


/**
 * @brief Function _unconditionallyUnstableThreePhase_ determines if a phase diagram test point is inside of a region where instability results in unconditional 3-phase behavior.
 *
 * It is assumed that the 3-phase region of a phase diagram is given by a triangular shape.  This "outer" triangle may have points that exhibit metastable 2-phase behavior; however,
 * another "inner" triangle exists where instabilities result in unconditional 3-phase behavior.  This unconditionally unstable region is defined by creating three tie lines that
 * connect the very last points of the binodal curves that define the binary liquid, component A liquid-vapor, and component B liquid-vapor regions.  The region inscribed by these
 * three lines defines the unconditionally unstable 3-phase region.
 *
 * @note A test point that falls on a line is assumed to be not inside of the unconditionally unstable 3-phase region.
 *
 * @param [in] particlesATotal Total number of A particles in the given phase diagram test point.
 * @param [in] particlesBTotal Total number of B particles in the given phase diagram test point.
 *
 * @return Boolean value indicating if the point is inside (1) or outside (0) the triangle.
 */
int unconditionallyUnstableThreePhase(double pointA, double pointB) {
	double position = 100.0;
	FILE *unstable_region_coords = openFile("unconditionally-unstable-region", ".dat", "r");
	if (unstable_region_coords) {
		struct BinodalPoint binaryLiquidLine1;
		struct BinodalPoint binaryLiquidLine2;
		struct BinodalPoint liquidVaporALine1;
		struct BinodalPoint liquidVaporALine2;
		struct BinodalPoint liquidVaporBLine1;
		struct BinodalPoint liquidVaporBLine2;
		fscanf(unstable_region_coords, "%8lf %8lf", &binaryLiquidLine1.rhoA, &binaryLiquidLine1.rhoB);
		fscanf(unstable_region_coords, "%8lf %8lf", &binaryLiquidLine2.rhoA, &binaryLiquidLine2.rhoB);
		fscanf(unstable_region_coords, "%8lf %8lf", &liquidVaporALine1.rhoA, &liquidVaporALine1.rhoB);
		fscanf(unstable_region_coords, "%8lf %8lf", &liquidVaporALine2.rhoA, &liquidVaporALine2.rhoB);
		fscanf(unstable_region_coords, "%8lf %8lf", &liquidVaporBLine1.rhoA, &liquidVaporBLine1.rhoB);
		fscanf(unstable_region_coords, "%8lf %8lf", &liquidVaporBLine2.rhoA, &liquidVaporBLine2.rhoB);
		fclose(unstable_region_coords);

		// (y-y1)*(x2-x1) - (x-x1)*(y2-y1)
		// A-rich point is always point 1; B-rich point is always point 2... "above" the line is + position; "below" the line is - position
		position = (pointB-binaryLiquidLine1.rhoB)*(binaryLiquidLine2.rhoA-binaryLiquidLine1.rhoA)
								- (pointA-binaryLiquidLine1.rhoA)*(binaryLiquidLine2.rhoB-binaryLiquidLine1.rhoB);
		if (position > 0.0) {
			position = (pointB-liquidVaporALine1.rhoB)*(liquidVaporALine2.rhoA-liquidVaporALine1.rhoA)
													- (pointA-liquidVaporALine1.rhoA)*(liquidVaporALine2.rhoB-liquidVaporALine1.rhoB);
			if (position < 0.0) {
				position = (pointB-liquidVaporBLine1.rhoB)*(liquidVaporBLine2.rhoA-liquidVaporBLine1.rhoA)
														- (pointA-liquidVaporBLine1.rhoA)*(liquidVaporBLine2.rhoB-liquidVaporBLine1.rhoB);
				if (position < 0.0) return 1; // return inside triangle only if below binary liquid line and above both liquid-vapor lines
			}
		}
	}
	else printf("Binodals for this phase diagram must be constructed before the unconditionally unstable 3-phase region can be defined.\n\n");
	return 0; // always return outside if the line checks above don't succeed
} // end function unconditionallyUnstableThreePhase()


/**
 * @brief Function _defineThreePhaseRegion_ uses coordinates from 3-phase free energy configurations to define the triangular bounds of the region.
 *
 * This function sets the coordinates that define a triangular 3-phase region of a phase diagram.  The region is determined by maximizing the area
 * of the triangle that bounds the region.  The global variable _maxArea_ tracks the maximum area, and the coordinates of the triangle are written
 * to the global array _threePhaseRegion[]_.
 *
 * @param [in] x1 X-coordinate of the triangle's 1st vertex.
 * @param [in] y1 Y-coordinate of the triangle's 1st vertex.
 * @param [in] x2 X-coordinate of the triangle's 2nd vertex.
 * @param [in] y2 Y-coordinate of the triangle's 2nd vertex.
 * @param [in] x3 X-coordinate of the triangle's 3rd vertex.
 * @param [in] y3 Y-coordinate of the triangle's 3rd vertex.
 */
int defineThreePhaseRegion(double x1, double y1, double x2, double y2, double x3, double y3) {
	double area = fabs(0.5 * (x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)));
	if (area > maxArea) {
		threePhaseRegion[0] = x1;
		threePhaseRegion[1] = y1;
		threePhaseRegion[2] = x2;
		threePhaseRegion[3] = y2;
		threePhaseRegion[4] = x3;
		threePhaseRegion[5] = y3;
		maxArea = area;
		threePhaseRegionExists = 1;
		return 1;
	}
	return 0;
} // end function defineThreePhaseRegion()


/**
 * @brief Function _tieLineNeeded_ determines if a theoretical tie line is to be drawn for a given phase diagram test point.
 *
 * This function compares a phase diagram test point to a list of points for which theoretical tie lines are to be drawn.  If the test point
 * is equal to a point in the list, a tie line is needed.
 *
 * @param [in] particlesA Number of A particles in the test point.
 * @param [in] particlesB Number of B particles in the test point.
 *
 * @return Boolean indicating if a tie line is to be drawn (1) or not (0).
 */
int tieLineNeeded(double particlesA, double particlesB) {

	// Static arrays that determine the A,B particle phase diagram test point for which theoretical tie lines are to be drawn
	static double A[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4};
	static double B[] = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
			1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.1, 1.2, 1.3, 1.4};

	int particlePairs = *(&A + 1) - A; // equivalent to sizeof(A)/sizeof(double);
	int tieLineNeeded = 0;
	for (int i = 0; i < particlePairs; i++) {
		double threshold = 1e-6;
		if (isEqual(A[i], particlesA, threshold)) {
			if (isEqual(B[i], particlesB, threshold)) {
				tieLineNeeded = 1;
				break;
			}
		}
	}
	return tieLineNeeded;
} // end function tieLineNeeded()


/**
 * @brief Function _backSideOfPhaseDiagram_ determines if a given phase diagram test point is in the binary liquid region of the phase diagram.
 *
 * This function determines if a test point is on the "back" of the phase diagram (binary liquid region) by comparing the point to a line
 * which defines the "front" and "back" parts of the phase diagram.  The line is is constructed to have a negative slope and is the diagonal
 * connecting the first and last points in the 3-phase data.  If the phase diagram has no 3-phase region, the line connects either the largest
 * A- and B-component axis intercepts (when at least one exists) or  the van der Waals singluarity points for each component.  The global flags
 *
 * When a 3-phase region is not well-behaved and the line is not entirely accurate, this function will also scan for 3-phase behavior by
 * scanning for 3-phase behavior adjacent to the test point.  If all 3-phase behavior is "below" the test point, the test point is in the binary
 * liquid region.  Otherwise, it is not.
 *
 * @param [in] pointA A-coordinate of the test point.
 * @param [in] pointB B-coordinate of the test point.
 *
 * @return Integer indicating if the test point is on the front (+1) or back (-1) of the phase diagram (0 is on the line).
 */
int backSideOfPhaseDiagram(double pointA, double pointB) {
	static double lineA1 = 0, lineB1 = 0, lineA2 = 0, lineB2 = 0;
	static double minB = 10.0; // smallest B-coordinate of a point with 3-phase behavior
	int readEOF = 0;
	double position = 100.0;
	double particlesATotal = 0.0, particlesBTotal = 0.0;
	double rhoA1 = 0, rhoA2 = 0, rhoA3 = 0;
	double rhoB1 = 0, rhoB2 = 0, rhoB3 = 0;
	double zeroThreshold = 1e-14;

	if (!setLineAPoint && !setLineBPoint) { // if the line is already set...
		position = (pointB-lineB1)*(lineA2-lineA1) - (pointA-lineA1)*(lineB2-lineB1); // (y-y1)*(x2-x1) - (x-x1)*(y2-y1)
		if (position < 0.0) {
			if (pointB < minB && threePhaseRegionExists) return 1; // helps keep spinodal of liquid/vapor regions from getting clipped by the line
			else return -1; // - slope on phase diagram, so this is the true/positive result
		}
		else if (position >= 0.0 && !insideThreePhaseRegion(pointA, pointB)) return 1; // skip column scanning below if not metastable point in 3-phase region
		// else if the line exists but may not be entirely accurate; scan for 3-phase behavior in the same column as the test point
		// else the line doesn't exist; build it
	}

	if (threePhaseRegionExists) { // use 3-phase data to either build the line or scan point-by-point to see when a metastable point is really on the back side
		FILE *three_phase_coords = openFile("threePhaseDiagram-densities-threePhases", ".dat", "r");
		readEOF = fscanf(three_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
				&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);

		if (setLineBPoint) { // set the first line point; first line in the file is one end point
			lineA2 = particlesATotal;
			lineB2 = particlesBTotal;
			setLineBPoint = 0;
		}

		// Scan 3-phase data to see if there is 3-phase behavior in the same column as the test point
		// TODO: optimize this just a bit by detecting when moved past the same column; if line is already set, break loop early
		double phaseDetectionThreshold = 1e-4;
		while (readEOF != EOF) {
			if (rhoB1 < minB) minB = rhoB1; // set a "global" lower cut-off for 3-phase behavior
			if (isEqual(pointA, particlesATotal, phaseDetectionThreshold)) { // if 3-phase point in the same column as test point...
				if (pointB > particlesBTotal) position = -1; // declare test point on back side if 3-phase behavior is below the test point
				else position = 1; // toggle to front of phase diagram if another 3-phase point is above the test point
			}
			readEOF = fscanf(three_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
					&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
		}

		if (setLineAPoint) { // set the second line point; last line in the file is other end point
			lineA1 = particlesATotal;
			lineB1 = particlesBTotal;
			setLineAPoint = 0;
		}
		fclose(three_phase_coords);
	}
	else { // if no 3-phase region, set the line to be either the axis intercepts or the singularity line
		FILE *two_phase_coords = openFile("threePhaseDiagram-densities-twoPhases", ".dat", "r");
		readEOF = fscanf(two_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
				&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
		while (readEOF != EOF && (setLineAPoint || setLineBPoint)) { // find axis intercepts (if any)
			if (isEqual(rhoA2, 0.0, zeroThreshold) && rhoB2 > lineB2) {
				lineB2 = rhoB2; // keep lineA2 = 0 from initialization
				setLineBPoint = 0;
			}
			if (isEqual(rhoB1, 0.0, zeroThreshold) && rhoA1 > lineA1) {
				lineA1 = rhoA1; // keep lineB1 = 0 from initialization
				setLineAPoint = 0;
			}
			readEOF = fscanf(two_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
					&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
		}
		fclose(two_phase_coords);

		if (setLineAPoint && setLineBPoint) { // use VDW singularities if no axis intercepts
			lineA1 = 1.0 / tcA; // keep lineB1 = 0 from initialization
			lineB2 = 1.0 / tcB; // keep lineA2 = 0 from initialization
			setLineAPoint = 0;
			setLineBPoint = 0;
		}
		if (setLineBPoint) { // if only 1 intercept is found, mirror the other
			lineB2 = lineA1;
			setLineBPoint = 0;
		}
		if (setLineAPoint) {
			lineA1 = lineB2;
			setLineAPoint = 0;
		}
	}

	// Determine position of the test point relative to the line and return the result
	if (isEqual(position, 100.0, zeroThreshold)) position = (pointB-lineB1)*(lineA2-lineA1) - (pointA-lineA1)*(lineB2-lineB1); // (y-y1)*(x2-x1) - (x-x1)*(y2-y1)
	if (position > 0.0) return 1;
	else if (position < 0.0) return -1; // - slope on phase diagram, so this is the true/positive result
	else return 0; // on the line
} // end function backSideOfPhaseDiagram()


/**
 * @brief Function _sortMetastableMinimizationResults_ evaluates points within the 3-phase region of a phase diagram and identifies if metastable behavior (either 2- or 3-phase)
 * or unconditionally unstable 3-phase behavior is expected.
 *
 * This function reads the log files for both the 2-phase and 3-phase results from the free energy minimization and sorts all points into one of three groups defined by the
 * phase separation behavior that is expected:
 * 1. unconditionally unstable phase separation into 3-phases
 * 2. metastable behavior that may result in 3-phases forming
 * 3. unconditionally unstable phase separation into 2-phases
 *
 * @todo Shouldn't happen, but detect errors where the minimization algorithm returns a 2-phase result in the unconditionally unstable 3-phase region OR a 3-phase result outside
 * of the full 3-phase region definition.
 */
void sortMetastableMinimizationResults() {
	if (setThreePhaseRegion) readThreePhaseRegionCoordinates(); // need 3-phase coords to define metastable points

	if (threePhaseRegionExists) {
		int readEOF = 0;
		double particlesA = 0.0, particlesB = 0.0;
		double rhoA1 = 0.0, rhoB1 = 0.0;
		double rhoA2 = 0.0, rhoB2 = 0.0;
		double rhoA3 = 0.0, rhoB3 = 0.0;
		FILE *phase_results = openFile("threePhaseDiagram-densities-twoPhases", ".dat", "r");
		FILE *metastable_results = openFile("minimization-metastable-densities", ".dat", "w");
		FILE *unconditional_2_results = openFile("minimization-unconditionally-2-phase", ".dat", "w");
		FILE *unconditional_3_results = openFile("minimization-unconditionally-3-phase", ".dat", "w");

		// Categorize minimization results with 2-phase behavior
		readEOF = fscanf(phase_results, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f", &particlesA, &particlesB, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
		while (readEOF != EOF) {
			//			if (unconditionallyUnstableThreePhase(particlesA, particlesB)) { // should never be any here...
			//				fprintf(unconditional_3_results,"%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", particlesA, particlesB, rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);
			//			}
			if (insideThreePhaseRegion(particlesA, particlesB)) {
				fprintf(metastable_results,"%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", particlesA, particlesB, rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);
			}
			else {
				fprintf(unconditional_2_results,"%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", particlesA, particlesB, rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);
			}

			readEOF = fscanf(phase_results, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f", &particlesA, &particlesB, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
		} // end while

		fclose(phase_results);
		phase_results = openFile("threePhaseDiagram-densities-threePhases", ".dat", "r");

		// Categorize minimization results with 3-phase behavior
		readEOF = fscanf(phase_results, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f", &particlesA, &particlesB, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
		while (readEOF != EOF) {
			if (unconditionallyUnstableThreePhase(particlesA, particlesB)) {
				fprintf(unconditional_3_results,"%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", particlesA, particlesB, rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);
			}
			else if (insideThreePhaseRegion(particlesA, particlesB)) {
				fprintf(metastable_results,"%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", particlesA, particlesB, rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);
			}
			else { // should never be any here...
				fprintf(unconditional_2_results,"%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", particlesA, particlesB, rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);
			}

			readEOF = fscanf(phase_results, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f", &particlesA, &particlesB, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
		} // end while

		fclose(phase_results);
		fclose(metastable_results);
		fclose(unconditional_2_results);
		fclose(unconditional_3_results);
		printf("Metastable minimization results sorted.\n\n");
	} // end if (threePhaseRegionExists)
	else printf("\nCannot classify metastable behavior without a defined 3-phase region.\n\n");
} // end function sortMetastableMinimizationResults()


/**
 * @brief Function _generatePhaseDiagramBinodalLines_ generates files that define the binodal lines of a phase diagram.
 *
 * This function reads the log file for the 2-phase points of the phase diagram and builds phase diagram binodals from the densities in phases
 * 1 and 2.  Binodals are divided into six segments, and densities are sorted and saved into six arrays (1 array per binodal segment).  The sorting
 * logic determines if a density is in the A-rich liquid-vapor region, B-rich liquid-vapor, or A/B liquid-liquid region of the phase diagram.
 * The arrays are then output to binodal log files suitable for plotting with gnuplot.
 */
void generatePhaseDiagramBinodalLines() {
	int readEOF = 0;
	int backOrFront = 0;
	int binodalArraySize = 0;
	int binodalInsertIndex = 0;
	int binodalAVaporTail = -1, binodalALiquidTail = -1, binodalARichTail = -1;
	int binodalBVaporTail = -1, binodalBLiquidTail = -1, binodalBRichTail = -1;
	int binodalAVaporBreak = -1, binodalALiquidBreak = -1, binodalARichBreak = -1;
	int binodalBVaporBreak = -1, binodalBLiquidBreak = -1, binodalBRichBreak = -1;
	double particlesATotal = 0.0, particlesBTotal = 0.0;
	double rhoA1 = 0.0, rhoA2 = 0.0, rhoA3 = 0.0;
	double rhoB1 = 0.0, rhoB2 = 0.0, rhoB3 = 0.0;
	double phaseDetectionThreshold = 1e-4;
	double zeroThreshold = 1e-14;
	double metastableFilterThreshold = 0.01;
	double binodalBreakStep = 2.0 * minimizationParticlesStepSize;

	printf("Building binodal lines...\n");
	binodalArraySize = 0.5 * pow(((bA<bB ? 1/bA : 1/bB)/minimizationParticlesStepSize), 2.0);  // full triangle of possible points that phase separate

	struct BinodalPoint binodalPoint1, binodalPoint2;
	struct BinodalPoint * binodalAVaporArray = calloc(binodalArraySize, sizeof(binodalPoint1));
	struct BinodalPoint * binodalALiquidArray = calloc(binodalArraySize, sizeof(binodalPoint1));
	struct BinodalPoint * binodalARichArray = calloc(binodalArraySize, sizeof(binodalPoint1));
	struct BinodalPoint * binodalBVaporArray = calloc(binodalArraySize, sizeof(binodalPoint1));
	struct BinodalPoint * binodalBLiquidArray = calloc(binodalArraySize, sizeof(binodalPoint1));
	struct BinodalPoint * binodalBRichArray = calloc(binodalArraySize, sizeof(binodalPoint1));

	if (setThreePhaseRegion) readThreePhaseRegionCoordinates();

	FILE *two_phase_coords = openFile("threePhaseDiagram-densities-twoPhases", ".dat", "r");
	readEOF = fscanf(two_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
			&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
	while (readEOF != EOF) {
		if (!isnan(rhoA1) && !isnan(rhoB1) && !isnan(rhoA2) && !isnan(rhoB2) && !isnan(rhoA3) && !isnan(rhoB3)) {
			checkVaporPhase(&rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3, particlesATotal, particlesBTotal);
			binodalPoint1.rhoA = rhoA1;
			binodalPoint1.rhoB = rhoB1;
			binodalPoint2.rhoA = rhoA2;
			binodalPoint2.rhoB = rhoB2;

			// Insert the binodal points in pairs depending on where in the phase diagram the parent/spinodal region point is located
			binodalInsertIndex = 0;
			backOrFront = backSideOfPhaseDiagram(particlesATotal, particlesBTotal);
			if (backOrFront == -1 || backOrFront == 0 ) { // binary liquid region
				for (int i = 0; i <= binodalARichTail; i++) {
					if (rhoA1 >= binodalARichArray[i].rhoA) {
						binodalInsertIndex = i + 1;
					}
				}
				arrayInsertBinodalPoint(binodalARichArray, binodalPoint1, binodalInsertIndex, &binodalARichTail);
				arrayInsertBinodalPoint(binodalBRichArray, binodalPoint2, binodalInsertIndex, &binodalBRichTail);
			}
			else if ( (particlesATotal > particlesBTotal) && (!isEqual(rhoA1, rhoA2, phaseDetectionThreshold)) ) { // A-liquid/vapor region
				for (int i = 0; i <= binodalAVaporTail; i++) {
					if (sortABinodalByIncreasingB) {
						if ((rhoB2 >= binodalAVaporArray[i].rhoB && vdwInteractionFactor <= 1.0) ||
								(rhoB1 >= binodalALiquidArray[i].rhoB && vdwInteractionFactor > 1.0)) binodalInsertIndex = i + 1;
					}
					else {
						if (rhoA2 >= binodalAVaporArray[i].rhoA) binodalInsertIndex = i + 1;
					}
				}

				if (rhoA1 > rhoB2 + metastableFilterThreshold) { // only add if it is not a metastable point where the phases get flipped
					arrayInsertBinodalPoint(binodalAVaporArray, binodalPoint2, binodalInsertIndex, &binodalAVaporTail);
					arrayInsertBinodalPoint(binodalALiquidArray, binodalPoint1, binodalInsertIndex, &binodalALiquidTail);
				}
			}
			else if ( (particlesBTotal > particlesATotal) && (!isEqual(rhoB1, rhoB2, phaseDetectionThreshold)) ) { // B-liquid vapor region
				for (int i = 0; i <= binodalBVaporTail; i++) {
					if (sortBBinodalByIncreasingA) {
						if ((rhoA1 >= binodalBVaporArray[i].rhoA && vdwInteractionFactor <= 1.0) ||
								(rhoA2 >= binodalBLiquidArray[i].rhoA && vdwInteractionFactor > 1.0)) binodalInsertIndex = i + 1;
					}
					else {
						if (rhoB1 >= binodalBVaporArray[i].rhoB) binodalInsertIndex = i + 1;
					}
				}

				if (rhoB2 > rhoA1 + metastableFilterThreshold) { // only add if it is not a metastable point where the phases get flipped
					arrayInsertBinodalPoint(binodalBVaporArray, binodalPoint1, binodalInsertIndex, &binodalBVaporTail);
					arrayInsertBinodalPoint(binodalBLiquidArray, binodalPoint2, binodalInsertIndex, &binodalBLiquidTail);
				}
			}
		}
		readEOF = fscanf(two_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
				&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
	} // end while
	fclose(two_phase_coords);

	// Scan binodal data arrays to see if there are large jumps (discontinuities in the binodal lines)... will only detect the first such one
	binodalAVaporBreak = binodalAVaporTail;
	binodalALiquidBreak = binodalALiquidTail;
	binodalARichBreak = binodalARichTail;
	binodalBVaporBreak = binodalBVaporTail;
	binodalBLiquidBreak = binodalBLiquidTail;
	binodalBRichBreak = binodalBRichTail;
	for (int i = 1; i < binodalAVaporTail; i++) {
		if (!isEqual(binodalAVaporArray[i].rhoB, binodalAVaporArray[i-1].rhoB, binodalBreakStep) &&
				!isEqual(binodalAVaporArray[i-1].rhoB, 0.0, zeroThreshold)) {
			binodalAVaporBreak = i;
			break;
		}
	}
	for (int i = 1; i < binodalALiquidTail; i++) {
		if (!isEqual(binodalALiquidArray[i].rhoB, binodalALiquidArray[i-1].rhoB, binodalBreakStep) &&
				!isEqual(binodalALiquidArray[i-1].rhoB, 0.0, zeroThreshold)) {
			binodalALiquidBreak = i;
			break;
		}
	}
	for (int i = 1; i < binodalARichTail; i++) {
		if (!isEqual(binodalARichArray[i].rhoA, binodalARichArray[i-1].rhoA, binodalBreakStep) &&
				!isEqual(binodalARichArray[i-1].rhoA, 0.0, zeroThreshold)) {
			binodalARichBreak = i;
			break;
		}
	}
	for (int i = 1; i < binodalBVaporTail; i++) {
		if (!isEqual(binodalBVaporArray[i].rhoA, binodalBVaporArray[i-1].rhoA, binodalBreakStep) &&
				!isEqual(binodalBVaporArray[i-1].rhoA, 0.0, zeroThreshold)) {
			binodalBVaporBreak = i;
			break;
		}
	}
	for (int i = 1; i < binodalBLiquidTail; i++) {
		if (!isEqual(binodalBLiquidArray[i].rhoA, binodalBLiquidArray[i-1].rhoA, binodalBreakStep) &&
				!isEqual(binodalBLiquidArray[i-1].rhoA, 0.0, zeroThreshold)) {
			binodalBVaporBreak = i;
			break;
		}
	}
	for (int i = 1; i < binodalBRichTail; i++) {
		if (!isEqual(binodalBRichArray[i].rhoB, binodalBRichArray[i-1].rhoB, binodalBreakStep) &&
				!isEqual(binodalBRichArray[i-1].rhoB, 0.0, zeroThreshold)) {
			binodalARichBreak = i;
			break;
		}
	}

	// Read out arrays to binodal file... done in 2 stages to account for line breaks that may have been detected above
	FILE *phaseDiagram_binodal = openFile("threePhaseDiagram-binodal", ".dat", "w");
	fprintf(phaseDiagram_binodal, "#! B-vapor binodal\n");
	for (int i = 0; i < binodalBVaporBreak; i++) fprintf(phaseDiagram_binodal, "%.15f %.15f\n", binodalBVaporArray[i].rhoA, binodalBVaporArray[i].rhoB);
	if (binodalBVaporBreak != binodalBVaporTail) {
		fprintf(phaseDiagram_binodal, "\n\n#! B-vapor binodal - post-break\n");
		for (int i = binodalBVaporBreak; i <= binodalBVaporTail; i++) {
			fprintf(phaseDiagram_binodal, "%.15f %.15f\n", binodalBVaporArray[i].rhoA, binodalBVaporArray[i].rhoB);
		}
	}
	fprintf(phaseDiagram_binodal, "\n\n#! A-vapor binodal\n");
	for (int i = 0; i < binodalAVaporBreak; i++) fprintf(phaseDiagram_binodal, "%.15f %.15f\n", binodalAVaporArray[i].rhoA, binodalAVaporArray[i].rhoB);
	if (binodalAVaporBreak != binodalAVaporTail) {
		fprintf(phaseDiagram_binodal, "\n\n#! A-vapor binodal - post-break\n");
		for (int i = binodalAVaporBreak; i <= binodalAVaporTail; i++) {
			fprintf(phaseDiagram_binodal, "%.15f %.15f\n", binodalAVaporArray[i].rhoA, binodalAVaporArray[i].rhoB);
		}
	}
	fprintf(phaseDiagram_binodal, "\n\n#! B-liquid binodal\n");
	for (int i = 0; i < binodalBLiquidBreak; i++) fprintf(phaseDiagram_binodal, "%.15f %.15f\n", binodalBLiquidArray[i].rhoA, binodalBLiquidArray[i].rhoB);
	if (binodalBLiquidBreak != binodalBLiquidTail) {
		fprintf(phaseDiagram_binodal, "\n\n#! B-liquid binodal - post-break\n");
		for (int i = binodalBLiquidBreak; i <= binodalBLiquidTail; i++) {
			fprintf(phaseDiagram_binodal, "%.15f %.15f\n", binodalBLiquidArray[i].rhoA, binodalBLiquidArray[i].rhoB);
		}
	}
	fprintf(phaseDiagram_binodal, "\n\n#! B-binary liquid binodal\n");
	for (int i = 0; i < binodalBRichBreak; i++) fprintf(phaseDiagram_binodal, "%.15f %.15f\n", binodalBRichArray[i].rhoA, binodalBRichArray[i].rhoB);
	if (binodalBRichBreak != binodalBRichTail) {
		fprintf(phaseDiagram_binodal, "\n\n#! B-binary liquid binodal - post-break\n");
		for (int i = binodalBRichBreak; i <= binodalBRichTail; i++) {
			fprintf(phaseDiagram_binodal, "%.15f %.15f\n", binodalBRichArray[i].rhoA, binodalBRichArray[i].rhoB);
		}
	}
	fprintf(phaseDiagram_binodal, "\n\n#! A-liquid binodal\n");
	for (int i = 0; i < binodalALiquidBreak; i++) fprintf(phaseDiagram_binodal, "%.15f %.15f\n", binodalALiquidArray[i].rhoA, binodalALiquidArray[i].rhoB);
	if (binodalALiquidBreak != binodalALiquidTail) {
		fprintf(phaseDiagram_binodal, "\n\n#! A-liquid binodal - post-break\n");
		for (int i = binodalALiquidBreak; i <= binodalALiquidTail; i++) {
			fprintf(phaseDiagram_binodal, "%.15f %.15f\n", binodalALiquidArray[i].rhoA, binodalALiquidArray[i].rhoB);
		}
	}
	fprintf(phaseDiagram_binodal, "\n\n#! A-binary liquid binodal\n");
	for (int i = 0; i < binodalARichBreak; i++) fprintf(phaseDiagram_binodal, "%.15f %.15f\n", binodalARichArray[i].rhoA, binodalARichArray[i].rhoB);
	if (binodalARichBreak != binodalARichTail) {
		fprintf(phaseDiagram_binodal, "\n\n#! A-binary liquid binodal - post-break\n");
		for (int i = binodalARichBreak; i <= binodalARichTail; i++) {
			fprintf(phaseDiagram_binodal, "%.15f %.15f\n", binodalARichArray[i].rhoA, binodalARichArray[i].rhoB);
		}
	}
	fclose(phaseDiagram_binodal);

	// Approximate the unconditionally unstable region by connecting end points of the binodal lines; A-rich point is always listed first in a pair
	FILE *unstable_region_coords = openFile("unconditionally-unstable-region", ".dat", "w");
	int minAIndex = 0, minBIndex = 0;
	double minRhoA = binodalARichArray[0].rhoA, minRhoB = binodalBRichArray[0].rhoB;
	for (int i = 0; i < binodalARichTail; i++) { // B-side binary liquid binodal "wraps" back on itself; need to find minA in A-rich and minB in B-rich to construct correct line here
		if (binodalARichArray[i].rhoA < minRhoA) {
			minAIndex = i;
			minRhoA = binodalARichArray[i].rhoA;
		}
		if (binodalBRichArray[i].rhoB < minRhoB) {
			minBIndex = i;
			minRhoB = binodalBRichArray[i].rhoB;
		}
	}
	fprintf(unstable_region_coords, "%f %f\n%f %f\n\n",
			binodalARichArray[minAIndex].rhoA, binodalARichArray[minAIndex].rhoB, binodalBRichArray[minBIndex].rhoA, binodalBRichArray[minBIndex].rhoB);
	fprintf(unstable_region_coords, "%f %f\n%f %f\n\n", binodalALiquidArray[binodalAVaporTail-1].rhoA, binodalALiquidArray[binodalAVaporTail-1].rhoB,
			binodalAVaporArray[binodalALiquidTail-1].rhoA, binodalAVaporArray[binodalALiquidTail-1].rhoB);
	fprintf(unstable_region_coords, "%f %f\n%f %f\n\n", binodalBVaporArray[binodalBVaporTail-1].rhoA, binodalBVaporArray[binodalBVaporTail-1].rhoB,
			binodalBLiquidArray[binodalBLiquidTail-1].rhoA, binodalBLiquidArray[binodalBLiquidTail-1].rhoB);
	fclose(unstable_region_coords);

	free(binodalAVaporArray);
	free(binodalALiquidArray);
	free(binodalARichArray);
	free(binodalBVaporArray);
	free(binodalBLiquidArray);
	free(binodalBRichArray);
	printf("Binodal lines complete\n");

	// With unconditionally unstable and metastable 3-phase regions defined, sift through minimization results to tag points inside those regions
	sortMetastableMinimizationResults();
} // end function generatePhaseDiagramBinodalLines()


/**
 * @brief Function _calculatePressure_ calculates the bulk pressure for a given density pair.
 *
 * Given an density pair, this function calculates the bulk pressure that density pair.
 * The function is written generically for a two-component mixture.  A single component chemical potential can be determined by setting one
 * of the density coordinates to zero (e.g. rhoA={non-zero}, rhoB=0.0 gives a single-component result).
 *
 * @param [in] rhoA The A-component coordinate of the mixture density.
 * @param [in] rhoB The B-component coordinate of the mixture density.
 *
 * @return The pressure value.
 */
double calculatePressure(double rhoA, double rhoB) {
	return (rhoA+rhoB)*theta*(1 + (bA*rhoA+bB*rhoB)/(1.-bA*rhoA-bB*rhoB)) - aA*rhoA*rhoA - 2*aAB*rhoA*rhoB - aB*rhoB*rhoB;
} // end function calculatePressure()


/**
 * @brief Function _calculateChemicalPotential_ calculates the bulk chemical potential for a given density pair.
 *
 * Given an density pair and a specified component, this function calculates the bulk chemical potential for that component at that density pair.
 * The function is written generically for a two-component mixture.  A single component chemical potential can be determined by setting the
 * cross-density coordinate to zero (e.g. *component="A", rhoA={non-zero}, rhoB=0.0 gives a single-component result).
 *
 * @note This function doesn't guard against calling with a component whose corresponding density is zero (e.g. *component="A", rhoA=0.0 is undefined).
 *
 * @param [in] *component A string literal equal to "A" or "B" to specify the desired chemical potential component. Values other than "A" or "B" are undefined.
 * @param [in] rhoA The A-component coordinate of the mixture density.
 * @param [in] rhoB The B-component coordinate of the mixture density.
 *
 * @return The chemical potential value.
 */
double calculateChemicalPotential(char *component, double rhoA, double rhoB) {
	double a = 0.1, b = 1./3.;
	double rho1 = 0.0, rho2 = 0.0;
	double excludedVolume = 0.0;
	double zeroThreshold = 1e-14;

	// Set parameters according to the desired component
	if (strcmp(component, "A")) {
		rho1 = rhoA;
		rho2 = rhoB;
		a = aA;
		b = bA;
	}
	else if (strcmp(component, "B")) {
		rho1 = rhoB;
		rho2 = rhoA;
		a = aB;
		b = bB;
	}

	if (isEqual(rho2, 0.0, zeroThreshold)) { // single component case
		excludedVolume = 1.0 - b*rho1;
		return theta*log(rho1/excludedVolume) + theta*b*rho1/excludedVolume + theta - 2.0*a*rho1;
	}
	else { // two component case
		excludedVolume = 1.0 - bA*rhoA - bB*rhoB;
		return theta*log(rho1/excludedVolume) + theta*b*(rho1+rho2)/excludedVolume + theta - 2.0*a*rho1 - 2.0*aAB*rho2;
	}
} // end function calculateChemicalPotential()


/**
 * @brief Function _generateMinimizedPressureAndMuDeviations_ calculates deviations among the theoretical pressure and chemical potential values for
 * each phase of a free energy minimization to confirm correct equilibrium behavior.
 *
 * This function is used to calculate the pressures and chemical potentials of the densities found in free energy minimization.
 * It calculates absolute value differences between the 2 or 3 phases and logs them for plotting.
 * A difference of zero means the 2 phases have the same pressure and/or chemical potential in equilibrium (way it should be).
 * Only calculates pressures and chemical potentials in the bulk phases.
 *
 * @note The negative of the maximum pressure/chemical potential deviation of the 2-phase data separates the 2- from 3-phase portions of each data set.
 *
 * @todo Could functionalize the logic common to the 2-phase and 3-phase data handling loops.
 */
void generateMinimizedPressureAndMuDeviations() {
	int readEOF = 0;
	double particlesATotal = 0.0, particlesBTotal = 0.0;
	double rhoA1 = 0.0, rhoB1 = 0.0, rhoA2 = 0.0, rhoB2 = 0.0, rhoA3 = 0.0, rhoB3 = 0.0;
	double pressure1 = 0.0, pressure2 = 0.0, pressure3 = 0.0;
	double muA1 = 0.0, muB1 = 0.0, muA2 = 0.0, muB2 = 0.0, muA3 = 0.0, muB3 = 0.0;
	double pressureDeviation = 0.0, muADeviation = 0.0, muBDeviation = 0.0;
	double maxPressureDeviation = 0.0, maxMuADeviation = 0.0, maxMuBDeviation = 0.0;
	FILE *phaseDiagram_densities_twoPhases = openFile("threePhaseDiagram-densities-twoPhases", ".dat", "r");
	FILE *phaseDiagram_densities_threePhases = openFile("threePhaseDiagram-densities-threePhases", ".dat", "r");
	FILE *pressure_deviations_theoretical = openFile("pressure-deviations", ".dat", "w");
	FILE *mu_deviations_theoretical = openFile("chemical-potential-deviations", ".dat", "w");

	printf("processing 2-phase densities...\n");
	readEOF = fscanf(phaseDiagram_densities_twoPhases, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
			&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
	while (readEOF != EOF) {

		//
		// Calculate bulk chemical potentials for each phase
		//
		if (particlesATotal == 0) {
			muA1 = 0;
			muA2 = 0;
			muA3 = 0;
			muB1 = calculateChemicalPotential("B", 0.0, rhoB1);
			muB2 = calculateChemicalPotential("B", 0.0, rhoB2);
			muB3 = calculateChemicalPotential("B", 0.0, rhoB3);
		}
		else if (particlesBTotal == 0) {
			muA1 = calculateChemicalPotential("A", rhoA1, 0.0);
			muA2 = calculateChemicalPotential("A", rhoA2, 0.0);
			muA3 = calculateChemicalPotential("A", rhoA3, 0.0);
			muB1 = 0;
			muB2 = 0;
			muB3 = 0;
		}
		else {
			muA1 = calculateChemicalPotential("A", rhoA1, rhoB1);
			muB1 = calculateChemicalPotential("B", rhoA1, rhoB1);
			muA2 = calculateChemicalPotential("A", rhoA2, rhoB2);
			muB2 = calculateChemicalPotential("B", rhoA2, rhoB2);
			muA3 = calculateChemicalPotential("A", rhoA3, rhoB3);
			muB3 = calculateChemicalPotential("B", rhoA3, rhoB3);
		}

		// Determine chemical potential deviations among phases
		muADeviation = fabs(muA1-muA2);
		if (fabs(muA1-muA3) > muADeviation) muADeviation = fabs(muA1-muA3);
		if (fabs(muA2-muA3) > muADeviation) muADeviation = fabs(muA2-muA3);
		if (muADeviation > maxMuADeviation) maxMuADeviation = muADeviation;

		muBDeviation = fabs(muB1-muB2);
		if (fabs(muB1-muB3) > muBDeviation) muBDeviation = fabs(muB1-muB3);
		if (fabs(muB2-muB3) > muBDeviation) muBDeviation = fabs(muB2-muB3);
		if (muBDeviation > maxMuBDeviation) maxMuBDeviation = muBDeviation;

		//
		// Determine bulk pressures
		//
		pressure1 = calculatePressure(rhoA1, rhoB1);
		pressure2 = calculatePressure(rhoA2, rhoB2);
		pressure3 = calculatePressure(rhoA3, rhoB3);

		// Determine pressure deviations among phases
		pressureDeviation = fabs(pressure1-pressure2);
		if (fabs(pressure1-pressure3) > pressureDeviation) pressureDeviation = fabs(pressure1-pressure3);
		if (fabs(pressure2-pressure3) > pressureDeviation) pressureDeviation = fabs(pressure2-pressure3);
		if (pressureDeviation > maxPressureDeviation) maxPressureDeviation = pressureDeviation;

		fprintf(pressure_deviations_theoretical,"%.15f %.15f %.15f %.15f %.15f %.15f\n",
				particlesATotal, particlesBTotal, pressureDeviation, pressure1, pressure2, pressure3);
		fprintf(mu_deviations_theoretical,"%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n",
				particlesATotal, particlesBTotal, muADeviation, muBDeviation, muA1, muB1, muA2, muB2, muA2, muB3);
		readEOF = fscanf(phaseDiagram_densities_twoPhases, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
				&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
	}

	// Mark the data sets with negative of the max deviations to divide the 2-phase from the 3-phase densities
	fprintf(pressure_deviations_theoretical,"%f %f %.15f\n", 0.0, 0.0, -maxPressureDeviation);
	fprintf(mu_deviations_theoretical,"%f %f %.15f %.15f\n", 0.0, 0.0, -maxMuADeviation, -maxMuBDeviation);

	printf("processing 3-phase densities...\n");
	readEOF = fscanf(phaseDiagram_densities_threePhases, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
			&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
	while (readEOF != EOF) {

		//
		// Calculate bulk chemical potentials for each phase
		//
		if (particlesATotal == 0) {
			muA1 = 0;
			muA2 = 0;
			muA3 = 0;
			muB1 = calculateChemicalPotential("B", 0.0, rhoB1);
			muB2 = calculateChemicalPotential("B", 0.0, rhoB2);
			muB3 = calculateChemicalPotential("B", 0.0, rhoB3);
		}
		else if (particlesBTotal == 0) {
			muA1 = calculateChemicalPotential("A", rhoA1, 0.0);
			muA2 = calculateChemicalPotential("A", rhoA2, 0.0);
			muA3 = calculateChemicalPotential("A", rhoA3, 0.0);
			muB1 = 0;
			muB2 = 0;
			muB3 = 0;
		}
		else {
			muA1 = calculateChemicalPotential("A", rhoA1, rhoB1);
			muB1 = calculateChemicalPotential("B", rhoA1, rhoB1);
			muA2 = calculateChemicalPotential("A", rhoA2, rhoB2);
			muB2 = calculateChemicalPotential("B", rhoA2, rhoB2);
			muA3 = calculateChemicalPotential("A", rhoA3, rhoB3);
			muB3 = calculateChemicalPotential("B", rhoA3, rhoB3);
		}

		// Determine chemical potential deviations among phases
		muADeviation = fabs(muA1-muA2);
		if (fabs(muA1-muA3) > muADeviation) muADeviation = fabs(muA1-muA3);
		if (fabs(muA2-muA3) > muADeviation) muADeviation = fabs(muA2-muA3);

		muBDeviation = fabs(muB1-muB2);
		if (fabs(muB1-muB3) > muBDeviation) muBDeviation = fabs(muB1-muB3);
		if (fabs(muB2-muB3) > muBDeviation) muBDeviation = fabs(muB2-muB3);

		//
		// Determine bulk pressures
		//
		pressure1 = calculatePressure(rhoA1, rhoB1);
		pressure2 = calculatePressure(rhoA2, rhoB2);
		pressure3 = calculatePressure(rhoA3, rhoB3);

		// Determine pressure deviations among phases
		pressureDeviation = fabs(pressure1-pressure2);
		if (fabs(pressure1-pressure3) > pressureDeviation) pressureDeviation = fabs(pressure1-pressure3);
		if (fabs(pressure2-pressure3) > pressureDeviation) pressureDeviation = fabs(pressure2-pressure3);

		fprintf(pressure_deviations_theoretical,"%.15f %.15f %.15f %.15f %.15f %.15f\n",
				particlesATotal, particlesBTotal, pressureDeviation, pressure1, pressure2, pressure3);
		fprintf(mu_deviations_theoretical,"%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n",
				particlesATotal, particlesBTotal, muADeviation, muBDeviation, muA1, muB1, muA2, muB2, muA2, muB3);
		readEOF = fscanf(phaseDiagram_densities_threePhases, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
				&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
	}

	fclose(phaseDiagram_densities_twoPhases);
	fclose(phaseDiagram_densities_threePhases);
	fclose(pressure_deviations_theoretical);
	fclose(mu_deviations_theoretical);
	printf("Pressure and chemical potential deviations complete\n\n");
} // end function generateMinimizedPressureAndMuDeviations()


/**
 * @brief Function _determinePhaseSeparation_ takes a single test point and initializes the free energy minimization for that point to determine its phase separated state.
 *
 * This function orchestrates the free energy minimization process for a single (A,B) test point.  It uses the results of a stability analysis at that point to initialize the
 * particle count and volumina for each of 3 possible phases and calls the minimizing function.  If the minimizing function returns a state where one phase is still in an unstable
 * region, this function will "reminimize" by reallocating particles and volumina such that the stable phases are packed into one and the unstable "stuck" phase is split into 2 along
 * its unstable eigenvector.  If the reminimization fails 100 times, this function will terminiate the minization and not log any information for the unsuccessful point.
 *
 * This function is run for two uses, and each use will write the results to a different location.  When minimizing all points of a phase diagram automatically, this function
 * writes to a log file.  When using a utility to minimize just a single test point (presumably for manual troubleshooting purposes), this function writes to stdout.
 *
 * @param [in] particlesA The A-component coordinate of the test point.
 * @param [in] particlesB The B-component coordinate of the test point.
 * @param [in/out] *densities_tmp A pointer to the location to write the results of the phase separation (stdout or "densities_tmp" log file).
 */
void determinePhaseSeparation(double particlesA, double particlesB, FILE *densities_tmp) {
	int phaseChangeOccurred = 0;
	int twoOrThreePhases = 0;
	int minimizationTake2 = 0, minimizationIteration = 0, maxMinimizationAttempts = 100;
	double eigenvectorA = 0.0, eigenvectorB = 0.0;
	double freeEnergy = 0.0, freeEnergy1 = 0.0, freeEnergy2 = 0.0, freeEnergy3 = 0.0;
	double particlesA1 = 0.0, particlesB1 = 0.0;
	double particlesA2 = 0.0, particlesB2 = 0.0;
	double particlesA3 = 0.0, particlesB3 = 0.0;
	double volume1 = 0.0, volume2 = 0.0, volume3 = 0.0, volumeTotal = 1.0;
	double rhoA1 = 0, rhoA2 = 0, rhoA3 = 0;
	double rhoB1 = 0, rhoB2 = 0, rhoB3 = 0;
	double tmpParticlesA = 0.0, tmpParticlesB = 0.0, tmpVolume = 0.0;
	double phaseDetectionThreshold = 5e-3; // also tested phaseDetection 1e-3, 1e-4, 0.01, 1e-9, 1e-5
	double trialStepSize = minimizationStepFactor * minimizationParticlesStepSize;

	volume1 = volumeTotal / 3.0;
	volume2 = volumeTotal / 3.0;

	// Only try to minimize if the stability analysis indicates it should minimize
	if (calculateMinimizationPath(particlesA, particlesB, &eigenvectorA, &eigenvectorB)) {
		particlesA1 = (particlesA + trialStepSize*eigenvectorA) * volume1;
		particlesB1 = (particlesB - trialStepSize*eigenvectorB) * volume1;
		particlesA2 = (particlesA - trialStepSize*eigenvectorA) * volume2;
		particlesB2 = (particlesB + trialStepSize*eigenvectorB) * volume2;

		phaseChangeOccurred = minimizeFreeEnergyTwoComponentsThreePhases(theta, particlesA, particlesB, volumeTotal,
				&particlesA1, &particlesB1, &volume1, &particlesA2, &particlesB2, &volume2, &freeEnergy, &freeEnergy1, &freeEnergy2, &freeEnergy3);
		minimizationIteration = 1;
		while (phaseChangeOccurred && minimizationIteration < maxMinimizationAttempts) {

			// Determine if the 1st minimization resulted in 2 or 3 phases, and write all to a temp file
			if (minimizationTake2 == 0) { // only write here on the 1st minimization attempt... subsequent attempts are handled below
				twoOrThreePhases = detectTwoOrThreePhaseRegion(phaseDetectionThreshold, densities_tmp, particlesA, particlesB, volumeTotal,
						particlesA1, particlesB1, volume1, particlesA2, particlesB2, volume2, freeEnergy, freeEnergy1, freeEnergy2, freeEnergy3);
			}

			// Repack particles/volumina and reminimize if the minimization went to an unstable section of the phase diagram
			// Otherwise set control flags to just move on to the next test point when a good 2- or 3-phase minimization has occurred
			if (twoOrThreePhases == -1) { // phase 1 still has a - eigenvalue
				printf("\nPhase 1 -eigenvalue: 1=(%f,%f) 2=(%f,%f) 3=(%f,%f)\n",
						particlesA1/volume1, particlesB1/volume1, particlesA2/volume2, particlesB2/volume2,
						(particlesA-particlesA1-particlesA2)/(volumeTotal-volume1-volume2),
						(particlesB-particlesB1-particlesB2)/(volumeTotal-volume1-volume2));
				minimizationTake2 = 1;
				tmpParticlesA = particlesA1;
				tmpParticlesB = particlesB1;
				tmpVolume = volume1;

				if (minimizationIteration < maxMinimizationAttempts-1) { // use phase 2 as the free phase; split phase 1 volume evenly between the 2
					volume2 = 0.5 * volume1;
					volume1 *= 0.5;
				}
				else { // these should be "roughly" equal; pack all into the smaller of the remaining volumes...
					if (volume2 < volume3) {
						particlesA2 += particlesA1;
						particlesB2 += particlesB1;
						volume2 += volume1;
					}
					else {
						particlesA3 += particlesA1;
						particlesB3 += particlesB1;
						volume3 += volume1;
					}
					particlesA1 = 0.0; // ... and zero out phase 1 to force a 2-phase result and terminate the loop
					particlesB1 = 0.0;
					volume1 = 0.0;
				}
			}
			else if (twoOrThreePhases == -2) { // phase 2 still has a - eigenvalue
				printf("\nPhase 2 -eigenvalue: 1=(%f,%f) 2=(%f,%f) 3=(%f,%f)\n",
						particlesA1/volume1, particlesB1/volume1, particlesA2/volume2, particlesB2/volume2,
						(particlesA-particlesA1-particlesA2)/(volumeTotal-volume1-volume2),
						(particlesB-particlesB1-particlesB2)/(volumeTotal-volume1-volume2));
				minimizationTake2 = 1;
				tmpParticlesA = particlesA2;
				tmpParticlesB = particlesB2;
				tmpVolume = volume2;

				if (minimizationIteration < maxMinimizationAttempts-1) { // use phase 1 as the free phase; split phase 2 volume evenly between the 2
					volume1 = 0.5 * volume2;
					volume2 *= 0.5;
				}
				else { // these should be "roughly" equal; pack all into the smaller of the remaining volumes...
					if (volume1 < volume3) {
						particlesA1 += particlesA2;
						particlesB1 += particlesB2;
						volume1 += volume2;
					}
					else {
						particlesA3 += particlesA2;
						particlesB3 += particlesB2;
						volume3 += volume2;
					}
					particlesA2 = 0.0; // ... and zero out phase 2 to force a 2-phase result and terminate the loop
					particlesB2 = 0.0;
					volume2 = 0.0;
				}
			}
			else if (twoOrThreePhases == -3) { // phase 3 still has a - eigenvalue
				printf("\nPhase 3 -eigenvalue: 1=(%f,%f) 2=(%f,%f) 3=(%f,%f)\n",
						particlesA1/volume1, particlesB1/volume1, particlesA2/volume2, particlesB2/volume2,
						(particlesA-particlesA1-particlesA2)/(volumeTotal-volume1-volume2),
						(particlesB-particlesB1-particlesB2)/(volumeTotal-volume1-volume2));
				minimizationTake2 = 1;
				tmpParticlesA = particlesA - particlesA1 - particlesA2;
				tmpParticlesB = particlesB - particlesB1 - particlesB2;
				tmpVolume = volumeTotal - volume1 - volume2;

				if (minimizationIteration < maxMinimizationAttempts-1) { // divide the phase 3 volume between phases 1 & 2
					volume1 = 0.5 * tmpVolume;
					volume2 = volume1;
				}
				else { // these should be "roughly" equal; pack all into the smaller of the remaining volumes...
					if (volume1 < volume2) {
						particlesA1 += particlesA3;
						particlesB1 += particlesB3;
						volume1 += volume3;
					}
					else {
						particlesA2 += particlesA3;
						particlesB2 += particlesB3;
						volume2 += volume3;
					}
					particlesA3 = 0.0; // ... and zero out phase 3 to force a 2-phase result and terminate the loop
					particlesB3 = 0.0;
					volume3 = 0.0;
				}
			}
			else if (twoOrThreePhases == 2) {
				minimizationTake2 = 0;
				break; // good minimization; results written on previous loop iteration, so break the reminimization while loop
			}
			else if (twoOrThreePhases == 3) {
				minimizationTake2 = 0;
				threePhaseRegionExists = 1;
				rhoA1 = particlesA1 / volume1;
				rhoB1 = particlesB1 / volume1;
				rhoA2 = particlesA2 / volume2;
				rhoB2 = particlesB2 / volume2;
				rhoA3 = (particlesA-particlesA1-particlesA2) / (volumeTotal-volume1-volume2);
				rhoB3 = (particlesB-particlesB1-particlesB2) / (volumeTotal-volume1-volume2);
				checkVaporPhase(&rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3, particlesA, particlesB);
				if (defineThreePhaseRegion(rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3)) {
					rhoAThreePhaseRegion = particlesA;
					rhoBThreePhaseRegion = particlesB;
				}
				break; // good minimization; results written on previous loop iteration, so break the reminimization while loop
			}

			if (minimizationTake2) {
				minimizationIteration++;
				printf("Take %i... AT=%f BT=%f: phase A=%f B=%f\n",
						minimizationIteration, particlesA, particlesB, tmpParticlesA/tmpVolume, tmpParticlesB/tmpVolume);

				calculateMinimizationPath(tmpParticlesA/tmpVolume, tmpParticlesB/tmpVolume, &eigenvectorA, &eigenvectorB);

				// The repacked volume has already been split between phases 1/2... just split the particles now and point in the right direction
				particlesA1 = (0.5*tmpParticlesA + trialStepSize*eigenvectorA);
				particlesB1 = (0.5*tmpParticlesB - trialStepSize*eigenvectorB);
				particlesA2 = (0.5*tmpParticlesA - trialStepSize*eigenvectorA);
				particlesB2 = (0.5*tmpParticlesB + trialStepSize*eigenvectorB);
				printf("... repacking: A1=%f B1=%f V1=%f A2=%f B2=%f V2=%f A3=%f B3=%f V3=%f\n",
						particlesA1, particlesB1, volume1, particlesA2, particlesB2, volume2,
						(particlesA-particlesA1-particlesA2), (particlesB-particlesB1-particlesB2), (volumeTotal-volume1-volume2));
				printf("... new densities: 1=(%f,%f) 2=(%f,%f) 3=(%f,%f)\n",
						particlesA1/volume1, particlesB1/volume1, particlesA2/volume2, particlesB2/volume2,
						(particlesA-particlesA1-particlesA2)/(volumeTotal-volume1-volume2),
						(particlesB-particlesB1-particlesB2)/(volumeTotal-volume1-volume2));

				// On last attempt only, to preserve the zeroing out of a phase above in case it's still a - eigenvalue, skip the reminimization...
				if (minimizationIteration < maxMinimizationAttempts) {
					phaseChangeOccurred = minimizeFreeEnergyTwoComponentsThreePhases(theta, particlesA, particlesB, volumeTotal,
							&particlesA1, &particlesB1, &volume1, &particlesA2, &particlesB2, &volume2,
							&freeEnergy, &freeEnergy1, &freeEnergy2, &freeEnergy3);
				}
				else phaseChangeOccurred = 0; // ... and force the printout below to say "no effect"

				// Determine if the 2nd minimization resulted in 2 or 3 phases, and write all to a temp file
				// When zeroing out a phase on last minimization attempt, zeroes will force writing a 2-phase result... may not be an optimal result
				twoOrThreePhases = detectTwoOrThreePhaseRegion(phaseDetectionThreshold, densities_tmp, particlesA, particlesB, volumeTotal,
						particlesA1, particlesB1, volume1, particlesA2, particlesB2, volume2, freeEnergy, freeEnergy1, freeEnergy2, freeEnergy3);

				if (phaseChangeOccurred) {
					if (twoOrThreePhases == 3) printf("Qa'plah!... 1=(%f,%f) 2=(%f,%f) 3=(%f,%f)\n\n",
							particlesA1/volume1, particlesB1/volume1, particlesA2/volume2, particlesB2/volume2,
							(particlesA-particlesA1-particlesA2)/(volumeTotal-volume1-volume2),
							(particlesB-particlesB1-particlesB2)/(volumeTotal-volume1-volume2));
					else if (twoOrThreePhases == 2) printf("\"metastable\" again... 1=(%f,%f) 2=(%f,%f) 3=(%f,%f)\n\n",
							particlesA1/volume1, particlesB1/volume1, particlesA2/volume2, particlesB2/volume2,
							(particlesA-particlesA1-particlesA2)/(volumeTotal-volume1-volume2),
							(particlesB-particlesB1-particlesB2)/(volumeTotal-volume1-volume2));
					else printf("hmm, another - eigenvalue... %i: 1=(%f,%f) 2=(%f,%f) 3=(%f,%f)\n",
							twoOrThreePhases, particlesA1/volume1, particlesB1/volume1, particlesA2/volume2, particlesB2/volume2,
							(particlesA-particlesA1-particlesA2)/(volumeTotal-volume1-volume2),
							(particlesB-particlesB1-particlesB2)/(volumeTotal-volume1-volume2));
				}
				else printf("... reminimizing had no effect: 1=(%f,%f) 2=(%f,%f) 3=(%f,%f)\n\n",
						particlesA1/volume1, particlesB1/volume1, particlesA2/volume2, particlesB2/volume2,
						(particlesA-particlesA1-particlesA2)/(volumeTotal-volume1-volume2),
						(particlesB-particlesB1-particlesB2)/(volumeTotal-volume1-volume2));
			} // end if (minimizationTake2)
		} // end while (phaseChangeOccurred)
	} // end if (calculateMinimizationPath)
	else if (densities_tmp == stdout) printf("(%f,%f) is not expected to phase separate.\n\n", particlesA, particlesB);
} // end function determinePhaseSeparation()


/**
 * @brief Function _calculatePhaseDiagramRhoAVsRhoBThreePhasesTheoretical_ calculates a theoretical binary phase diagram that allows for 3 phases to co-exist.
 *
 * This function generates a theoretical phase diagram by looping over all possible (A,B) density pairs and performing a minimization of the
 * Helmholtz free energy for a binary van der Waals fluid mixture.  The limits of the density pairs to test are defined by the global parameter _b_,
 * which is proportional to the diameter of a hard-core particle (i.e. related to excluded volume).  Using this as the limit avoids a singularity in
 * the free energy calculation.
 *
 * @note The global parameter _b_ is defaults to equal to 1/3.
 * @note The local parameter _phaseDetectionThreshold_ controls the sensitivity that determines if a test point is classified as a 2-phase or 3-phase
 * minimization (tighter thresholds will classify more points as 3-phase).
 */
void calculatePhaseDiagramRhoAVsRhoBThreePhasesTheoretical(){
	int twoOrThreePhases = 0;
	int readEOF = 0;
	double freeEnergy = 0.0, freeEnergy1 = 0.0, freeEnergy2 = 0.0, freeEnergy3 = 0.0;
	double particlesATotal = 0.0, particlesBTotal = 0.0;
	double particlesA1 = 0.0, particlesB1 = 0.0;
	double particlesA2 = 0.0, particlesB2 = 0.0;
	double particlesA3 = 0.0, particlesB3 = 0.0;
	double volume1 = 0.0, volume2 = 0.0, volume3 = 0.0, volumeTotal = 1.0;
	double rhoA1 = 0, rhoA2 = 0, rhoA3 = 0;
	double rhoB1 = 0, rhoB2 = 0, rhoB3 = 0;
	double particleSizeThreshold = 1e-6;
	double currentColumn = -minimizationParticlesStepSize;
	double rhoAThreePhaseRegion = 0.0, rhoBThreePhaseRegion = 0.0; // track the A,B point that defines the 3-phase region
	FILE *phaseDiagram_densities_twoPhases = openFile("threePhaseDiagram-densities-twoPhases", ".dat", "w");
	FILE *phaseDiagram_densities_threePhases = openFile("threePhaseDiagram-densities-threePhases", ".dat", "w");
	FILE *phaseDiagram_tieLines = openFile("threePhaseDiagram-tie-lines-theoretical", ".dat", "w");
	FILE *densities_tmp = openFile("densities-tmp", ".dat", "w+");

	time_t clockRunTime;
	clock_t processorRunTime;
	clockRunTime = time(NULL);
	processorRunTime = clock();
	printf("\n-----------------------------------------------------------------------------------------\n"
			"Minimizing free energy... 2 components, 3 phases allowed...\n");
	printf("tcA=%f\ttcB=%f\tncB=%f\tVDW interaction=%f\n\n", tcA, tcB, ncB, vdwInteractionFactor);

	// Write line coordinates of points connecting component A and component B excluded volume b parameters
	FILE *singularity_coords = openFile("singularity-line-coords", ".dat", "w");
	fprintf(singularity_coords, "%.15f %.15f\n%.15f %.15f",	0.0, 1/bB, 1/bA, 0.0);
	fclose(singularity_coords);

	for (particlesATotal = 0; particlesATotal < (1/bA-particleSizeThreshold); particlesATotal += minimizationParticlesStepSize) { // columns
		for (particlesBTotal = 0; (particlesBTotal+particlesATotal) < (1/bB-particleSizeThreshold); particlesBTotal += minimizationParticlesStepSize) { // rows
			if (!isEqual(currentColumn, particlesATotal, particleSizeThreshold)) {
				currentColumn = particlesATotal;
				printf("Scanning column nA0=%f\n", particlesATotal);
			}
			determinePhaseSeparation(particlesATotal, particlesBTotal, densities_tmp);
		} // end for loop (B particles)
	} // end for loop (A particles)

	//
	// Write the phase-separated points to log files specific to the 2- or 3-phase regions of the phase diagram
	// ... could functionalize this easy, but I'm too lazy right now to re-define all of the local variables in a new function that are already here  :-)
	//
	fseek(densities_tmp, 0, SEEK_SET); // reset the file pointer to the start
	readEOF = fscanf(densities_tmp, "%1i %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf "
			"%9lf %9lf %9lf %9lf",
			&twoOrThreePhases,
			&particlesATotal, &particlesBTotal, &volumeTotal,
			&particlesA1, &particlesB1, &volume1, &rhoA1, &rhoB1,
			&particlesA2, &particlesB2, &volume2, &rhoA2, &rhoB2,
			&particlesA3, &particlesB3, &volume3, &rhoA3, &rhoB3,
			&freeEnergy, &freeEnergy1, &freeEnergy2, &freeEnergy3);
	while (readEOF != EOF) {
		checkVaporPhase(&rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3, particlesATotal, particlesBTotal);

		if (twoOrThreePhases == 2) {
			fprintf(phaseDiagram_densities_twoPhases, "%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n",
					particlesATotal, particlesBTotal, rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3, volume1, volume2, volume3);
		}
		else if (twoOrThreePhases == 3) {
			if (rhoA3 == 0 && rhoB3 == 0) {
				printf("(%.15f,%.15f) identified 3-phase; appears to be meta-stable with 2 phases\n", particlesATotal, particlesBTotal);
				printf("rhoA1 = %.18f\trhoB1 = %.18f\n", rhoA1, rhoB1);
				printf("rhoA2 = %.18f\trhoB2 = %.18f\n", rhoA2, rhoB2);
				printf("rhoA3 = %.18f\trhoB3 = %.18f\n", rhoA3, rhoB3);
				printf("... assigning 3 = 1\n\n");
				rhoA3 = rhoA1;
				rhoB3 = rhoB1;
			}
			fprintf(phaseDiagram_densities_threePhases, "%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n",
					particlesATotal, particlesBTotal, rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3, volume1, volume2, volume3);
		}
		else printf("Unphysical phase detected... %i %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n",
				twoOrThreePhases, particlesATotal, particlesBTotal, rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);

		// Automatically generate tie lines and write
		if (tieLineNeeded(particlesATotal, particlesBTotal)) {
			fprintf(phaseDiagram_tieLines, "%.15f %.15f\n%.15f %.15f\n%.15f %.15f\n\n", rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);
		}

		#ifdef DEBUG_TIELINES_ON
		FILE *single_tieLine_points;
		FILE *single_tieLine;
		char single_tieLine_points_name[255];
		sprintf(single_tieLine_points_name, "/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/tie-lines/tieLines-points_%f_%f.dat",
				particlesATotal, particlesBTotal);
		single_tieLine_points = fopen(single_tieLine_points_name, "w");
		fprintf(single_tieLine_points, "%e %e %e %e %e %e %e %e\n", particlesATotal, particlesBTotal, theoreticalDensityA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);
		fclose(single_tieLine_points);

		char single_tieLine_name[255];
		sprintf(single_tieLine_name, "/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/tie-lines/tieLines_%f_%f.dat",
				particlesATotal, particlesBTotal);
		single_tieLine = fopen(single_tieLine_name, "w");
		fprintf(single_tieLine, "%e %e\n%e %e\n%e %e\n", theoreticalDensityA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);
		fclose(single_tieLine);
		#endif

		readEOF = fscanf(densities_tmp, "%1i %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf "
				"%9lf %9lf %9lf %9lf",
				&twoOrThreePhases,
				&particlesATotal, &particlesBTotal, &volumeTotal,
				&particlesA1, &particlesB1, &volume1, &rhoA1, &rhoB1,
				&particlesA2, &particlesB2, &volume2, &rhoA2, &rhoB2,
				&particlesA3, &particlesB3, &volume3, &rhoA3, &rhoB3,
				&freeEnergy, &freeEnergy1, &freeEnergy2, &freeEnergy3);
	} // end while

	fclose(phaseDiagram_densities_twoPhases);
	fclose(phaseDiagram_densities_threePhases);
	fclose(phaseDiagram_tieLines);
	fclose(densities_tmp);

	writeThreePhaseRegionCoordinates(); // write 3-phase coordinates to a file so future program runs that do not minimize free energy don't need to re-define them
	printf("\n\n3-phase region defined by (%f,%f)\n", rhoAThreePhaseRegion, rhoBThreePhaseRegion);
	printf("%f %f %f %f %f %f\n", threePhaseRegion[0], threePhaseRegion[1], threePhaseRegion[2], threePhaseRegion[3], threePhaseRegion[4], threePhaseRegion[5]);
	printMinimizationResults(rhoAThreePhaseRegion, rhoBThreePhaseRegion);

	generateMinimizedPressureAndMuDeviations();
	generatePhaseDiagramBinodalLines();
	gnuplotTwoComponentThreePhase();
	logVDWParameters();

	printRunTime( (clock()-processorRunTime), difftime(time(NULL),clockRunTime) );
	printf("\nDone! Free energy minimization complete.\n"
			"--------------------------------------------------------------------------------------------------------------\n\n");
} // end function calculatePhaseDiagramRhoAVsRhoBThreePhasesTheoretical()


/**
 * @brief Function _minimizeSinglePoint_ is a wrapper function used to minimize a single user-specified test point.  The A- and B-coordinates of the point to minimize are manually
 * specified in the GUI.
 */
void minimizeSinglePoint() {
	determinePhaseSeparation(singlePointRhoA, singlePointRhoB, stdout);
} // end function minimizeSinglePoint()


/**
 * @brief Function _exploreVDWShieldRegion_ controls the creation of several phase diagrams around and within the so-called "shield region" of the binary VDW global phase diagram.
 */
void exploreVDWShieldRegion() {
//	// When bigLambda = {0.4817, 0.4666, 0.4515, 0.4364, 0.3921, 0.3478, 0.3035};
//	double nuArray[7] = {0.5183, 0.5334, 0.5485, 0.5636, 0.6079, 0.6522, 0.6965};
//
	ncB = 1.0;		// keep bA/bB symmetric; zeta = 0.0
//	tcA = 0.475;	// tc ratio = 1.0 when zeta = 0.0
//	tcB = 0.475;

	tcA = 0.4;		// tc ratio = 1.0 when zeta = 0.0
	tcB = 0.4;
	vdwInteractionFactor = 0.5636; // right in the middle
	for (int i = 0; i <= 50; i++) {
		//		vdwInteractionFactor = nuArray[i];

		setInitializeRandom(); // quick random init to make sure the log files have the correct variables in the filenames
		calculateKandSParameters();
		calculatePhaseDiagramRhoAVsRhoBThreePhasesTheoretical();

		tcA += 0.001;
		tcB += 0.001;
	}
} // end function exploreVDWShieldRegion()


/**
 * @brief Function _determineLBMetastableResults_ evaluates lattice Boltzmann simulation results for a single test point to determine if the test point
 * exhibits 2-phase or 3-phase behavior.
 *
 * This function calculates root sum square (RSS) error statistics for each of the three possible phases emerging from a LB simulation.  The RSS values
 * are compared to thresholds derived during testing to classify the results as exhibiting three types of 2-phase behavior (unstable, metastable,
 * binary liquid) or 3-phase behavior.
 *
 * @note This function assumes the global _threePhaseRegion[]_ coordinates are written in A-rich (indices 0,1), B-rich (2,3), vapor order (4,5).
 *
 * @return Integer value indicating if the simulation exhibits 2-phase liquid-vapor (0 or 1), 2-phase binary liquid (2), or 3-phase behavior (-1).
 *
 * @todo Metastable points with small kappa values on the very asymmetric diagrams may not fully form new phases... results get "stuck" inside the spinodal.
 */
int determineLBMetastableResults(double rhoAMaxA, double rhoBMaxA, double rhoAMaxB, double rhoBMaxB, double rhoAMinN, double rhoBMinN) {
	double phase1Separation = sqrt(fabs(pow(threePhaseRegion[0]-rhoAMaxA,2) + pow(threePhaseRegion[1]-rhoBMaxA,2)));
	double phase2Separation = sqrt(fabs(pow(threePhaseRegion[2]-rhoAMaxB,2) + pow(threePhaseRegion[3]-rhoBMaxB,2)));
	double phase3Separation = sqrt(fabs(pow(threePhaseRegion[4]-rhoAMinN,2) + pow(threePhaseRegion[5]-rhoBMinN,2)));
	printf("p1sep=%f p2sep=%f p3sep=%f\n", phase1Separation, phase2Separation, phase3Separation);

	// 2 of 3 phases at least must be below a separation threshold from the 3-phase theoretical values for 3-phase behavior to be possible
	// Even if all 3 are below the threshold, could still be a binary liquid case
	if ((phase1Separation < metastableThreshold || phase2Separation < metastableThreshold) && phase3Separation < metastableThreshold) {
		if (isEqual(phase1Separation, phase2Separation, 1e-3)) return 2; // binary liquid
		else if (phase1Separation > metastableThresholdFactor*metastableThreshold ||
				phase2Separation > metastableThresholdFactor*metastableThreshold) return 1; // 2-phase
		else return -1; // 3-phase
	}
	else return 0; // 2-phase behavior that is probably closer to unstable than metastable
} // end function determineLBMetastableResults()


/**
 * @brief Function _sortMetastableLBResults_ reviews lattice Boltzmann simulation results for all metastable test points to classify 2-phase or 3-phase behavior.
 *
 * This is the controlling function that reads the results of LB simulations of test points in the 3-phase region that are expected to exhibit
 * metastable behavior.  The results are input to function _determineLBMetastableResults_ to be classified accordingly.
 */
void sortMetastableLBResults() {
	int readEOF = 0;
	double particlesA = 0.0, particlesB = 0.0;
	double rhoAMaxA = 0.0, rhoBMaxA = 0.0;
	double rhoAMaxB = 0.0, rhoBMaxB = 0.0;
	double rhoAMinN = 0.0, rhoBMinN = 0.0;
	FILE *phase_results = openFile("LB-metastable-region", ".dat", "r");
	FILE *metastable_2_results = openFile("LB-metastable-2-phase-region", ".dat", "w");
	FILE *metastable_3_results = openFile("LB-metastable-3-phase-region", ".dat", "w");

	if (setThreePhaseRegion) readThreePhaseRegionCoordinates(); // need 3-phase coords to define metastable points

	printf("\nSorting metastable LB results based on 2-phase or 3-phase behavior...\n");
	if (threePhaseRegionExists) {

		// Read lines from the metastable LB results file in a triplet: 1. maxA point, 2. maxB point, 3. minN point
		readEOF = fscanf(phase_results, "%17lf %17lf %8lf %8lf", &rhoAMaxA, &rhoBMaxA, &particlesA, &particlesB);
		if (readEOF != EOF) {
			readEOF = fscanf(phase_results, "%17lf %17lf %8lf %8lf", &rhoAMaxB, &rhoBMaxB, &particlesA, &particlesB);
			if (readEOF != EOF) readEOF = fscanf(phase_results, "%17lf %17lf %8lf %8lf", &rhoAMinN, &rhoBMinN, &particlesA, &particlesB); // this line has EOF
			else printf("\n\nError reading from metastable LB results file... may not have been written in the expected triple-line format.\n\n");
		}
		else printf("\n\nError reading from metastable LB results file... may not have been written in the expected triple-line format.\n\n");

		while (readEOF != EOF) {
			if (isEqual(rhoAMinN, 0.0, 1e-6) && isEqual(rhoBMinN, 0.0, 1e-6)) { // 3rd phase zeroed out if metastable point held 2 phases
				fprintf(metastable_2_results,"%f %f %f %f\n", rhoAMaxA, rhoBMaxA, particlesA, particlesB); // phase1Index read into rhoA
				fprintf(metastable_2_results,"%f %f %f %f\n", rhoAMaxB, rhoBMaxB, particlesA, particlesB); // phase2Index read into rhoB
			}
			else { // the metastable 2-phase initialization tried to form a 3rd phase
				fprintf(metastable_3_results,"%f %f %f %f\n", rhoAMaxA, rhoBMaxA, particlesA, particlesB);
				fprintf(metastable_3_results,"%f %f %f %f\n", rhoAMaxB, rhoBMaxB, particlesA, particlesB);
				fprintf(metastable_3_results,"%f %f %f %f\n", rhoAMinN, rhoBMinN, particlesA, particlesB);
			}

			// Read the next line triplet
			readEOF = fscanf(phase_results, "%17lf %17lf %8lf %8lf", &rhoAMaxA, &rhoBMaxA, &particlesA, &particlesB);
			if (readEOF != EOF) {
				readEOF = fscanf(phase_results, "%17lf %17lf %8lf %8lf", &rhoAMaxB, &rhoBMaxB, &particlesA, &particlesB);
				if (readEOF != EOF) readEOF = fscanf(phase_results, "%17lf %17lf %8lf %8lf", &rhoAMinN, &rhoBMinN, &particlesA, &particlesB); // this line has EOF
				else {
					printf("\n\nError reading from metastable LB results file... may not have been written in the expected triple-line format.\n\n");
					break;
				}
			}
			else {
				printf("\n\nEnd of file.\n\n");
				break;
			}
		} // end while
	}
	else printf("\nCannot classify metastable behavior without a defined 3-phase region.\n\n");

	fclose(phase_results);
	fclose(metastable_2_results);
	fclose(metastable_3_results);
	printf("Metastable LB results sorted.\n\n");
} // end function sortMetastableLBResults()


/**
 * @brief Function _getPhaseDiagramRhoAVsRhoBTwoPhaseRegion_ runs LB simulations in the 2-phase regions of a phase diagram.
 *
 * This function runs lattice Boltzmann simulations in the 2-phase regions of a phase diagram at test points defined by the global values
 * _nA0_, _nB0_ for A- and B-component total particles.  All test points are automatically looped through, and a test point is excluded if it
 * falls within a defined 3-phase region of the phase diagram.  The LB simulation points are along the lines defined as follows:
 * - The A-axis at [1,0] (the value 1e-12 is used instead of 0 to avoid discontinuities).
 * - The B-axis at [0,1] (again using 1e-12 instead of 0).
 * - Vertically from the A-axis (1,0) to the point [1,1] inclusive in steps of +0.1.
 * - Horizontally from the B-axis (0,1) to the point (1,1) exclusive in steps of +0.1.
 * - Diagonally from the B-axis (0,1) to the A-axis (1,0) exclusive in steps of (+0.1,-0.1).
 * - Diagonally from the point (1,1) exclusive to the point [1.4,1.4] inclusive in steps of (+0.1,+0.1).
 *
 * Test points on the A- and B-axes are initialized with a randomized density profile centered on the global values _nA0_ and _nB0_.
 * All other test points are initialized by default in a 2-phase step profile with densities given by the theoretical minimization.
 *
 * @note The loop for the last diagonal line of test points automatically adjusts the global values _kappa_ and _gammaMu_ to "soften"
 * the simulation at the interfaces and help keep the simulation stable.
 * @note The global parameter _useTwoPhaseStepInitialization_ can change the method used to initialize the LB simulations (1 = random,
 * 2 = 2-phase steps, 3 = 3-phase steps).
 */
void getPhaseDiagramRhoAVsRhoBTwoPhaseRegion() {
	int keepManualKappaFit = 0;
	double lbStepSize = 0.1;
	double particleThreshold = 1e-3;
	double tmpKappa = kappa;
	double evA, evB; // throw-away variables

	// Start benchmarking to compare the processor time used to the wall clock time taken to execute
	time_t clockRunTime;
	clock_t processorRunTime;
	clockRunTime = time(NULL);
	processorRunTime = clock();

	setInitializeRandom(); // quick random init to make sure the log file has the correct variables in the filename
	FILE *phase_results = openFile("LB-two-phase-region", ".dat", "w");
	if (useStepOrRandomInitialization) useTwoOrThreePhaseInitialization = 2;
	else useTwoOrThreePhaseInitialization = 1; // random profile

	// Set manual kappa selection for the duration of this series of LB sims
	determineKappaAutomatically ? (determineKappaAutomatically = 0, keepManualKappaFit = 0) : (keepManualKappaFit = 1);

	//
	// Hack the A-component axis
	//
	nA0 = 1.0;
	nB0 = 1e-12;
	if (!calculateMinimizationPath(nA0, nB0, &evA, &evB)) printf("(%f,%f) is not expected to phase separate... skipping.\n\n", nA0, nB0);
	else {
		printf("\nnA0=%f\tnB0=%f\n", nA0, nB0);
		setLBInitializationProfile();
		for (int i = 0; i < phase_iterations; i++) iteration();
		findIndicesOfMinMaxValues1DArray(n1, &phase1Index, &phase2Index);
		if (!isEqual(n1[phase1Index], nA0, particleThreshold) && !isEqual(n1[phase2Index], nA0, particleThreshold)) {
			fprintf(phase_results,"%f %f %f %f\n", n1[phase1Index], n2[phase1Index], nA0, nB0);
			fprintf(phase_results,"%f %f %f %f\n", n1[phase2Index], n2[phase2Index], nA0, nB0);
		}
	}

	//
	// Hack the B-component axis
	//
	nA0 = 1e-12;
	nB0 = 1.0;
	if (!calculateMinimizationPath(nA0, nB0, &evA, &evB)) printf("(%f,%f) is not expected to phase separate... skipping.\n\n", nA0, nB0);
	else {
		printf("\nnA0=%f\tnB0=%f\n", nA0, nB0);
		determineKappaAutomatically = 0; // reset in case the last test point didn't have a stable simulation
		kappa = tmpKappa; // reset kappa; instabilities in previous runs may alter it from the desired value
		setLBInitializationProfile();
		for (int i = 0; i < phase_iterations; i++) iteration();
		findIndicesOfMinMaxValues1DArray(n2, &phase1Index, &phase2Index);
		if (!isEqual(n2[phase1Index], nB0, particleThreshold) && !isEqual(n2[phase2Index], nB0, particleThreshold)) {
			fprintf(phase_results,"%f %f %f %f\n", n1[phase1Index], n2[phase1Index], nA0, nB0);
			fprintf(phase_results,"%f %f %f %f\n", n1[phase2Index], n2[phase2Index], nA0, nB0);
		}
	}

	setThreePhaseRegion = 1; // ensure the correct boundaries are being used

	//
	// Move vertically from A-component axis to [1,1]
	//
	nA0 = 1.0;
	for (nB0 = 0.1; nB0 < 1.0-particleThreshold; nB0 += lbStepSize) {

		// Only handle expected 2-phase region points here... separate function to do the 3-phase LB simulations
		if (insideThreePhaseRegion(nA0, nB0)) continue;
		if (!calculateMinimizationPath(nA0, nB0, &evA, &evB)) {
			printf("(%f,%f) is not expected to phase separate... skipping.\n", nA0, nB0);
			continue;
		}
		else printf("\nnA0=%f\tnB0=%f\n", nA0, nB0);

		determineKappaAutomatically = 0; // reset in case the last test point didn't have a stable simulation
		kappa = tmpKappa; // reset kappa; instabilities in previous runs may alter it from the desired value
		setLBInitializationProfile(); // do this every time in case no theoretical data forces random initialization but step is wanted
		for (int i = 0; i < phase_iterations; i++) iteration();
		if (!isEqual(n1[phase1Index], nA0, particleThreshold) && !isEqual(n2[phase1Index], nB0, particleThreshold)) {
			fprintf(phase_results,"%f %f %f %f\n", n1[phase1Index], n2[phase1Index], nA0, nB0);
			fprintf(phase_results,"%f %f %f %f\n", n1[phase2Index], n2[phase2Index], nA0, nB0);
		}
	}

	//
	// Move horizontally from B-component axis to (1.0,1.0)
	//
	nB0 = 1.0;
	for (nA0 = 0.1; nA0 < 1.0-particleThreshold; nA0 += lbStepSize) {

		// Only handle expected 2-phase region points here... separate function to do the 3-phase LB simulations
		if (insideThreePhaseRegion(nA0, nB0)) continue;
		if (!calculateMinimizationPath(nA0, nB0, &evA, &evB)) {
			printf("(%f,%f) is not expected to phase separate... skipping.\n", nA0, nB0);
			continue;
		}
		else printf("\nnA0=%f\tnB0=%f\n", nA0, nB0);

		determineKappaAutomatically = 0; // reset in case the last test point didn't have a stable simulation
		kappa = tmpKappa; // reset kappa; instabilities in previous runs may alter it from the desired value
		setLBInitializationProfile(); // do this every time in case no theoretical data forces random initialization but step is wanted
		for (int i = 0; i < phase_iterations; i++) iteration();
		if (!isEqual(n1[phase1Index], nA0, particleThreshold) && !isEqual(n2[phase1Index], nB0, particleThreshold)) {
			fprintf(phase_results,"%f %f %f %f\n", n1[phase1Index], n2[phase1Index], nA0, nB0);
			fprintf(phase_results,"%f %f %f %f\n", n1[phase2Index], n2[phase2Index], nA0, nB0);
		}
	}

	//
	// Move diagonally from B-axis (0.0,1.0) to A-axis (1.0,0.0)
	//
	nB0 = 1.0;
	for (nA0 = 0.1; nA0 < 1.0-particleThreshold; nA0 += lbStepSize) {
		nB0 -= lbStepSize; // doing it this way to keep the loop on track if an iteration is skipped below

		// Only handle expected 2-phase region points here... separate function to do the 3-phase LB simulations
		if (insideThreePhaseRegion(nA0, nB0)) continue;
		if (!calculateMinimizationPath(nA0, nB0, &evA, &evB)) {
			printf("(%f,%f) is not expected to phase separate... skipping.\n", nA0, nB0);
			continue;
		}
		else printf("\nnA0=%f\tnB0=%f\n", nA0, nB0);

		determineKappaAutomatically = 0; // reset in case the last test point didn't have a stable simulation
		kappa = tmpKappa; // reset kappa; instabilities in previous runs may alter it from the desired value
		setLBInitializationProfile(); // do this every time in case no theoretical data forces random initialization but step is wanted
		for (int i = 0; i < phase_iterations; i++) iteration();

		// This region may have maxA/maxB and minA/minB paired together instead of the usual maxA/minB, minA/maxB
		if (!isEqual(n1[phase1Index], nA0, particleThreshold) && !isEqual(n2[phase1Index], nB0, particleThreshold)) {
			fprintf(phase_results,"%f %f %f %f\n", n1[phase1Index], n2[phase1Index], nA0, nB0);
			fprintf(phase_results,"%f %f %f %f\n", n1[phase2Index], n2[phase2Index], nA0, nB0);
		}
	}

	//
	// Move diagonally from [1.0,1.0] to (1.5,1.5) - expect only 2-phase behavior in this region
	//
	nA0 = 0.9;
	nB0 = 0.9;
	int j = -1;
	while ( (nA0 <= 1.5-particleThreshold) && (nB0 < 1.5-particleThreshold) ) {
		nA0 += 0.1; // doing it this way to keep the loop on track if an iteration is skipped below
		nB0 += 0.1;
		j++;
		if (!calculateMinimizationPath(nA0, nB0, &evA, &evB)) {
			printf("(%f,%f) is not expected to phase separate... skipping.\n", nA0, nB0);
			continue;
		}
		else printf("\nnA0=%f\tnB0=%f\n", nA0, nB0);

		(kappa < 0.15 + j*0.0875) ? (kappa = 0.15 + j*0.0875) : 0; // set minimum kappa values for the binary liquid region; linearly increase from .15 to .51
		(tmpKappa > kappa) ? (kappa = tmpKappa) : 0; // now use larger of the minimum allowed kappa or the manually chosen one (also fixes prior instabilities)
		determineKappaAutomatically = 0; // reset in case the last test point didn't have a stable simulation
		setLBInitializationProfile(); // do this every time in case no theoretical data forces random initialization but step is wanted
		for (int i = 0; i < phase_iterations; i++) iteration();

		if (!isEqual(n1[phase1Index], nA0, particleThreshold) && !isEqual(n2[phase1Index], nB0, particleThreshold)) {
			fprintf(phase_results,"%f %f %f %f\n", n1[phase1Index], n2[phase1Index], nA0, nB0);
			fprintf(phase_results,"%f %f %f %f\n", n1[phase2Index], n2[phase2Index], nA0, nB0);
		}
	}

	clockRunTime = time(NULL) - clockRunTime;
	processorRunTime = clock() - processorRunTime;
	useTheoreticalDensities = 1;
	fclose(phase_results);

	// Reset kappa to automatic fit if started out that way; keep manual fit otherwise
	keepManualKappaFit ? (determineKappaAutomatically = 0) : (determineKappaAutomatically = 1);

	printRunTime(processorRunTime, clockRunTime);
	printf("\nDone! 2-phase region LB simulations complete.\n\n");

} // end function getPhaseDiagramRhoAVsRhoBTwoPhaseRegion()


/**
 * @brief Function _getPhaseDiagramRhoAVsRhoBThreePhaseRegion_ runs LB simulations in the 3-phase region of a phase diagram.
 *
 * This function exhaustively tests the 3-phase region of a phase diagram with lattice Boltzmann simulations.  Test points are defined by the
 * global values _nA0_, _nB0_ for A- and B-component total particles.  Tests include both "full" 3-phase points and 2-phase points that exist
 * in the 3-phase region and phase separate into metastable densities.
 * - Test points classified as "full" 3-phase points are initialized by default in a 3-phase step profile.
 * - Test points that are 2-phase metastable within the 3-phase region are initialized by default in a 2-phase step profile.
 *
 * @note The global boolean _useThreePhaseStepInitialization_ can toggle the initialization of the LB simulations between a step profile
 * (on) and one with random noise about the A- and B-component test point values (off).
 * @note The global boolean _checkOnlyTwoPhasePoints_ can be used to test only the 2-phase metastable points and exclude all "full" 3-phase
 * points from the simulations.
 */
void getPhaseDiagramRhoAVsRhoBThreePhaseRegion() {
	int readEOF = 0;
	int keepManualKappaFit = 0;
	int minNIndex = 0, maxAIndex = 0, maxBIndex = 0;
	double particlesATotal = 0.0, particlesBTotal = 0.0;
	double tmpKappa = kappa;
	FILE *phaseDiagram_densities_theoretical;
	FILE *phase_results;

	// Set manual kappa selection for the duration of this series of LB sims
	determineKappaAutomatically ? (determineKappaAutomatically = 0, keepManualKappaFit = 0) : (keepManualKappaFit = 1);

	// Start "benchmarking" to compare the processor time used to the wall clock time taken to execute
	time_t clockRunTime;
	clock_t processorRunTime;
	clockRunTime = time(NULL);
	processorRunTime = clock();

	setInitializeRandom(); // quick random init to make sure the log files have the correct variables in the filenames
	if (checkThreePhaseLBPoints) { // can skip LB sims of 3-phase points to save time
		printf("\nChecking points that minimized to 3 phases in the 3-phase region...\n\n");
		phaseDiagram_densities_theoretical = openFile("threePhaseDiagram-densities-threePhases", ".dat", "r");
		phase_results = openFile("LB-three-phase-region", ".dat", "w");

		readEOF = fscanf(phaseDiagram_densities_theoretical, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
				&particlesATotal, &particlesBTotal, &theoreticalDensityA1, &theoreticalDensityB1,
				&theoreticalDensityA2, &theoreticalDensityB2, &theoreticalDensityA3, &theoreticalDensityB3);
		while (readEOF != EOF) {
			nA0 = particlesATotal;
			nB0 = particlesBTotal;
			determineKappaAutomatically = 0; // reset in case the last test point didn't have a stable simulation
			kappa = tmpKappa; // reset kappa; instabilities in previous runs may alter it from the desired value
			if (useStepOrRandomInitialization) useTwoOrThreePhaseInitialization = 3;
			setLBInitializationProfile();

			printf("nA0=%f\tnB0=%f\n", nA0, nB0);
			for (int i = 0; i < phase_iterations; i++) iteration();

			// After phases settle, identify the min/max values across the lattice
			minNIndex = 0;
			maxAIndex = 0;
			maxBIndex = 0;
			for (int i = 0; i < XDIM; i++) {
				if (n1[i] > n1[maxAIndex]) maxAIndex = i;
				if (n2[i] > n2[maxBIndex]) maxBIndex = i;
				if (n[i] < n[minNIndex]) minNIndex = i;
			}

			fprintf(phase_results,"%.15f %.15f %f %f\n", n1[maxAIndex], n2[maxAIndex], nA0, nB0);
			fprintf(phase_results,"%.15f %.15f %f %f\n", n1[maxBIndex], n2[maxBIndex], nA0, nB0);
			fprintf(phase_results,"%.15f %.15f %f %f\n", n1[minNIndex], n2[minNIndex], nA0, nB0);
			readEOF = fscanf(phaseDiagram_densities_theoretical, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
					&particlesATotal, &particlesBTotal, &theoreticalDensityA1, &theoreticalDensityB1,
					&theoreticalDensityA2, &theoreticalDensityB2, &theoreticalDensityA3, &theoreticalDensityB3);
		} // end while

		fclose(phase_results);
		fclose(phaseDiagram_densities_theoretical);
	} // end if (checkThreePhaseLBPoints)

	//
	// Check the 2-phase points in case any are within the 3-phase region
	//

	if (checkTwoPhaseMetastableLBPoints) { // can skip 2-phase metastable LB sims to save time
		printf("\nChecking 2-phase points that may be metastable in the 3-phase region...\n\n");
		phaseDiagram_densities_theoretical = openFile("threePhaseDiagram-densities-twoPhases", ".dat", "r"); // recycle the same file handles as before
		phase_results = openFile("LB-metastable-region", ".dat", "w");

		setThreePhaseRegion = 1; // ensure the correct boundaries are being used
		readEOF = fscanf(phaseDiagram_densities_theoretical, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
				&particlesATotal, &particlesBTotal, &theoreticalDensityA1, &theoreticalDensityB1,
				&theoreticalDensityA2, &theoreticalDensityB2, &theoreticalDensityA3, &theoreticalDensityB3);
		while (readEOF != EOF) {
			if (insideThreePhaseRegion(particlesATotal, particlesBTotal)) {
				nA0 = particlesATotal;
				nB0 = particlesBTotal;
				determineKappaAutomatically = 0; // reset in case the last test point didn't have a stable simulation
				kappa = tmpKappa; // reset kappa; instabilities in previous runs may alter it from the desired value
				if (useStepOrRandomInitialization) useTwoOrThreePhaseInitialization = 2;
				setLBInitializationProfile();

				printf("nA0=%f\tnB0=%f\n", nA0, nB0);
				for (int i = 0; i < phase_iterations; i++) iteration();

				// After phases settle, identify the min/max values across the lattice
				minNIndex = 0;
				maxAIndex = 0;
				maxBIndex = 0;
				for (int i = 0; i < XDIM; i++) {
					if (n1[i] > n1[maxAIndex]) maxAIndex = i;
					if (n2[i] > n2[maxBIndex]) maxBIndex = i;
					if (n[i] < n[minNIndex]) minNIndex = i;
				}

				int metastableBehavior = determineLBMetastableResults(n1[maxAIndex], n2[maxAIndex], n1[maxBIndex], n2[maxBIndex], n1[minNIndex], n2[minNIndex]);
				printf("metastableBehavior = %i\n", metastableBehavior);
				if (metastableBehavior == -1) { // 3rd phase has emerged
					fprintf(phase_results,"%.15f %.15f %f %f\n", n1[maxAIndex], n2[maxAIndex], nA0, nB0);
					fprintf(phase_results,"%.15f %.15f %f %f\n", n1[maxBIndex], n2[maxBIndex], nA0, nB0);
					fprintf(phase_results,"%.15f %.15f %f %f\n", n1[minNIndex], n2[minNIndex], nA0, nB0);
				}
				else { // sample from the center of the domains if metastable point holds 2 phases
					fprintf(phase_results,"%.15f %.15f %f %f\n", n1[phase1Index], n2[phase1Index], nA0, nB0);
					fprintf(phase_results,"%.15f %.15f %f %f\n", n1[phase2Index], n2[phase2Index], nA0, nB0);
					fprintf(phase_results,"%.15f %.15f %f %f\n", 0.0, 0.0, nA0, nB0); // zero this out if metastable holds the 2 phases it was initialized with
				}
				readEOF = fscanf(phaseDiagram_densities_theoretical, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
						&particlesATotal, &particlesBTotal, &theoreticalDensityA1, &theoreticalDensityB1,
						&theoreticalDensityA2, &theoreticalDensityB2, &theoreticalDensityA3, &theoreticalDensityB3);
			}
			else {
				printf("Skipping %f, %f...\n", particlesATotal, particlesBTotal);
				readEOF = fscanf(phaseDiagram_densities_theoretical, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
						&particlesATotal, &particlesBTotal, &theoreticalDensityA1, &theoreticalDensityB1,
						&theoreticalDensityA2, &theoreticalDensityB2, &theoreticalDensityA3, &theoreticalDensityB3);
			}
		} // end while

		fclose(phase_results);
		fclose(phaseDiagram_densities_theoretical);
		sortMetastableLBResults();
	} // end if (checkTwoPhaseMetastableLBPoints)

	// Reset kappa to automatic fit if started out that way; keep manual fit otherwise
	keepManualKappaFit ? (determineKappaAutomatically = 0) : (determineKappaAutomatically = 1);

	clockRunTime = time(NULL) - clockRunTime;
	processorRunTime = clock() - processorRunTime;
	printRunTime(processorRunTime, clockRunTime);
	printf("\nDone! 3-phase region LB simulations complete.\n\n");
} // end function getPhaseDiagramRhoAVsRhoBThreePhaseRegion()


/**
 * @brief Function _generatePhaseDiagramRhoAVsRhoB_ controls the lattice Boltzmann simulations to generate a phase diagram in the 2- and 3-phase regions.
 *
 * This function is a wrapper function to control the lattice Boltzmann simulations in the 2- and 3-phase regions of a theoretical phase diagram.
 * It allows testing of only the 2-phase, only the 3-phase, or both phase diagram regions.  It provides the ability to select a single chemical potential
 * forcing method or all three for use in the simulation.  However, note that only the default forcing method (2) is allowed for testing the 3-phase
 * region.
 *
 * @note The global boolean _testAllForcingMethods_ controls if all chemical potential forcing methods or only the manually selected one are tested.
 * @note The global boolean _setPhaseDiagramTwoPhaseRegion_ controls if the 2-phase LB simulations are to be run.
 * @note The global boolean _setPhaseDiagramThreePhaseRegion_ controls if the 3-phase LB simulations are to be run.
 */
void generatePhaseDiagramRhoAVsRhoB() {
	double tmpKappa = kappa;

	// Generate LB phase diagram for the 2-phase region
	if (setPhaseDiagramTwoPhaseRegion) getPhaseDiagramRhoAVsRhoBTwoPhaseRegion();

	// Generate LB phase diagram for the 3-phase region (including metastable 2-phase points)
	kappa = tmpKappa; // reset kappa; the binary liquids or instabilities in the 2-phase LB sims may alter it from the desired value
	if (setPhaseDiagramThreePhaseRegion && useChemicalPotentialForcingMethod == 2) getPhaseDiagramRhoAVsRhoBThreePhaseRegion(); // only allow method 2 (other is horrible)
	else if (setPhaseDiagramThreePhaseRegion && useChemicalPotentialForcingMethod != 2) printf("\nNot testing 3-phase region with a non-optimal forcing method.\n\n");
} // end function generatePhaseDiagramRhoAVsRhoB()


/**
 * @brief Function _generateInterfaceWidthFitValues_ is a test driver used to verify the relationship between kappa and interface widths.
 *
 * This function loops over 6 density test pairs: (0.8,0.8), (1.0,0.05), (1.0,0.3), (1.0,0.55), (1.0,1.0), (1.4,1.4).  The pairs were chosen to provide
 * test coverage in the 2-phase liquid-vapor, 2-phase binary liquid, and 3-phase regions of most phase diagrams.  For each pair, the kappa parameter is
 * varied from 0.1 to 5.0 to generate a fit for the LB simulation interface width.  At each kappa value, the LB simulation is initialized and allowed to
 * iterate 50,000 iterations so a quasi-equilibrium interface width can be numerically measured (i.e. the normal initialization routine).  A factor is
 * derived by assuming the interface width is proportional to the square root of kappa (width = kappaFactor * sqrt(kappa)).  The kappa, interface width,
 * and kappaFactor (known as alpha in some literature) are logged, and an associated gnuplot script is written to .
 */
void generateInterfaceWidthFitValues() {
	int numberOfTestPoints = 6;
	int keepManualKappaFit = 0;
	double testAParticles[6] = {0.8, 1.0, 1.0, 1.0, 1.0, 1.4};
	double testBParticles[6] = {0.8, 0.05, 0.3, 0.55, 1.0, 1.4};
	double minKappaValue = 0.1;
	double maxKappaValue = 5.0;
	double kappaStep = 0.1;
	double tmpKappa = 0.1;
	double threshold = 1e-6;

	for (int i = 0; i < numberOfTestPoints; i++) {
		nA0 = testAParticles[i];
		nB0 = testBParticles[i];

		char width_fit_results_name[255];
		sprintf(width_fit_results_name, "interface-width-fit_nA0%f_nB0%f", nA0, nB0);
		FILE *width_fit_results = openFile(width_fit_results_name, ".dat", "w");

		for (kappa = minKappaValue; kappa <= maxKappaValue+threshold; kappa += kappaStep) {
			determineKappaAutomatically ? (determineKappaAutomatically = 0, keepManualKappaFit = 0) : (keepManualKappaFit = 1); // ensure manual kappa is used
			tmpKappa = kappa; // save current kappa to keep the loop going correctly if the gamma fit is unstable and an automatic kappa is used at any point
			kappaFactorDetermined = 0; // set this flag to force a fresh interface fit for each kappa
			setLBInitializationProfile(); // interface width fit happens in this function
			kappaFactor = interfaceWidth / sqrt(kappa);
			fprintf(width_fit_results, "%f %f %f %f %f\n", kappa, kappaFactor, interfaceWidthForGivenKappa, gammaMu, gammaFactor);
			kappa = tmpKappa; // restore the loop entry value of kappa to stay on track in case it was automatically overridden
		}

		fclose(width_fit_results);
		gnuplotInterfaceWidthFit();
	}
	keepManualKappaFit ?: (determineKappaAutomatically = 1);
} // end function generateInterfaceWidthFitValues()


/**
 * @brief Function _calculateKandSParameters_ uses the values of the van der Waals constants for each component in the mixture to calculate the values of
 * the global phase diagram parameters (zeta, Lambda) for a binary van der Waals mixture as defined by van Konynenburg and Scott.
 *
 * @see P. H. van Konynenburg and R. L. Scott, Philosophical Transactions of the Royal Society of London A: Mathematical, Physical and
 * Engineering Sciences 298, 495 (1980)
 */
void calculateKandSParameters() {
	double a1 = 0.0, a2 = 0.0, a12 = aAB;
	double b1 = 0.0, b2 = 0.0;
	if (tcB < tcA) {
		a1 = aB;
		a2 = aA;
		b1 = bB;
		b2 = bA;
	}
	else {
		a1 = aA;
		a2 = aB;
		b1 = bA;
		b2 = bB;
	}

	double xi = (b2 - b1) / (b1 + b2);
	double zeta = (a2/(b2*b2) - a1/(b1*b1)) / (a1/(b1*b1) + a2/(b2*b2));
	double bigLambda = (a1/(b1*b1) - 2.0*a12/(b1*b2) + a2/(b2*b2)) / (a1/(b1*b1) + a2/(b2*b2));

	printf("\nvan K and Scott parameters:\nxi = %f\nzeta = %f\nLambda = %f\n\n", xi, zeta, bigLambda);
} // end function calculateKandSParameters()


/**
 * @brief Function _printComponentDensities_ prints the density values for both A- and B-components to the console.
 */
void printComponentDensities() {
	printf("\n\n");
	for (int i = 0; i < XDIM; i++) printf("%f ", n1[i]);
	printf("\n");
	for (int i = 0; i < XDIM; i++) printf("%f ", n2[i]);
	printf("\n\n");
} // end function printComponentDensities()


/**
 * @brief Function _printComponentMaxMin_ prints the maximum and minimum density values for both A- and B-components to the console.
 */
void printComponentMaxMin() {
	double maxA = n1[0];
	double minA = n1[0];
	double maxB = n2[0];
	double minB = n2[0];

	for (int i = 0; i < XDIM; i++) {
		if (n1[i] > maxA) maxA = n1[i];
		if (n1[i] < minA) minA = n1[i];
		if (n2[i] > maxB) maxB = n2[i];
		if (n2[i] < minB) minB = n2[i];
	}

	printf("\n");
	printf("maxA=%f\tminA=%f\n", maxA, minA);
	printf("maxB=%f\tminB=%f\n", maxB, minB);
	printf("\n\n");
} // end function printComponentMaxMin()


/**
 * @brief Function _printLBResults_ prints the lattice Boltzmann simulation results for a given test point to the console.
 *
 * @note This function will also print out the minimization results for comparison.
 *
 * @param [in] phase The region of the phase diagram (2 or 3) in which the given test point is located.
 * @param [in] particlesA The number of A particles in the given test point.
 * @param [in] particlesB the number of B particles in the given test point.
 */
void printLBResults(int phase, double particlesA, double particlesB) {
	int readEOF = 0;
	int lbFound = 0;
	int writePhase = 1;
	double threshold = 1e-6;
	double particlesATotal = 0.0, particlesBTotal = 0.0;
	double rhoA1 = 0.0, rhoA2 = 0.0, rhoA3 = 0.0;
	double rhoB1 = 0.0, rhoB2 = 0.0, rhoB3 = 0.0;
	double tmpRhoA = 0.0, tmpRhoB = 0.0;
	char *LB_coords_name;

	// Decide which LB results file to look in
	if (phase == 3) LB_coords_name = "LB-three-phase-region";
	else if (phase == 2 && insideThreePhaseRegion(particlesA, particlesB)) LB_coords_name = "LB-metastable-region";
	else LB_coords_name = "LB-two-phase-region";
	FILE *LB_coords = openFile(LB_coords_name, ".dat", "r");

	readEOF = fscanf(LB_coords, "%17lf %17lf %8lf %8lf", &tmpRhoA, &tmpRhoB, &particlesATotal, &particlesBTotal);
	while (readEOF != EOF) {
		if ( isEqual(particlesA, particlesATotal, threshold) && isEqual(particlesB, particlesBTotal, threshold) ) {
			lbFound = 1;
			switch (writePhase) {
			case 1:
				rhoA1 = tmpRhoA;
				rhoB1 = tmpRhoB;
				writePhase = 2;
				break;
			case 2:
				rhoA2 = tmpRhoA;
				rhoB2 = tmpRhoB;
				writePhase = 3;
				break;
			case 3:
				rhoA3 = tmpRhoA;
				rhoB3 = tmpRhoB;
				writePhase = -1;
				break;
			default:
				break;
			}
		}
		readEOF = fscanf(LB_coords, "%17lf %17lf %8lf %8lf", &tmpRhoA, &tmpRhoB, &particlesATotal, &particlesBTotal);
	}

	if (lbFound) {
		printf("LB results found... %s\n", LB_coords_name);
		printf("particlesA=%.15f\tparticlesB=%.15f\n", particlesA, particlesB);
		printf("%.15f %.15f %.15f %.15f %.15f %.15f\n\n", rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);
		printf("Compared to theory...\n");
		printMinimizationResults(particlesA, particlesB);
		determineLBMetastableResults(rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3); // used here to just print the 2 test statistics
	}

	fclose(LB_coords);
} // end function printLBResults()


/**
 * @brief Function _printMinimizationParentPoint_ prints out the parent test point for a given density pair that resulted from a free energy minimization.
 *
 * @note The child point coordinates are given by the global parameters _childRhoA_ and _childRhoB_.
 */
void printMinimizationParentPoint() {
	int readEOF = 0;
	int parentFound = 0;
	double threshold = 1e-6;
	double particlesATotal = 0.0, particlesBTotal = 0.0;
	double rhoA1 = 0.0, rhoA2 = 0.0, rhoA3 = 0.0;
	double rhoB1 = 0.0, rhoB2 = 0.0, rhoB3 = 0.0;

	printf("Searching for the parent of the minimization point %f, %f...\n\n", childRhoA, childRhoB);

	FILE *two_phase_coords = openFile("threePhaseDiagram-densities-twoPhases", ".dat", "r");
	readEOF = fscanf(two_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
			&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
	while (readEOF != EOF) {
		if ( (isEqual(rhoA1, childRhoA, threshold) && isEqual(rhoB1, childRhoB, threshold)) ||
				(isEqual(rhoA2, childRhoA, threshold) && isEqual(rhoB2, childRhoB, threshold)) ||
				(isEqual(rhoA3, childRhoA, threshold) && isEqual(rhoB3, childRhoB, threshold)) ) {
			printf("\n2-phase point found...\n");
			printf("particlesATotal=%.15f\tparticlesBTotal=%.15f\n", particlesATotal, particlesBTotal);
			printf("%.15f %.15f %.15f %.15f %.15f %.15f\n\n", rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);
			parentFound = 1;
			printLBResults(2, particlesATotal, particlesBTotal);
			break;
		}
		else {
			readEOF = fscanf(two_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
					&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
		}
	}
	fclose(two_phase_coords);

	FILE *three_phase_coords = openFile("threePhaseDiagram-densities-threePhases", ".dat", "r");
	readEOF = fscanf(three_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
			&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
	while (!parentFound && readEOF != EOF) {
		if ( (isEqual(rhoA1, childRhoA, threshold) && isEqual(rhoB1, childRhoB, threshold)) ||
				(isEqual(rhoA2, childRhoA, threshold) && isEqual(rhoB2, childRhoB, threshold)) ||
				(isEqual(rhoA3, childRhoA, threshold) && isEqual(rhoB3, childRhoB, threshold)) ) {
			printf("\n3-phase point found...\n");
			printf("particlesATotal=%.15f\tparticlesBTotal=%.15f\n", particlesATotal, particlesBTotal);
			printf("%.15f %.15f %.15f %.15f %.15f %.15f\n\n", rhoA1, rhoB1, rhoA2, rhoB2, rhoA3, rhoB3);
			parentFound = 1;
			printLBResults(3, particlesATotal, particlesBTotal);
			break;
		}
		else {
			readEOF = fscanf(three_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf %17lf %17lf %*f %*f %*f",
					&particlesATotal, &particlesBTotal, &rhoA1, &rhoB1, &rhoA2, &rhoB2, &rhoA3, &rhoB3);
		}
	}
	fclose(three_phase_coords);
} // end function printMinimizationParentPoint()


/**
 * @brief Function _printLBParentPoint_ prints out the parent test point for a given density pair that resulted from a lattice Boltzmann simulation.
 *
 * @note The child point coordinates are given by the global parameters _childRhoA_ and _childRhoB_.
 */
void printLBParentPoint() {
	int readEOF = 0;
	int parentFound = 0;
	double threshold = 1e-6;
	double particlesATotal = 0.0, particlesBTotal = 0.0;
	double tmpRhoA = 0.0, tmpRhoB = 0.0;

	printf("Searching for the parent of the LB result point of (%f,%f)...\n\n", childRhoA, childRhoB);

	FILE *LB_coords = openFile("LB-two-phase-region", ".dat", "r");
	if (LB_coords != NULL) {
		readEOF = fscanf(LB_coords, "%17lf %17lf %8lf %8lf", &tmpRhoA, &tmpRhoB, &particlesATotal, &particlesBTotal);
		while (readEOF != EOF) {
			if ( isEqual(tmpRhoA, childRhoA, threshold) && isEqual(tmpRhoB, childRhoB, threshold) ) {
				parentFound = 1;
				printf("\nFound 2-phase parent of LB results point (%f,%f)...\n", childRhoA, childRhoB);
				printLBResults(2, particlesATotal, particlesBTotal);
				break;
			}
			else readEOF = fscanf(LB_coords, "%17lf %17lf %8lf %8lf", &tmpRhoA, &tmpRhoB, &particlesATotal, &particlesBTotal);
		}
		fclose(LB_coords);
	}
	else printf("No 2-phase data to search.\n");

	if (!parentFound) {
		printf("Not in 2-phase data...\n");
		LB_coords = openFile("LB-metastable-region", ".dat", "r");
		if (LB_coords != NULL) {
			readEOF = fscanf(LB_coords, "%17lf %17lf %8lf %8lf", &tmpRhoA, &tmpRhoB, &particlesATotal, &particlesBTotal);
			while (readEOF != EOF) {
				if ( isEqual(tmpRhoA, childRhoA, threshold) && isEqual(tmpRhoB, childRhoB, threshold) ) {
					parentFound = 1;
					printf("\nFound metastable parent of LB results point (%f,%f)...\n", childRhoA, childRhoB);
					printLBResults(2, particlesATotal, particlesBTotal);
					break;
				}
				else readEOF = fscanf(LB_coords, "%17lf %17lf %8lf %8lf", &tmpRhoA, &tmpRhoB, &particlesATotal, &particlesBTotal);
			}
			fclose(LB_coords);
		}
		else printf("No metastable data to search.\n");
	}

	if (!parentFound) {
		printf("Not in metastable data...\n");
		LB_coords = openFile("LB-three-phase-region", ".dat", "r");
		if (LB_coords != NULL) {
			readEOF = fscanf(LB_coords, "%17lf %17lf %8lf %8lf", &tmpRhoA, &tmpRhoB, &particlesATotal, &particlesBTotal);
			while (readEOF != EOF) {
				if ( isEqual(tmpRhoA, childRhoA, threshold) && isEqual(tmpRhoB, childRhoB, threshold) ) {
					printf("\nFound 3-phase parent of LB results point (%f,%f)...\n", childRhoA, childRhoB);
					printLBResults(3, particlesATotal, particlesBTotal);
					break;
				}
				else readEOF = fscanf(LB_coords, "%17lf %17lf %8lf %8lf", &tmpRhoA, &tmpRhoB, &particlesATotal, &particlesBTotal);
			}
			fclose(LB_coords);
		}
		else printf("No 3-phase data to search.\n");
	}

	if (!parentFound) printf("\nUnable to find the parent of LB result (%f,%f)... "
			"ensure data exists, >= 6 digits of precision for the child point, and that rounding errors are not an issue.\n\n", childRhoA, childRhoB);
} // end function printLBParentPoint()


//
// Commented out functions below were used in development and are saved for posterity
//


//double calculateFreeEnergy(double particlesA, double particlesB, double volume) {
//	double invalidFreeEnergy = 1000.* ncA;
//	double freeEnergy = invalidFreeEnergy;
//	double excludedVolume = particlesA*bA + particlesB*bB;
//	double threshold = 1e-14;
//
//	excludedVolume = particlesA*bA + particlesB*bB;
//	if (isEqual(particlesA, 0.0, threshold) && isEqual(particlesB, 0.0, threshold) && isEqual(volume, 0.0, threshold)) {
//		freeEnergy = 0;
//	}
//	else if (!isEqual(volume, 0.0, threshold) && volume > 0 && !isEqual(volume, excludedVolume, threshold) && volume > excludedVolume) {
//		if (isEqual(particlesA, 0.0, threshold) && !isEqual(particlesB, 0.0, threshold) && particlesB > 0) {
//			freeEnergy = particlesB*theta*log(particlesB/(volume-excludedVolume)) - aB*particlesB*particlesB/volume - theta*particlesB;
//		}
//		else if (isEqual(particlesB, 0.0, threshold) && !isEqual(particlesA, 0.0, threshold) && particlesA > 0) {
//			freeEnergy = particlesA*theta*log(particlesA/(volume-excludedVolume)) - aA*particlesA*particlesA/volume - theta*particlesA;
//		}
//		else if (!isEqual(particlesA, 0.0, threshold) && particlesA > 0 && !isEqual(particlesB, 0.0, threshold) && particlesB > 0) {
//			freeEnergy = particlesA*theta*log(particlesA/(volume-excludedVolume)) + particlesB*theta*log(particlesB/(volume-excludedVolume))
//			- aA*particlesA*particlesA/volume - 2*aAB*particlesB*particlesA/volume - aB*particlesB*particlesB/volume - theta*(particlesA+particlesB);
//		}
//	}
//
//	return freeEnergy;
//}
//
//
//double calculateFreeEnergyDensity(double rhoA, double rhoB) {
//	double invalidFreeEnergy = 1000. * ncA;
//	double freeEnergy = invalidFreeEnergy;
//	double excludedVolume = rhoA*bA + rhoB*bB;
//	double zeroThreshold = 1e-14;
//
//	if (isEqual(rhoA, 0.0, zeroThreshold)) {
//		freeEnergy = rhoB*theta*log(rhoB/(1.0-excludedVolume)) - aB*rhoB*rhoB - theta*rhoB;
//	}
//	else if (isEqual(rhoB, 0.0, zeroThreshold)) {
//		freeEnergy = rhoA*theta*log(rhoA/(1.0-excludedVolume)) - aA*rhoA*rhoA - theta*rhoA;
//	}
//	else {
//		freeEnergy = rhoA*theta*log(rhoA/(1.0-excludedVolume)) + rhoB*theta*log(rhoB/(1.0-excludedVolume)) - aA*rhoA*rhoA -
//				2*aAB*rhoB*rhoA - aB*rhoB*rhoB - theta*(rhoA+rhoB);
//	}
//
//	return freeEnergy;
//}
//
//
//void generateFreeEnergySlices() {
//	//	int testPoints = 2;
//	int readEOF = 0;
//	double particlesA = 0.0, particlesB = 0.0;
//	double rhoA = 0.0, rhoB = 0.0;
//	double A1 = 0.0, B1 = 0.0, A2 = 0.0, B2 = 0.0;
//	double sStep = 0.001;
//
//	printf("Calculating free energy tie line slices to generate a surface...\n");
//	printf("tcA=%f\ttcB=%f\tncB=%f\tVDW interaction=%f\n", tcA, tcB, ncB, vdwInteractionFactor);
//
//	FILE *twoPhase_densities = openFile("threePhaseDiagram-densities-twoPhases", ".dat", "r");
//	FILE *free_energy_curves = openFile("free-energy-curves", ".dat", "w");
//
//	readEOF = fscanf(twoPhase_densities, "%17lf %17lf %17lf %17lf %17lf %17lf %*17f %*17f %*17f %*17f %*17f", &particlesA, &particlesB, &A1, &B1, &A2, &B2);
//	while (readEOF != EOF) {
//		if (isEqual(particlesA, 1.0, 1e-6) || isEqual(particlesB, 1.0, 1e-6) || isEqual(particlesA, particlesB, 1e-6)) {
//			for (double s = -100.0*sStep; s < (1.0 + 100.0*sStep); s += sStep) { // add extra steps on either side of the actual minimum densities to visualize curvature better
//				rhoA = A1 + s*(A2-A1);
//				rhoB = B1 + s*(B2-B1);
//				fprintf(free_energy_curves, "%f %f %f %f %f %.15f\n", particlesA, particlesB, s, rhoA, rhoB, calculateFreeEnergyDensity(rhoA, rhoB));
//			}
//		}
//		readEOF = fscanf(twoPhase_densities, "%17lf %17lf %17lf %17lf %17lf %17lf %*17f %*17f %*17f %*17f %*17f", &particlesA, &particlesB, &A1, &B1, &A2, &B2);
//	}
//
//	fclose(free_energy_curves);
//	fclose(twoPhase_densities);
//
////	//	Old test code to generate individual slices for troubleshooting... modified it to the above to generate a series to create a surface
////	static double AT[] = {1.3, 1.2, 1.4, 1.0, 1.3, 1.4};
////	static double BT[] = {1.3, 1.2, 1.4, 1.3, 2.0, 1.6};
////	static double A1[] = {1.877294, 1.523002, 2.165767, 1.171777, 2.694563, 2.400318}; // phase 1 is A-rich; pulled these from minimization results
////	static double B1[] = {0.615274, 0.814340, 0.492974, 1.093126, 0.329706, 0.412494};
////	static double A2[] = {0.440905, 0.595267, 0.344743, 0.813309, 0.212135, 0.280339}; // phase 2 is B-rich
////	static double B2[] = {2.318969, 1.922043, 2.649918, 1.524835, 3.302956, 2.929182};
////
////	FILE *free_energy_curves;
////	char free_energy_curves_name[255];
////
////	for (int i = 0; i < testPoints; i++) {
////		sprintf(free_energy_curves_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/"
////				"free-energy-curves-test_%f-%f_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat", AT[i], BT[i], tcA, tcB, aA, aAB, aB);
////		free_energy_curves = fopen(free_energy_curves_name, "w");
////		for (double s = -100.0*sStep; s < (1.0 + 100.0*sStep); s += sStep) { // add extra steps on either side of the actual minimum densities to visualize curvature better
////			rhoA = A1[i] + s*(A2[i]-A1[i]);
////			rhoB = B1[i] + s*(B2[i]-B1[i]);
////			excludedVolume = rhoA*bA + rhoB*bB;
////			freeEnergyDensity = rhoA*theta*log(rhoA/(1.0-excludedVolume)) + rhoB*theta*log(rhoB/(1.0-excludedVolume)) -
////					aA*rhoA*rhoA - 2*aAB*rhoB*rhoA - aB*rhoB*rhoB - theta*(rhoA+rhoB);
////
////			fprintf(free_energy_curves, "%f %f %f %f\n", rhoA, rhoB, s, freeEnergyDensity);
////		}
////		fclose(free_energy_curves);
////	}
//	printf("\nFree energy density slices generated.\n");
//} // end function generateFreeEnergySlices()
//
//
//void generateFreeEnergyData(){
//	int readEOF = 0;
//	double particlesA = 0.0, particlesB= 0.0;
//	double freeEnergyUnminimized = 0.0, freeEnergyMinimized = 0.0;
//
//	printf("Calculating free energy data set...\n");
//	printf("tcA=%f\ttcB=%f\tncB=%f\tVDW interaction=%f\n", tcA, tcB, ncB, vdwInteractionFactor);
//
//	FILE *densities_tmp = openFile("densities-tmp", ".dat", "r");
//	FILE *free_energy_data = openFile("free-energy-differences", ".dat", "w");
//	readEOF = fscanf(densities_tmp, "%*1i %17lf %17lf %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f "
//			"%9lf %*9f %*9f %*9f", &particlesA, &particlesB, &freeEnergyMinimized);
//	while (readEOF != EOF) {
//		freeEnergyUnminimized = calculateFreeEnergyDensity(particlesA, particlesB);
//		fprintf(free_energy_data, "%f %f %.15f\n", particlesA, particlesB, freeEnergyUnminimized-freeEnergyMinimized);
//		readEOF = fscanf(densities_tmp, "%*1i %17lf %17lf %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f %*17f "
//				"%9lf %*9f %*9f %*9f", &particlesA, &particlesB, &freeEnergyMinimized);
//	} // end while
//	fclose(free_energy_data);
//	fclose(densities_tmp);
//
//	printf("\nDone! Free energy data calculated.\n");
//} // end function generateFreeEnergyData()
////
////
////
////
////void fitKappa(FILE *kappa_results) {
////	int decreaseKappa = 1;
////	double theoreticalADensityProfile[XDIM];
////	double theoreticalBDensityProfile[XDIM];
////	double aSSE = 0.0;
////	double bSSE = 0.0;
////	double sseThreshold = 1e-12;
////	double minSSE = 1000.0;
////	double minKappa = 1000.0;
////	double minKappaFactor = 1.0;
////	double decreaseKappaFactor = 1.0;
////	double increaseKappaFactor = -0.1;
////	double n1Phase1Min = 0.0;
////	double n2Phase1Min = 0.0;
////	double n1Phase2Min = 0.0;
////	double n2Phase2Min = 0.0;
////
////	minKappa = 1000.0;
////	minSSE = 1000.0;
////	kappaFactor = 1.0;
////	minKappaFactor = 1.0;
////	decreaseKappa = 1;
////	decreaseKappaFactor = 1.0;
////	increaseKappaFactor = -0.2;
////	while (increaseKappaFactor < -0.01) {
////		setLBInitializationProfile();
////		aSSE = 0.0;
////		bSSE = 0.0;
////		for (int i = 0; i < XDIM; i++) {
////			theoreticalADensityProfile[i] = n1[i];
////			theoreticalBDensityProfile[i] = n2[i];
////		}
////		for (int i = 0; i < phase_iterations; i++) {
////			iteration();
////		}
////
////		if (isnan(n1[0]) || isnan(n2[0])) { // check the first element of each component to make sure the simulation was stable
////			printf("simulation was unstable...\n");
////			aSSE = 1000.0;
////			bSSE = 1000.0;
////		}
////		else {
////			for (int i = 0; i < XDIM; i++) {
////				aSSE += pow(n1[i]-theoreticalADensityProfile[i], 2);
////				bSSE += pow(n2[i]-theoreticalBDensityProfile[i], 2);
////			}
////			printf("aSSE=%.15f bSSE=%.15f\n", aSSE, bSSE);
////		}
////
////		if (aSSE + bSSE + sseThreshold < minSSE) {
////			printf("replacing minSSE=%.15f with aSSE+bSSE=%.15f\n", minSSE, aSSE+bSSE);
////			minSSE = aSSE + bSSE;
////			minKappa = kappa;
////			minKappaFactor = kappaFactor;
////			n1Phase1Min = n1[phase1Index];
////			n2Phase1Min = n2[phase1Index];
////			n1Phase2Min = n1[phase2Index];
////			n2Phase2Min = n2[phase2Index];
////		}
////		else if (decreaseKappa) {
////			decreaseKappaFactor *= 0.5;
////		}
////		else {
////			increaseKappaFactor *= 0.5;
////		}
////		kappaFactor = minKappaFactor + (decreaseKappa ? decreaseKappaFactor : increaseKappaFactor);
////
////		if (isEqual(kappaFactor, minKappaFactor, 0.01) || decreaseKappaFactor < 0.1) {
////			decreaseKappa = 0;
////		}
////	} // end while
////
////	printf("%f %f %f %f %.15f %.15f %.15f %.15f\n\n", nA0, nB0, minKappa, minKappaFactor, n1Phase1Min, n2Phase1Min, n1Phase2Min, n2Phase2Min);
////	fprintf(kappa_results,"%f %f %f %f %.15f %.15f %.15f %.15f\n", nA0, nB0, minKappa, minKappaFactor, n1Phase1Min, n2Phase1Min, n1Phase2Min, n2Phase2Min);
////} // end function fitKappa()
////
////
////void generateKappaFitValues() {
////	double particleThreshold = 1e-6;
////	double lbStepSize = 0.01;
////	double evA;
////	double evB;
////
////	FILE *kappa_results;
////	char kappa_results_name[255];
////	sprintf(kappa_results_name, "/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/kappa-gamma-fit_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat",
////			tcA, tcB, aA, aAB, aB);
////	kappa_results = fopen(kappa_results_name,"w");
////
////	// Start benchmarking to compare the processor time used to the wall clock time taken to execute
////	time_t clockRunTime;
////	clock_t processorRunTime;
////	clockRunTime = time(NULL);
////	processorRunTime = clock();
////	setThreePhaseRegion = 1; // ensure the correct boundaries are being used
////
////	//
////	// Move vertically from A-component axis to [1,1]
////	//
////	nA0 = 1.0;
////	for (nB0 = 0.01; nB0 <= 1.0; nB0 += lbStepSize) {
////		printf("\nFitting best kappa value for (%f,%f) with initialized interface width=%f\n\n", nA0, nB0, interfaceWidth);
////
////		// Only use 2-phase region points below the 3-phase region, and ensure a - eigenvalue
////		if (insideThreePhaseRegion(nA0, nB0) || !calculateMinimizationPath(nA0, nB0, &evA, &evB)) {
////			break;
////		}
////		fitKappa(kappa_results);
////	}
////
//////	//
//////	// Move horizontally from B-component axis to [1,1]
//////	//
//////	nB0 = 1.0;
//////	for (nA0 = 0.01; nA0 < 1.0-particleThreshold; nA0 += lbStepSize) {
//////
//////		// Only use 2-phase region points and those expected to phase separate to generate kappa fit data
//////		if (insideThreePhaseRegion(nA0, nB0) || !calculateMinimizationPath(nA0, nB0, &evA, &evB)) {
//////			break;
//////		}
//////		fitKappa(kappa_results);
//////	}
////
////	//
////	// Move diagonally from B-axis (0.0,1.0) to A-axis (1.0,0.0)
////	//
////	nB0 = 0.99;
////	for (nA0 = 0.01; nA0 < 1.0-particleThreshold; nA0 += lbStepSize) {
////		printf("\nFitting best kappa value for (%f,%f) with initialized interface width=%f\n\n", nA0, nB0, interfaceWidth);
////
////		// Only use 2-phase region points below the 3-phase region, and ensure a - eigenvalue
////		if (insideThreePhaseRegion(nA0, nB0) || !calculateMinimizationPath(nA0, nB0, &evA, &evB)) {
////			nB0 -= lbStepSize;
////			continue;
////		}
////		fitKappa(kappa_results);
////		nB0 -= lbStepSize;
////	}
////
//////	//
//////	// Move diagonally from [1.01,1.01] to (1.5,1.5) - expect only 2-phase behavior in this region
//////	//
//////	nA0 = 1.01;
//////	nB0 = 1.01;
//////	while ( (nA0 <= 1.5-particleThreshold) && (nB0 < 1.5-particleThreshold) ) {
//////
//////		// Only use 2-phase region points above the 3-phase region, and ensure a - eigenvalue
//////		if (insideThreePhaseRegion(nA0, nB0) || !calculateMinimizationPath(nA0, nB0, &evA, &evB)) {
//////			nA0 += lbStepSize;
//////			nB0 += lbStepSize;
//////			continue;
//////		}
//////		fitKappa(kappa_results);
//////
//////		nA0 += lbStepSize;
//////		nB0 += lbStepSize;
//////	}
////
//////	//
//////	// Use (0.8,0.8) as the center point of the 3-phase region
//////	//
//////	nA0 = 0.8;
//////	nB0 = 0.8;
//////	useTwoPhaseStepInitialization = 3;
//////	fitKappa(kappa_results);
////
////	fclose(kappa_results);
////	clockRunTime = time(NULL) - clockRunTime;
////	processorRunTime = clock() - processorRunTime;
////	printRunTime(processorRunTime, clockRunTime);
////	printf("\nDone with kappa convergence data\n\n");
////} // end function generateKappaFitValues()
////
////
////void threePhaseRibbonTest() {
////	double lbStepSize = 0.001;
////	double particleThreshold = 1e-6;
////
////	FILE *phase_results;
////	char phase_results_name[255];
////	sprintf(phase_results_name, "/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/"
////			"test-3-phase-ribbon-point_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat", tcA, tcB, aA, aAB, aB);
////	phase_results = fopen(phase_results_name,"w");
////
////	// Start benchmarking to compare the processor time used to the wall clock time taken to execute
////	time_t clockRunTime;
////	clock_t processorRunTime;
////	clockRunTime = time(NULL);
////	processorRunTime = clock();
////
////	printf("Testing points across a suspected ribbon of 3-phase behavior...\n\n");
////
////	//
////	// Move horizontally from [0.8,0.67] to (0.82,0.67)
////	//
////	nA0 = 0.80;
////	nB0 = 0.67;
////	for (; nA0 < 0.82-particleThreshold; nA0 += lbStepSize) {
////		printf("\nnA0=%f\tnB0=%f\n", nA0, nB0);
////		determineKappaAutomatically = 0; // reset in case the last test point didn't have a stable simulation
////		setLBInitializationProfile(); // do this every time in case no theoretical data forces random initialization but step is wanted
////		for (int i=0; i < phase_iterations; i++) {
////			iteration();
////		}
////
////		findIndicesOfMinMaxValues1DArray(n1, &phase1Index, &phase2Index);
////		if (!isEqual(n1[phase1Index], nA0, particleThreshold) && !isEqual(n2[phase1Index], nB0, particleThreshold)) {
////			fprintf(phase_results,"%f %f %f %f\n", n1[phase1Index], n2[phase1Index], nA0, nB0);
////			fprintf(phase_results,"%f %f %f %f\n\n", n1[phase2Index], n2[phase2Index], nA0, nB0);
////		}
////		Events(1);
////		DrawGraphs();
////	}
////
////	//
////	// Move down from [0.82,0.67] to (0.82,0.66)
////	//
////	for (; nB0 > 0.66+particleThreshold; nB0 -= lbStepSize) {
////		printf("\nnA0=%f\tnB0=%f\n", nA0, nB0);
////		determineKappaAutomatically = 0; // reset in case the last test point didn't have a stable simulation
////		setLBInitializationProfile(); // do this every time in case no theoretical data forces random initialization but step is wanted
////		for (int i=0; i < phase_iterations; i++) {
////			iteration();
////		}
////
////		findIndicesOfMinMaxValues1DArray(n2, &phase1Index, &phase2Index);
////		if (!isEqual(n1[phase1Index], nA0, particleThreshold) && !isEqual(n2[phase1Index], nB0, particleThreshold)) {
////			fprintf(phase_results,"%f %f %f %f\n", n1[phase1Index], n2[phase1Index], nA0, nB0);
////			fprintf(phase_results,"%f %f %f %f\n\n", n1[phase2Index], n2[phase2Index], nA0, nB0);
////		}
////		Events(1);
////		DrawGraphs();
////	}
////
////	//
////	// Move horizontally from [0.82,0.66] to [0.84,0.66]
////	//
////	for (; nA0 < 0.84-particleThreshold; nA0 += lbStepSize) {
////		printf("\nnA0=%f\tnB0=%f\n", nA0, nB0);
////		determineKappaAutomatically = 0; // reset in case the last test point didn't have a stable simulation
////		setLBInitializationProfile(); // do this every time in case no theoretical data forces random initialization but step is wanted
////		for (int i=0; i < phase_iterations; i++) {
////			iteration();
////		}
////
////		findIndicesOfMinMaxValues1DArray(n1, &phase1Index, &phase2Index);
////		if (!isEqual(n1[phase1Index], nA0, particleThreshold) && !isEqual(n2[phase1Index], nB0, particleThreshold)) {
////			fprintf(phase_results,"%f %f %f %f\n", n1[phase1Index], n2[phase1Index], nA0, nB0);
////			fprintf(phase_results,"%f %f %f %f\n\n", n1[phase2Index], n2[phase2Index], nA0, nB0);
////		}
////		Events(1);
////		DrawGraphs();
////	}
////
////	clockRunTime = time(NULL) - clockRunTime;
////	processorRunTime = clock() - processorRunTime;
////	fclose(phase_results);
////	printRunTime(processorRunTime, clockRunTime);
////	printf("\nDone! 3-phase ribbon test complete.\n\n");
////} // end function threePhaseRibbonTest()
////
//
///**
// * @brief Function _calculatePhaseDiagramRhoAVsRhoBTwoPhasesTheoretical_ calculates a theoretical binary phase diagram that allows for 2 phases to co-exist.
// *
// * @note This function was used as an intermediate step only while developing the 3-phase algorithm.  It has not been maintained and is kept for posterity.
// */
//void calculatePhaseDiagramRhoAVsRhoBTwoPhasesTheoretical(){
//	int freeEnergyMinimumIndexA = 1, freeEnergyMinimumIndexB = 1, freeEnergyMinimumIndexV = 1;
//	double particlesStepSize = 0.05;
//	double freeEnergyThreshold = 1e-6;
//	double freeEnergyTrialArray[3][3][3];
//
//	FILE *phaseDiagram_densities_theoretical;
//	FILE *phaseDiagram_spinodal;
//	FILE *phaseDiagram_tieLines;
//	char phaseDiagram_densities_theoretical_name[255];
//	char phaseDiagram_spinodal_name[255];
//	char phaseDiagram_tieLines_name[255];
//	sprintf(phaseDiagram_densities_theoretical_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/"
//			"twoPhaseDiagram-densities-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat", tcA, tcB, aA, aAB, aB);
//	sprintf(phaseDiagram_spinodal_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/"
//			"twoPhaseDiagram-spinodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat", tcA, tcB, aA, aAB, aB);
//	sprintf(phaseDiagram_tieLines_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/"
//			"twoPhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat", tcA, tcB, aA, aAB, aB);
//	phaseDiagram_densities_theoretical = fopen(phaseDiagram_densities_theoretical_name, "w");
//	phaseDiagram_spinodal = fopen(phaseDiagram_spinodal_name, "w");
//	phaseDiagram_tieLines = fopen(phaseDiagram_tieLines_name, "w");
//
//	#ifdef DEBUG_MINIMIZATION_ON
//	double debugValue = 2.5;
//	double * debugVariable = &particlesBTotal;
//	FILE *debugData;
//	char * debugData_name = "/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/twoPhaseDiagram-debugData.dat";
//	debugData = fopen(debugData_name, "w");
//	#endif
//
//	for (double particlesATotal = 0; particlesATotal < (1/bA-freeEnergyThreshold); particlesATotal += particlesStepSize) {
//		for (double particlesBTotal = 0; (particlesBTotal+particlesATotal) < (1/bB-freeEnergyThreshold); particlesBTotal += particlesStepSize) {
//			double volumeTotal = 1.0;
//			double volume1 = 0.5 * volumeTotal;
//			double particlesA1 = particlesATotal * volume1;
//			double particlesB1 = particlesBTotal * volume1;
//			double trialStepSize=0.1;
//			int phaseChangeOccurred = 0;
//
//			// initialize the minimum free energy to be that of the current test point and see if the while loop converges elsewhere
//			double freeEnergyMinimum = calculateFreeEnergyMinimumTwoPhaseTrial(theta, particlesA1, particlesB1, volume1,
//					particlesATotal, particlesBTotal, volumeTotal);
//
//			#ifdef DEBUG_MINIMIZATION_ON
//			if ((*debugVariable < debugValue+freeEnergyThreshold) && (*debugVariable > debugValue-freeEnergyThreshold)) {
//				fprintf(debugData, "Starting free energy minimum = %.15f\t%.15f\t%.15f\t%f\t%f\t%f\t%f\n",
//						freeEnergyMinimum,
//						particlesATotal,
//						particlesBTotal,
//						volumeTotal,
//						particlesA1,
//						particlesB1,
//						volume1);
//				fprintf(debugData, "Total particles = %.15f\t(loop limit %.15f)\n\n", (particlesATotal+particlesBTotal), (1/b));
//			}
//			#endif
//
//			while (trialStepSize > freeEnergyThreshold) {
//				for (int i = 0; i < 3; i++) {
//					for (int j = 0; j < 3; j++) {
//						for (int k = 0; k < 3; k++){
//							freeEnergyTrialArray[i][j][k] = calculateFreeEnergyMinimumTwoPhaseTrial(theta,
//																									particlesA1+trialStepSize*(i-1),
//																									particlesB1+trialStepSize*(j-1),
//																									volume1+trialStepSize*(k-1),
//																									particlesATotal,
//																									particlesBTotal,
//																									volumeTotal);
//
//							#ifdef DEBUG_MINIMIZATION_ON
//							if ((*debugVariable < debugValue+freeEnergyThreshold) && (*debugVariable > debugValue-freeEnergyThreshold)) {
//								fprintf(debugData, "%.15f %.15f %.15f %.15f %i %i %i %.15f\n",
//										freeEnergyTrialArray[i][j][k],
//										particlesA1+trialStepSize*(i-1),
//										particlesB1+trialStepSize*(j-1),
//										volume1+trialStepSize*(k-1),
//										i,
//										j,
//										k,
//										trialStepSize);
//							}
//							#endif
//						}
//					}
//				}
//
//				#ifdef DEBUG_MINIMIZATION_ON
//				if ((*debugVariable < debugValue+freeEnergyThreshold) && (*debugVariable > debugValue-freeEnergyThreshold)) {
//					fprintf(debugData, "\n");
//				}
//				#endif
//
//				findIndexOfMinimumValue3DArray(freeEnergyTrialArray, &freeEnergyMinimumIndexA, &freeEnergyMinimumIndexB, &freeEnergyMinimumIndexV);
//
//				#ifdef DEBUG_MINIMIZATION_ON
//				if ((*debugVariable < debugValue+freeEnergyThreshold) && (*debugVariable > debugValue-freeEnergyThreshold)) {
//					fprintf(debugData, "Current free energy minimum = %.15f\tCompared to... %.15f\t%i %i %i\n\n",
//							freeEnergyMinimum,
//							freeEnergyTrialArray[freeEnergyMinimumIndexA][freeEnergyMinimumIndexB][freeEnergyMinimumIndexV],
//							freeEnergyMinimumIndexA,
//							freeEnergyMinimumIndexB,
//							freeEnergyMinimumIndexV);
//				}
//				#endif
//
//				if (freeEnergyTrialArray[freeEnergyMinimumIndexA][freeEnergyMinimumIndexB][freeEnergyMinimumIndexV] < freeEnergyMinimum) {
//					phaseChangeOccurred = 1;
//
//					#ifdef DEBUG_MINIMIZATION_ON
//					if ((*debugVariable < debugValue+freeEnergyThreshold) && (*debugVariable > debugValue-freeEnergyThreshold)) {
//						fprintf(debugData, "Phase change occurred... replacing %.15f with %.15f\t%i %i %i\n\n",
//								freeEnergyMinimum,
//								freeEnergyTrialArray[freeEnergyMinimumIndexA][freeEnergyMinimumIndexB][freeEnergyMinimumIndexV],
//								freeEnergyMinimumIndexA,
//								freeEnergyMinimumIndexB,
//								freeEnergyMinimumIndexV);
//					}
//					#endif
//
//					particlesA1 += (freeEnergyMinimumIndexA-1) * trialStepSize;
//					particlesB1 += (freeEnergyMinimumIndexB-1) * trialStepSize;
//					volume1 += (freeEnergyMinimumIndexV-1) * trialStepSize;
//					freeEnergyMinimum = freeEnergyTrialArray[freeEnergyMinimumIndexA][freeEnergyMinimumIndexB][freeEnergyMinimumIndexV];
//				}
//				else {
//					trialStepSize /= 2;
//				}
//			} // end while loop
//
//			if (phaseChangeOccurred) {
//				fprintf(phaseDiagram_densities_theoretical, "%e %e %e %e\n", particlesA1/volume1, particlesB1/volume1,
//						(particlesATotal-particlesA1)/(1-volume1), (particlesBTotal-particlesB1)/(1-volume1));
//				fprintf(phaseDiagram_spinodal, "%e %e\n", particlesATotal, particlesBTotal);
//				if (tieLineNeeded(particlesATotal, particlesBTotal)) {
//					fprintf(phaseDiagram_tieLines, "%e %e\n%e %e\n\n", particlesA1/volume1, particlesB1/volume1,
//							(particlesATotal-particlesA1)/(1-volume1), (particlesBTotal-particlesB1)/(1-volume1));
//				}
//			}
//
//		} // end for loop (B particles)
//	} // end for loop (A particles)
//
//	printf("Done!\n");
//	fclose(phaseDiagram_densities_theoretical);
//	fclose(phaseDiagram_spinodal);
//	fclose(phaseDiagram_tieLines);
//
//	#ifdef DEBUG_MINIMIZATION_ON
//	fclose(debugData);
//	#endif
//
//	gnuplotTwoComponentTwoPhase(theta);
//} // end function calculatePhaseDiagramRhoAVsRhoBTwoPhasesTheoretical()
