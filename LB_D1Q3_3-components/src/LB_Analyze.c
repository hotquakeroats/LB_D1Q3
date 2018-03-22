/*
 * LB_Analyze.c
 *
 *  Created on: Jul 7, 2016
 *      Author: clark
 */


#include "LB_D1Q3_3-components.h"


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
	int ipp = i + 2;
	int ip = i + 1;
	int im = i - 1;
	int imm = i - 2;

	// Gradient for a lattice with periodic BCs
	if (useBoundaryConditionsPeriodic) {
		if (ip == XDIM) {
			ip = 0;
		}
		else if (im == -1) {
			im = XDIM - 1;
		}
		result = 1./2. * (var[ip]-var[im]);
	}
	// Gradient for a lattice with walls at the boundaries
	// Adjust for existence of the walls with forward/backward differencing
	else {
		if (ip == XDIM) { // backward difference
			result = 1./2. * (var[imm] - 4.*var[im] + 3.*var[i]);
		}
		else if (im == -1) { // forward difference
			result = 1./2. * (-3.*var[i] + 4.*var[ip] - var[ipp]);
		}
		else { // central difference... normal case
			result = 1./2. * (var[ip]-var[im]);
		}
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
	int ippp = i + 3;
	int ipp = i + 2;
	int ip = i + 1;
	int im = i - 1;
	int imm = i - 2;
	int immm = i - 3;

	// Gradient for a lattice with periodic BCs
	// Instead of a halo, re-set the out-of-bounds index to the opposite end of the lattice
	if (useBoundaryConditionsPeriodic) {
		if (ip == XDIM) {
			ip = 0;
		}
		else if (im == -1) {
			im = XDIM - 1;
		}
		result = var[ip] - 2.*var[i] + var[im];
	}
	// Gradient for a lattice with walls at the boundaries
	// Adjust for existence of the walls with forward/backward differencing
	else {
		if (ip == XDIM) { // backward difference
			result = -1.*var[immm] + 4.*var[imm] -5.*var[im] + 2.*var[i];
		}
		else if (im == -1) { // forward difference
			result = 2.*var[i] - 5.*var[ip] + 4.*var[ipp] - var[ippp];
		}
		else { // central difference... normal case
			result = var[ip] - 2.*var[i] + var[im];
		}
	}

	return result;
} // end function laplace()


double filterArrayData(double *var, int i, int arrayNeighbors) {
	double result = 0;
	int windowWidth = arrayNeighbors*2+1;
	int index = i;

	for (index = i-arrayNeighbors; index <= i + arrayNeighbors; index++) {
		if (index < 0) {
			result += var[index+XDIM];
		}
		else if (index > XDIM-1) {
			result += var[index-XDIM];
		}
		else {
			result += var[index];
		}
	}

	return (result / windowWidth);
}


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
	int hours = 0;
	int minutes = 0;
	int seconds = 0;

	processorRunTime /= CLOCKS_PER_SEC;

	hours = processorRunTime / 3600;
	processorRunTime = processorRunTime % 3600;
	minutes = processorRunTime / 60;
	processorRunTime = processorRunTime % 60;
	seconds = processorRunTime;

	printf("Processor time... %.2d:%.2d:%.2d\n", hours, minutes, seconds);

	hours = clockRunTime / 3600;
	clockRunTime = (int) clockRunTime % 3600;
	minutes = clockRunTime / 60;
	clockRunTime = (int) clockRunTime % 60;
	seconds = clockRunTime;

	printf("Clock time... %.2d:%.2d:%.2d\n", hours, minutes, seconds);
} // end function printRunTime()


/**
 * @brief Function _logVDWParameters_ logs the van der Waals parameters that were used to create a phase diagram.
 */
void logVDWParameters() {
	FILE * vdwParametersLog;
	char vdwParametersLog_name[255];
	sprintf(vdwParametersLog_name, "/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/vdw-parameters-log_n0%f,%f,%f_tc%f,%f,%f", nA0, nB0, nC0, tcA, tcB, tcC);
	vdwParametersLog = fopen(vdwParametersLog_name, "w");
	fprintf(vdwParametersLog, "Total material:\n");
	fprintf(vdwParametersLog, "nA0 = %f\nnB0 = %f\nnC0 = %f\n\n", nA0, nB0, nC0);
	fprintf(vdwParametersLog, "Degrees of freedom:\n");
	fprintf(vdwParametersLog, "tcA = %f\ntcB = %f\ntcC = %f\nncB = %f\nncC = %f\nVDW interaction factor = %f\n\n", tcA, tcB, tcC, ncB, ncC, vdwInteractionFactor);
	fprintf(vdwParametersLog, "Determined values:\n");
	fprintf(vdwParametersLog, "ncA = %f\npcA = %f\npcB = %f\npcC = %f\naA = %f\naB = %f\naC = %f\naAB = %f\naAC = %f\naBC = %f\n\n", ncA, pcA, pcB, pcC, aA, aB, aC, aAB, aAC, aBC);
	fprintf(vdwParametersLog, "Force/Interface coefficients:\n");
	fprintf(vdwParametersLog, "gamma-Mu = %f\nlambda = %f\nkappa = %f\n\n", gammaMu, lambda, kappa);
	fclose(vdwParametersLog);
} // end function logVDWParameters()


/**
 * @brief Function _getDensityProfile_ logs the density profiles for each component of a LB simulation for later plotting.
 */
void getDensityProfile() {
	int i = 0;

	FILE *density_profile;
	char density_profile_name[255];

	sprintf(density_profile_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/density-profile-A_n0%f,%f,%f_tc%f,%f,%f", nA0, nB0, nC0, tcA, tcB, tcC);
	density_profile = fopen(density_profile_name,"w");
	for (i = 0; i < XDIM; i++) {
		fprintf(density_profile, "%.15f\n", nA[i]);
	}
	fclose(density_profile);

	sprintf(density_profile_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/density-profile-B_n0%f,%f,%f_tc%f,%f,%f", nA0, nB0, nC0, tcA, tcB, tcC);
	density_profile = fopen(density_profile_name,"w");
	for (i = 0; i < XDIM; i++) {
		fprintf(density_profile, "%.15f\n", nB[i]);
	}
	fclose(density_profile);

	sprintf(density_profile_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/density-profile-C_n0%f,%f,%f_tc%f,%f,%f", nA0, nB0, nC0, tcA, tcB, tcC);
	density_profile = fopen(density_profile_name,"w");
	for (i = 0; i < XDIM; i++) {
		fprintf(density_profile, "%.15f\n", nC[i]);
	}
	fclose(density_profile);

	printf("Density profiles saved...\n\n");
} // end function getDensityProfile()


/**
 * @brief Function _getPressureProfile_ logs the pressure profile of a LB simulation for later plotting.
 * It also saves the associated theoretical pressure.
 */
void getPressureProfile() {
	int i = 0;

	FILE *pressure_profile;
	char pressure_profile_name[255];

	sprintf(pressure_profile_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/pressure-profile_n0%f,%f,%f_tc%f,%f,%f", nA0, nB0, nC0, tcA, tcB, tcC);
	pressure_profile = fopen(pressure_profile_name,"w");
	for (i = 0; i < XDIM; i++) {
		fprintf(pressure_profile, "%.15f\n", pressure[i]);
	}
	fclose(pressure_profile);

	sprintf(pressure_profile_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/pressure-filtered_n0%f,%f,%f_tc%f,%f,%f", nA0, nB0, nC0, tcA, tcB, tcC);
	pressure_profile = fopen(pressure_profile_name,"w");
	for (i = 0; i < XDIM; i++) {
		fprintf(pressure_profile, "%.15f\n", pressureFiltered[i]);
	}
	fclose(pressure_profile);

	//	sprintf(pressure_profile_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/pressure-profile-theoretical_nA0%f_nB0%f.dat", nA0, nB0);
	//	pressure_profile = fopen(pressure_profile_name,"w");
	//	for (i = 0; i < XDIM; i++) {
	//		fprintf(pressure_profile, "%.15f\n", theoreticalPressureArray[i]);
	//	}
	//	fclose(pressure_profile);

	printf("Pressure profile saved...\n\n");
} // end function getDensityProfile()


/**
 * @brief Function _getChemicalPotentialProfile_ logs the chemical potential profiles for each component of a LB simulation for later plotting.
 * It also saves the associated theoretical chemical potentials.
 */
void getChemicalPotentialProfile() {
	int i = 0;

	FILE *chemical_potential_profile;
	char chemical_potential_profile_name[255];

	sprintf(chemical_potential_profile_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/chemical-potential-profile-A_n0%f,%f,%f_tc%f,%f,%f", nA0, nB0, nC0, tcA, tcB, tcC);
	chemical_potential_profile = fopen(chemical_potential_profile_name,"w");
	for (i = 0; i < XDIM; i++) {
		fprintf(chemical_potential_profile, "%.15f\n", muA[i]);
	}
	fclose(chemical_potential_profile);

	sprintf(chemical_potential_profile_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/chemical-potential-profile-B_n0%f,%f,%f_tc%f,%f,%f", nA0, nB0, nC0, tcA, tcB, tcC);
	chemical_potential_profile = fopen(chemical_potential_profile_name,"w");
	for (i = 0; i < XDIM; i++) {
		fprintf(chemical_potential_profile, "%.15f\n", muB[i]);
	}
	fclose(chemical_potential_profile);

	sprintf(chemical_potential_profile_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/chemical-potential-profile-C_n0%f,%f,%f_tc%f,%f,%f", nA0, nB0, nC0, tcA, tcB, tcC);
	chemical_potential_profile = fopen(chemical_potential_profile_name,"w");
	for (i = 0; i < XDIM; i++) {
		fprintf(chemical_potential_profile, "%.15f\n", muC[i]);
	}
	fclose(chemical_potential_profile);

	sprintf(chemical_potential_profile_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/chemical-potential-filtered-A_n0%f,%f,%f_tc%f,%f,%f", nA0, nB0, nC0, tcA, tcB, tcC);
	chemical_potential_profile = fopen(chemical_potential_profile_name,"w");
	for (i = 0; i < XDIM; i++) {
		fprintf(chemical_potential_profile, "%.15f\n", muAFiltered[i]);
	}
	fclose(chemical_potential_profile);

	sprintf(chemical_potential_profile_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/chemical-potential-filtered-B_n0%f,%f,%f_tc%f,%f,%f", nA0, nB0, nC0, tcA, tcB, tcC);
	chemical_potential_profile = fopen(chemical_potential_profile_name,"w");
	for (i = 0; i < XDIM; i++) {
		fprintf(chemical_potential_profile, "%.15f\n", muBFiltered[i]);
	}
	fclose(chemical_potential_profile);

	sprintf(chemical_potential_profile_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/chemical-potential-filtered-C_n0%f,%f,%f_tc%f,%f,%f", nA0, nB0, nC0, tcA, tcB, tcC);
	chemical_potential_profile = fopen(chemical_potential_profile_name,"w");
	for (i = 0; i < XDIM; i++) {
		fprintf(chemical_potential_profile, "%.15f\n", muCFiltered[i]);
	}
	fclose(chemical_potential_profile);

	//	sprintf(chemical_potential_profile_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/chemical-potential-profile-theoretical-A_nA0%f_nB0%f.dat", nA0, nB0);
	//	chemical_potential_profile = fopen(chemical_potential_profile_name,"w");
	//	for (i = 0; i < XDIM; i++) {
	//		fprintf(chemical_potential_profile, "%.15f\n", theoreticalMuAArray[i]);
	//	}
	//	fclose(chemical_potential_profile);
	//
	//	sprintf(chemical_potential_profile_name,"/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/chemical-potential-profile-theoretical-B_nA0%f_nB0%f.dat", nA0, nB0);
	//	chemical_potential_profile = fopen(chemical_potential_profile_name,"w");
	//	for (i = 0; i < XDIM; i++) {
	//		fprintf(chemical_potential_profile, "%.15f\n", theoreticalMuBArray[i]);
	//	}
	//	fclose(chemical_potential_profile);

	printf("Chemical potential profiles saved...\n\n");
} // end function getDensityProfile()


/**
 * @brief Function _printComponentDensities_ prints the density values for both A- and B-components to the console.
 */
void printComponentDensities() {
	int i = 0;

	printf("\n\n");

	for (i = 0; i < XDIM; i++) {
		printf("%f ", nA[i]);
	}

	printf("\n");

	for (i = 0; i < XDIM; i++) {
		printf("%f ", nB[i]);
	}

	printf("\n");

	for (i = 0; i < XDIM; i++) {
		printf("%f ", nC[i]);
	}

	printf("\n\n");
} // end function printComponentDensities()


/**
 * @brief Function _printComponentMaxMin_ prints the maximum and minimum density values for both A- and B-components to the console.
 */
void printComponentMaxMin() {
	int i = 0;
	double maxA = nA[0];
	double minA = nA[0];
	double maxB = nB[0];
	double minB = nB[0];
	double maxC = nC[0];
	double minC = nC[0];

	for (i = 0; i < XDIM; i++) {
		if (nA[i] > maxA) {
			maxA = nA[i];
		}
		if (nA[i] < minA) {
			minA = nA[i];
		}
		if (nB[i] > maxB) {
			maxB = nB[i];
		}
		if (nB[i] < minB) {
			minB = nB[i];
		}
		if (nC[i] > maxC) {
			maxC = nC[i];
		}
		if (nC[i] < minC) {
			minC = nC[i];
		}
	}

	printf("\n");
	printf("maxA=%f\tminA=%f\tratio=%f\n", maxA, minA, (maxA/minA));
	printf("maxB=%f\tminB=%f\tratio=%f\n", maxB, minB, (maxB/minB));
	printf("maxC=%f\tminC=%f\tratio=%f\n", maxC, minC, (maxC/minC));
	printf("\n\n");
} // end function printComponentMaxMin()


/**
 * @brief Function _printChemicalPotentials_ prints the chemical potentials across the lattice for all components to the console.
 */
void printChemicalPotentials() {
	int i = 0;

	printf("\n\n");

	for (i = 0; i < XDIM; i++) {
		printf("%f ", muA[i]);
	}

	printf("\n");

	for (i = 0; i < XDIM; i++) {
		printf("%f ", muB[i]);
	}

	printf("\n");

	for (i = 0; i < XDIM; i++) {
		printf("%f ", muC[i]);
	}

	printf("\n\n");
} // end function printChemicalPotentials()


/**
 * @brief Function _printVolumeExclusion_ prints the excluded volume across the lattice to the console.
 */
void printVolumeExclusion() {
	int i = 0;
	int maxIndex = 0;
	double maxExclusion = volumeExclusion[0];

	printf("\n\nVolume exclusion:\n");

	for (i = 0; i < XDIM; i++) {
		printf("%f ", volumeExclusion[i]);
		if (volumeExclusion[i] > maxExclusion) {
			maxExclusion = volumeExclusion[i];
			maxIndex = i;
		}
	}

	printf("\nmaxExclusion (%i): %f \n\n", maxIndex, maxExclusion);
} // end function printVolumeExclusion()
