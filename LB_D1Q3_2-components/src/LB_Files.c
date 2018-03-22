/**
 * @file LB_Files.c
 * @author Kent S. Ridl
 * @date 26 November 2017
 *
 * The module _LB_Files.c_ contains functionality to control opening log files (e.g. minimization results) and the creation of files that are constant
 * (e.g. Gnuplot scripts).
 */

#include "LB_D1Q3_2-components.h"


/**
 * @brief Function _openFile_ is a wrapper for fopen used to prepend directory location and append simulation-specific information to a file name.
 *
 * This function appends the critical temperature for each component (tcA, tcB) and the van Der Waals energetic interaction parameters (aA, aAB, aB)
 * to the file name of a log file.
 *
 * @param [in] *fileName A string literal holding the base name for the log file.
 * @param [in] *extension A string literal holding the file extension to use.
 * @param [in] *attributes A string literal holding the attributes (e.g. "r", "w", "w+", etc.) with which to open the file.
 *
 * @return Pointer to the file (FILE *).
 *
 * @see fopen documentation
 */
FILE * openFile(char *fileName, char *extension, char *attributes) {
	int currentLength = 0;
	char file_name[MAXSTRLEN];
	currentLength = snprintf(file_name, MAXSTRLEN, "%s/", dataDirectory);
	currentLength += snprintf(file_name+currentLength, MAXSTRLEN-currentLength, "%s", fileName);
	currentLength += snprintf(file_name+currentLength, MAXSTRLEN-currentLength, "_tcA%f_tcB%f_aA%f_aAB%f_aB%f", tcA, tcB, aA, aAB, aB);
	currentLength += snprintf(file_name+currentLength, MAXSTRLEN-currentLength, "%s", extension);

	FILE * file_handle;
	file_handle = fopen(file_name, attributes);
	return file_handle;
} // end function openFile()


int setDataDirectory() {
	dataDirectory = getenv("LB_MULTI");
	if (dataDirectory == NULL) {
		printf("\n\n\"$LB_MULTI\" is not set as a directory to write log files... defaulting to \"$HOME\"\n");
		dataDirectory = getenv("HOME");
		if (dataDirectory == NULL) {
			printf("\n\n!!!!!!!!!!\n\n");
			printf("WARNING -- No directory is set to which log files can be written!\n\n");
			printf("Please ensure $HOME is set, or alternatively define an environment variable \"$LB_MULTI\"");
			printf("\n\n!!!!!!!!!!\n\n");
			printf("Program terminating...\n");
			return 0;
		}
	}
	printf("\nLog data will be written to directory: %s\n\n", dataDirectory);
	return 1;
} // end function setDataDirectory()


/**
 * @brief Function _gnuplotTwoComponentTwoPhase_ creates a gnuplot script to display phase diagrams with 2 components and 2 phases allowed.
 *
 * @see gnuplot documentation for additional details.
 */
void gnuplotTwoComponentTwoPhase() {
	FILE * gnuplotScript;
	char *gnuplotScript_name = "/home/clark/school/Lattice Boltzmann/Maxwell construction/multicomponent/plot-multicomponent-twoPhase-diagram.gnuplot";
	gnuplotScript = fopen(gnuplotScript_name, "w");

	fprintf(gnuplotScript, "set xrange [0:3]\n");
	fprintf(gnuplotScript, "set yrange [0:3]\n");
	fprintf(gnuplotScript, "set size square");
	fprintf(gnuplotScript, "set key box\n");
	fprintf(gnuplotScript, "plot \"twoPhaseDiagram-densities-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 2:3\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"twoPhaseDiagram-densities-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 4:5\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"twoPhaseDiagram-spinodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"twoPhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'red'\n", tcA, tcB, aA, aAB, aB);

	fclose(gnuplotScript);
} // end function gnuplotTwoComponentTwoPhase()


/**
 * @brief Function _gnuplotTwoComponentThreePhase_ creates gnuplot scripts to display free energy minimization results with 2 components and 3 phases allowed.
 *
 * This function creates 11 gnuplot scripts.  Each script is named according to the information displayed and appended by van der Waals parameters and
 * critical temperatures that characterize the components in the simulation.
 * - plot-multicomponent-density-threePhase-diagram: Plots full phase diagram and binodal lines, theoretical tie lines, and 2- and 3-phase LB simulation results.
 * - plot-threePhase: Plots phase diagram and binodal points for 3-phase region and 3-phase LB simulation results.
 * - plot-twoPhase: Plots phase diagram and binodal points for 2-phase region, theoretical tie lines, and 2-phase LB simulation results.
 * - plot-metastable-region: Plots phase diagram and binodal points LB simulation results for metastable points inside the 3-phase region of the phase diagram.
 * - plot-metastable-2: Plots phase diagram and binodal points LB simulation results for metastable points that exhibit 2-phase behavior.
 * - plot-metastable-3: Plots phase diagram and binodal points LB simulation results for metastable points that exhibit 3-phase behavior.
 * - plot-pressure-deviations: Plots theoretical pressure deviations among 3 phases of a free energy minimization to confirm equilibrium.
 * - plot-pressure-map: Plots pressure deviations as a heat map on the density-density phase diagram.
 * - plot-chemical-potential-deviations: Plots theoretical chemical potential deviations among 3 phases of a free energy minimization to confirm equilibrium.
 * - plot-chemical-potential-A-map: Plots chemical potential deviations for the A component as a heat map on the density-density phase diagram.
 * - plot-chemical-potential-B-map: Plots chemical potential deviations for the B component as a heat map on the density-density phase diagram.
 *
 * @see gnuplot documentation for additional details.
 */
void gnuplotTwoComponentThreePhase() {
	FILE * gnuplotScript;

	// Plot the entire phase diagram
	gnuplotScript = openFile("plot-multicomponent-threePhase-diagram", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set xlabel \"Component A Density\"\n");
	fprintf(gnuplotScript, "set ylabel \"Component B Density\"\n");
	fprintf(gnuplotScript, "set size square\n");
	fprintf(gnuplotScript, "set grid\n\n");
	fprintf(gnuplotScript, "plot \"threePhaseDiagram-densities-twoPhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 13 linecolor rgb 'light-grey' title \"Minimization 2-phase Test Point\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-densities-threePhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 13 linecolor rgb 'dark-grey' title \"Minimization 3-phase Test Point\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with points pointtype 7 linecolor rgb 'black' title \"3-phase region coordinates\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'brown' title \"Theoretical Tie Line\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-binodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' title \"Theoretical Binodal\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-two-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'blue' title \"LB 2-phase region (log)\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-two-phase-region-nid_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'turquoise' title \"LB 2-phase region (nid)\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-metastable-2-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'green' title \"LB metastable region (2-phase)\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-metastable-3-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'red' title \"LB metastable region (3-phase)\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-three-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'gold' title \"LB 3-phase region\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"singularity-line-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);

	// Plot only the 3-phase region
	gnuplotScript = openFile("plot-threePhase", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set xlabel \"Component A Density\"\n");
	fprintf(gnuplotScript, "set ylabel \"Component B Density\"\n");
	fprintf(gnuplotScript, "set size square\n");
	fprintf(gnuplotScript, "set grid\n");
	fprintf(gnuplotScript, "set key box\n");
	fprintf(gnuplotScript, "plot \"threePhaseDiagram-densities-threePhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 3 linecolor rgb 'grey' title \"3-phase test point\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-metastable-3-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3:4 with points pointtype 3 linecolor rgb 'dark-grey' title \"Metastable point (3-phase)\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-binodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' title \"Binodal\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'web-green' title \"3-phase region\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"unconditionally-unstable-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'blue' title \"Unconditionally unstable\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-densities-threePhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3:4 with points pointtype 1 linecolor rgb 'red' title \"Phase 1\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-densities-threePhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 5:6 with points pointtype 2 linecolor rgb 'green' title \"Phase 2\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-densities-threePhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 7:8 with points pointtype 25 linecolor rgb 'black' title \"Phase 3\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-metastable-3-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'red' title \"LB metastable region (3-phase)\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-three-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'gold' title \"LB 3-phase region\"\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);

	// Plot only the 2-phase region
	gnuplotScript = openFile("plot-twoPhase", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set xlabel \"Component A Density\"\n");
	fprintf(gnuplotScript, "set ylabel \"Component B Density\"\n");
	fprintf(gnuplotScript, "set size square\n");
	fprintf(gnuplotScript, "set grid\n");
	fprintf(gnuplotScript, "set key box\n");
	fprintf(gnuplotScript, "plot \"threePhaseDiagram-densities-twoPhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 3 linecolor rgb 'grey' title \"2-phase test point\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-metastable-2-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3:4 with points pointtype 3 linecolor rgb 'dark-grey' title \"Metastable point (2-phase)\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-metastable-3-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3:4 with points pointtype 3 linecolor rgb 'white' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'red' title \"Tie line\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-binodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' title \"Binodal\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'web-green' title \"3-phase region\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-densities-twoPhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3:4 with points pointtype 1 linecolor rgb 'red' title \"Phase 1\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-densities-twoPhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 5:6 with points pointtype 2 linecolor rgb 'black' title \"Phase 2\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-metastable-2-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'green' title \"LB metastable region (2-phase)\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-two-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'blue' title \"LB 2-phase region\"\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);

	// Plot the metastable region
	gnuplotScript = openFile("plot-metastable-region", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set xlabel \"Component A Density\"\n");
	fprintf(gnuplotScript, "set ylabel \"Component B Density\"\n");
	fprintf(gnuplotScript, "set size square\n");
	fprintf(gnuplotScript, "set grid\n\n");
	fprintf(gnuplotScript, "plot \"LB-metastable-2-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3:4 with points pointtype 13 linecolor rgb 'light-grey' title \"Metastable 2-phase Test Point\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-metastable-3-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3:4 with points pointtype 13 linecolor rgb 'dark-grey' title \"Metastable 3-phase Test Point\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with points pointtype 7 linecolor rgb 'black' title \"3-phase region coordinates\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'brown' title \"Theoretical Tie Line\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-binodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' title \"Theoretical Binodal\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-metastable-2-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'green' title \"LB metastable region (2-phase)\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-metastable-3-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'red' title \"LB metastable region (3-phase)\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"singularity-line-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);

	// Plot the metastable region that exhibits 2-phase (-ish) behavior
	gnuplotScript = openFile("plot-metastable-2", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set xlabel \"Component A Density\"\n");
	fprintf(gnuplotScript, "set ylabel \"Component B Density\"\n");
	fprintf(gnuplotScript, "set size square\n");
	fprintf(gnuplotScript, "set grid\n\n");
	fprintf(gnuplotScript, "plot \"LB-metastable-2-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3:4 with points pointtype 13 linecolor rgb 'light-grey' title \"Metastable 2-phase Test Point\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with points pointtype 7 linecolor rgb 'black' title \"3-phase region coordinates\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'brown' title \"Theoretical Tie Line\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-binodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' title \"Theoretical Binodal\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-metastable-2-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'green' title \"LB metastable region (2-phase)\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"singularity-line-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);

	// Plot the metastable region that exhibits 3-phase (-ish) behavior
	gnuplotScript = openFile("plot-metastable-3", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set xlabel \"Component A Density\"\n");
	fprintf(gnuplotScript, "set ylabel \"Component B Density\"\n");
	fprintf(gnuplotScript, "set size square\n");
	fprintf(gnuplotScript, "set grid\n\n");
	fprintf(gnuplotScript, "plot \"LB-metastable-3-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3:4 with points pointtype 13 linecolor rgb 'dark-grey' title \"Metastable 3-phase Test Point\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with points pointtype 7 linecolor rgb 'black' title \"3-phase region coordinates\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'brown' title \"Theoretical Tie Line\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-binodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' title \"Theoretical Binodal\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"LB-metastable-3-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'red' title \"LB metastable region (3-phase)\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"singularity-line-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);

	// Plot the pressure differences among the three phases
	gnuplotScript = openFile("plot-pressure-deviations", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set key box\n");
	fprintf(gnuplotScript, "plot \"pressure-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3 with lines title \"Pressure deviations\"\n\n", tcA, tcB, aA, aAB, aB);
	//	fprintf(gnuplotScript, "# Actual pressures for each phase\n");
	//	fprintf(gnuplotScript, "replot \"pressure-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 4 with lines title \"pressure1\"\n", tcA, tcB, aA, aAB, aB);
	//	fprintf(gnuplotScript, "replot \"pressure-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 5 with lines title \"pressure2\"\n", tcA, tcB, aA, aAB, aB);
	//	fprintf(gnuplotScript, "replot \"pressure-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 6 with lines title \"pressure3\"", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);

	gnuplotScript = openFile("plot-pressure-map", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set size square\n");
	fprintf(gnuplotScript, "set grid\n");
	fprintf(gnuplotScript, "stats \"pressure-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "set cbrange [STATS_max/100:STATS_max]\n");
	fprintf(gnuplotScript, "set logscale cb\n");
	fprintf(gnuplotScript, "set format cb \"%%.1le%%L\"\n\n");
	fprintf(gnuplotScript, "plot \"pressure-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2:3 with points pointtype 7 pointsize 0.5 palette title \"Pressure Deviations\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-binodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with points pointtype 7 linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"singularity-line-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);

	// Plot the chemical potential difference among the three phases for each component
	gnuplotScript = openFile("plot-chemical-potential-deviations", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set key box\n");
	fprintf(gnuplotScript, "plot \"chemical-potential-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3 with lines title \"A deviations\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"chemical-potential-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 4 with lines title \"B deviations\"\n\n", tcA, tcB, aA, aAB, aB);
	//	fprintf(gnuplotScript, "# Actual chemical potentials\n");
	//	fprintf(gnuplotScript, "replot \"chemical-potential-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 5 with lines title \"A1\"\n", tcA, tcB, aA, aAB, aB);
	//	fprintf(gnuplotScript, "replot \"chemical-potential-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 6 with lines title \"B1\"\n", tcA, tcB, aA, aAB, aB);
	//	fprintf(gnuplotScript, "replot \"chemical-potential-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 7 with lines title \"A2\"\n", tcA, tcB, aA, aAB, aB);
	//	fprintf(gnuplotScript, "replot \"chemical-potential-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 8 with lines title \"B2\"\n", tcA, tcB, aA, aAB, aB);
	//	fprintf(gnuplotScript, "replot \"chemical-potential-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 9 with lines title \"A3\"\n", tcA, tcB, aA, aAB, aB);
	//	fprintf(gnuplotScript, "replot \"chemical-potential-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 10 with lines title \"B3\"\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);

	gnuplotScript = openFile("plot-chemical-potential-A-map", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set size square\n");
	fprintf(gnuplotScript, "set grid\n");
	fprintf(gnuplotScript, "stats \"chemical-potential-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "set cbrange [STATS_max/1000:STATS_max]\n");
	fprintf(gnuplotScript, "set logscale cb\n");
	fprintf(gnuplotScript, "set format cb \"%%.1le%%L\"\n\n");
	fprintf(gnuplotScript, "plot \"chemical-potential-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2:3 with points pointtype 7 pointsize 0.5 palette title \"Chemical Potential Deviations - A\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-binodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with points pointtype 7 linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"singularity-line-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);

	gnuplotScript = openFile("plot-chemical-potential-B-map", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set size square\n");
	fprintf(gnuplotScript, "set grid\n");
	fprintf(gnuplotScript, "stats \"chemical-potential-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 4\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "set cbrange [STATS_max/1000:STATS_max]\n");
	fprintf(gnuplotScript, "set logscale cb\n");
	fprintf(gnuplotScript, "set format cb \"%%.1le%%L\"\n\n");
	fprintf(gnuplotScript, "plot \"chemical-potential-deviations_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2:4 with points pointtype 7 pointsize 0.5 palette title \"Chemical Potential Deviations - B\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-binodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with points pointtype 7 linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"singularity-line-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);

	printf("Gnuplot scripts created!\n");
} // end function gnuplotTwoComponentThreePhase()


/**
 * @brief Function gnuplotEigenstuffData creates gnuplot scripts to display stability analysis results for a 2-component system.
 *
 * This function creates 4 gnuplot scripts.  Each script is named according to the information displayed and appended by van der Waals parameters and
 * critical temperatures that characterize the components in the simulation.
 * - plot-determinant-map: Plots the 2x2 Hessian determinant as a heat map on the density-density phase diagram.
 * - plot-eigenvalue-1-map: Plots the eigenvalue associated with the first eigenvector as a heat map on the density-density phase diagram.
 * - plot-eigenvalue-2-map: Plots the eigenvalue associated with the second eigenvector as a heat map on the density-density phase diagram.
 * - plot-eigenvector-map: Plots the two eigenvectors of the 2x2 Hessian matrix.  Eigenvectors are normalized, scaled by magnitude, and color coded to show
 * stable (blue) or unstable (red) behavior.
 *
 * @see gnuplot documentation for additional details.
 */
void gnuplotEigenstuffData() {
	FILE * gnuplotScript;

	gnuplotScript = openFile("plot-determinant-map", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set datafile missing \"?\"\n");
	fprintf(gnuplotScript, "set size square\n");
	fprintf(gnuplotScript, "set grid\n");
	fprintf(gnuplotScript, "stats \"eigenstuff_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "set cbrange [-1*log(abs(STATS_min)+1):log(abs(STATS_max)+1)]\n");
	fprintf(gnuplotScript, "set palette defined (-1*log(abs(STATS_min)+1) \"orange\", -1e-6 \"red\", 0 \"white\", 1e-6 \"blue\", log(abs(STATS_max)+1) \"black\")\n\n");
	fprintf(gnuplotScript, "plot \"eigenstuff_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2:($3<0 ? -1*log(abs($3)+1) : log(abs($3)+1)) with points pointtype 7 pointsize 1.5 palette title \"Determinant\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-densities-twoPhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 13 linecolor rgb 'light-grey' title \"Minimization 2-phase Test Point\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-densities-threePhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 13 linecolor rgb 'dark-grey' title \"Minimization 3-phase Test Point\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-binodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3:4 with points pointtype 7 linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 5:6 with points pointtype 7 linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"singularity-line-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);

	gnuplotScript = openFile("plot-eigenvalue-1-map", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set datafile missing \"?\"\n");
	fprintf(gnuplotScript, "set size square\n");
	fprintf(gnuplotScript, "set grid\n");
	fprintf(gnuplotScript, "stats \"eigenstuff_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 4\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "set cbrange [-1*log(abs(STATS_min)+1):log(abs(STATS_max)+1)]\n");
	fprintf(gnuplotScript, "set palette defined (-1*log(abs(STATS_min)+1) \"orange\", -1e-6 \"red\", 0 \"white\", 1e-6 \"blue\", log(abs(STATS_max)+1) \"black\")\n\n");
	fprintf(gnuplotScript, "plot \"eigenstuff_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2:($4<0 ? -1*log(abs($4)+1) : log(abs($4)+1)) with points pointtype 7 pointsize 1.5 palette title \"Eigenvalue 1\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-densities-twoPhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 13 linecolor rgb 'light-grey' title \"Minimization 2-phase Test Point\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-densities-threePhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 13 linecolor rgb 'dark-grey' title \"Minimization 3-phase Test Point\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-binodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3:4 with points pointtype 7 linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 5:6 with points pointtype 7 linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"singularity-line-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);

	gnuplotScript = openFile("plot-eigenvalue-2-map", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set datafile missing \"?\"\n");
	fprintf(gnuplotScript, "set size square\n");
	fprintf(gnuplotScript, "set grid\n");
	fprintf(gnuplotScript, "stats \"eigenstuff_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 5\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "set cbrange [-1*log(abs(STATS_min)+1):log(abs(STATS_max)+1)]\n");
	fprintf(gnuplotScript, "set palette defined (-1*log(abs(STATS_min)+1) \"orange\", -1e-6 \"red\", 0 \"white\", 1e-6 \"blue\", log(abs(STATS_max)+1) \"black\")\n\n");
	fprintf(gnuplotScript, "plot \"eigenstuff_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2:($5<0 ? -1*log(abs($5)+1) : log(abs($5)+1)) with points pointtype 7 pointsize 1.5 palette title \"Eigenvalue 2\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-densities-twoPhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 13 linecolor rgb 'light-grey' title \"Minimization 2-phase Test Point\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-densities-threePhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 13 linecolor rgb 'dark-grey' title \"Minimization 3-phase Test Point\"\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"threePhaseDiagram-binodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3:4 with points pointtype 7 linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 5:6 with points pointtype 7 linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot \"singularity-line-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);

	// For some unknown reason, the "replot" syntax will crash gnuplot for some mixtures... oh well.
	gnuplotScript = openFile("plot-eigenvector-map", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set datafile missing \"?\"\n");
	fprintf(gnuplotScript, "set size square\n");
	fprintf(gnuplotScript, "set grid\n\n");
	fprintf(gnuplotScript, "scale = %f\n", 0.5*minimizationParticlesStepSize);
	fprintf(gnuplotScript, "vnorm(x,y) = sqrt(x**2 + y**2)\n\n");
	fprintf(gnuplotScript, "plot \"eigenstuff_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using ($4<=0 ? $1 : 1/0):2:(scale*$6/vnorm($6,1)):(scale/vnorm($6,1)) with vectors linecolor rgb \"red\" notitle, \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"eigenstuff_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using ($4>0 ? $1 : 1/0):2:(scale*$6/vnorm($6,1)):(scale/vnorm($6,1)) with vectors linecolor rgb \"blue\" notitle, \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"eigenstuff_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using ($5<=0 ? $1 : 1/0):2:(scale*$7/vnorm($7,1)):(scale/vnorm($7,1)) with vectors linecolor rgb \"red\" notitle, \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"eigenstuff_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using ($5>0 ? $1 : 1/0):2:(scale*$7/vnorm($7,1)):(scale/vnorm($7,1)) with vectors linecolor rgb \"blue\" notitle, \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"threePhaseDiagram-densities-twoPhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 13 linecolor rgb 'light-grey' title \"Minimization 2-phase Test Point\", \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"threePhaseDiagram-densities-threePhases_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 13 linecolor rgb 'dark-grey' title \"Minimization 3-phase Test Point\", \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"threePhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle, \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"threePhaseDiagram-binodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle, \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'black' notitle, \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 3:4 with points pointtype 7 linecolor rgb 'black' notitle, \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"three-phase-region-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 5:6 with points pointtype 7 linecolor rgb 'black' notitle, \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"singularity-line-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);
} // end function gnuplotEigenstuffData()


/**
 * @brief Function _gnuplotInterfaceWidthFit_ creates gnuplot scripts to analyze the relationship between interface width and kappa.
 *
 * This function creates one gnuplot script to plot interface widths and associated kappa correction factors against a range of kappa values.
 * Also plotted are curve fits to see if the proportionality between interface width and sqrt(kappa) is satisfied.
 *
 * @see gnuplot documentation for additional details.
 */
void gnuplotInterfaceWidthFit() {
	char gnuplotScript_name[255];
	sprintf(gnuplotScript_name, "plot-interface-width-fit_nA0%f_nB0%f", nA0, nB0);
	FILE *gnuplotScript = openFile(gnuplotScript_name, ".gnuplot", "w");

	fprintf(gnuplotScript, "set xlabel \"Kappa\"\n");
	fprintf(gnuplotScript, "set key box\n");
	fprintf(gnuplotScript, "plot \"interface-width-fit_nA0%f_nB0%f_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:3 with points pointtype 7 linecolor rgb 'blue' title \"Interface width\", \\\n", nA0, nB0, tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"interface-width-fit_nA0%f_nB0%f_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 13 linecolor rgb 'red' title \"Kappa factor\"\n\n", nA0, nB0, tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "a1 = 1\n");
	fprintf(gnuplotScript, "b1 = 1\n");
	fprintf(gnuplotScript, "f(x) = a1*sqrt(x) + b1\n");
	fprintf(gnuplotScript, "fit f(x) \"interface-width-fit_nA0%f_nB0%f_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:3 via a1,b1\n", nA0, nB0, tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot f(x) linecolor rgb 'blue' title \"Interface width sqrt(x) fit\"\n\n");
	fprintf(gnuplotScript, "a3 = 1\n");
	fprintf(gnuplotScript, "h(x) = a3*sqrt(x)\n");
	fprintf(gnuplotScript, "fit h(x) \"interface-width-fit_nA0%f_nB0%f_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:3 via a3\n", nA0, nB0, tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot h(x) linecolor rgb 'green' title \"Interface width sqrt(x) fit - w/out intercept\"\n\n");
	fprintf(gnuplotScript, "a2 = 1\n");
	fprintf(gnuplotScript, "b2 = 1\n");
	fprintf(gnuplotScript, "g(x) = a2*x + b2\n");
	fprintf(gnuplotScript, "fit g(x) \"interface-width-fit_nA0%f_nB0%f_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 via a2,b2\n", nA0, nB0, tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, "replot g(x) linecolor rgb 'red' title \"Kappa factor linear fit\"\n\n");
	fclose(gnuplotScript);
} // end function gnuplotInterfaceWidthFit()


/**
 * @brief Function _gnuplotPublicationGraphics_ creates alternate versions of selected Gnuplot scripts that are more suitable for publication.
 *
 * The Gnuplot scripts that are created by other functions are meant primarily for analysis purposes.  Particular to the graphs of entire phase diagrams with both minimization and
 * lattice Boltzmann results, the colors, symbols, and line types are not optimal for publication (i.e. if a paper is printed in black and white).  This function creates a Gnuplot
 * script that accounts for such concerns.
 *
 * It also prints alternate free energy minimization data sets that emphasize ("filter") the unconditionally unstable 3-phase region of the phase diagram; all other points outside
 * of the unconditionally unstable 3-phase region are depicted the same.  There is no distinction made at this level to highlight metastable points, even for points that start in the
 * metastable region where the minimization resulted in 3-phase behavior.
 */
void gnuplotPublicationGraphics() {

	// Plot the entire phase diagram
	FILE *gnuplotScript = openFile("plot-multicomponent-threePhase-filtered", ".gnuplot", "w");
	fprintf(gnuplotScript, "set datafile commentschars \"#!\"\n");
	fprintf(gnuplotScript, "set xlabel \"Component A Density\"\n");
	fprintf(gnuplotScript, "set ylabel \"Component B Density\"\n");
	fprintf(gnuplotScript, "set size square\n");
	fprintf(gnuplotScript, "plot \"minimization-unconditionally-2-phase_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 13 linecolor rgb 'light-grey' title \"Two-phase/Metastable Test Point\", \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"minimization-metastable-densities_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 13 linecolor rgb 'light-grey' notitle, \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"minimization-unconditionally-3-phase_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 13 linecolor rgb 'dark-grey' title \"Three-phase Test Point\", \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"threePhaseDiagram-binodal_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' title \"Theoretical Binodal\", \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"threePhaseDiagram-tie-lines-theoretical_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linetype 0 linewidth 2 linecolor rgb 'brown' title \"Theoretical Tie Line\", \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"singularity-line-coords_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" with lines linecolor rgb 'black' notitle, \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"LB-two-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 7 linecolor rgb 'blue' title \"LB 2-phase region (log)\", \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"LB-two-phase-region-nid_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 9 linecolor rgb 'light-blue' title \"LB 2-phase region (nid)\", \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"LB-metastable-2-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 every 15 with points pointtype 10 linecolor rgb 'web-green' title \"LB metastable region (2-phase)\", \\\n", tcA, tcB, aA, aAB, aB);
	fprintf(gnuplotScript, " \"LB-three-phase-region_tcA%f_tcB%f_aA%f_aAB%f_aB%f.dat\" using 1:2 with points pointtype 5 linecolor rgb 'black' title \"LB 3-phase region\", \\\n", tcA, tcB, aA, aAB, aB);
	fclose(gnuplotScript);
} // end function gnuplotPublicationGraphics()


/**
 * @brief Function _logVDWParameters_ logs the van der Waals parameters that were used to create a phase diagram.
 */
void logVDWParameters() {
	FILE * vdwParametersLog = openFile("vdw-parameters-log", ".dat", "w");
	fprintf(vdwParametersLog, "Degrees of freedom:\n");
	fprintf(vdwParametersLog, "tcA = %f\ntcB = %f\nncB = %f\nVDW interaction factor = %f\n\n", tcA, tcB, ncB, vdwInteractionFactor);
	fprintf(vdwParametersLog, "Determined values:\n");
	fprintf(vdwParametersLog, "ncA = %f\npcA = %f\npcB = %f\naA = %f\naAB = %f\naB = %f\nbA=%f\nbB=%f\n\n", ncA, pcA, pcB, aA, aAB, aB, bA, bB);
	fclose(vdwParametersLog);
} // end function logVDWParameters()


/**
 * @brief Function _writeThreePhaseRegionCoordinates_ writes the coordinates of the points that bound the 3-phase region of a phase diagram.
 *
 * @note The coordinates for the first point are repeated to allow gnuplot to plot a closed polygon.
 */
void writeThreePhaseRegionCoordinates() {
	FILE * three_phase_coords = openFile("three-phase-region-coords", ".dat", "w");
	fprintf(three_phase_coords, "%.15f %.15f\n%.15f %.15f\n%.15f %.15f\n%.15f %.15f\n", // repeat the 1st point so gnuplot can close a polygon if desired
			threePhaseRegion[0], threePhaseRegion[1], threePhaseRegion[2], threePhaseRegion[3], threePhaseRegion[4], threePhaseRegion[5], threePhaseRegion[0], threePhaseRegion[1]);
	fclose(three_phase_coords);
	setThreePhaseRegion = 0;
} // end function writeThreePhaseRegionCoordinates()


/**
 * @brief Function _readThreePhaseRegionCoordinates_ reads the coordinates defining the 3-phase region of a phase diagram into a global array.
 *
 * The coordinates defining the 3-phase region of a phase diagram are stored in the global array _threePhaseRegion[]_ for runtime use.  When a
 * function needs this information but a phase diagram has not been created, this function will read in the coordinates that were logged from
 * a previous phase diagram creation and store them in _threePhaseRegion[]_.
 */
void readThreePhaseRegionCoordinates() {
	double zeroThreshold = 1e-14;

	FILE * three_phase_coords = openFile("three-phase-region-coords", ".dat", "r");
	fscanf(three_phase_coords, "%17lf %17lf %17lf %17lf %17lf %17lf", &threePhaseRegion[0], &threePhaseRegion[1], &threePhaseRegion[2], &threePhaseRegion[3], &threePhaseRegion[4], &threePhaseRegion[5]);

	if (isEqual(threePhaseRegion[0], 0.0, zeroThreshold) && isEqual(threePhaseRegion[1], 0.0, zeroThreshold) &&
			isEqual(threePhaseRegion[2], 0.0, zeroThreshold) && isEqual(threePhaseRegion[3], 0.0, zeroThreshold) &&
			isEqual(threePhaseRegion[4], 0.0, zeroThreshold) && isEqual(threePhaseRegion[5], 0.0, zeroThreshold)) {
		threePhaseRegionExists = 0;
	}
	else threePhaseRegionExists = 1;

	fclose(three_phase_coords);
	setThreePhaseRegion = 0; // coordinates no longer need to be set
	printf("\nRead three-phase region coordinates: (%f,%f) (%f,%f) (%f,%f)\n\n", threePhaseRegion[0], threePhaseRegion[1], threePhaseRegion[2], threePhaseRegion[3], threePhaseRegion[4], threePhaseRegion[5]);
} // end function readThreePhaseRegionCoordinates()


/**
 * @brief Function _getDensityProfile_ logs the density profiles for each component of a LB simulation for later plotting.
 */
void getDensityProfile() {
	FILE *density_profile;
	char density_profile_name[255];

	sprintf(density_profile_name,"density-profile-A_%i-iterations_nA0%f_nB0%f", iterations, nA0, nB0);
	density_profile = openFile(density_profile_name, ".dat", "w");
	for (int i = 0; i < XDIM; i++) fprintf(density_profile, "%.15f\n", n1[i]);
	fclose(density_profile);

	sprintf(density_profile_name,"density-profile-B_%i-iterations_nA0%f_nB0%f", iterations, nA0, nB0);
	density_profile = openFile(density_profile_name, ".dat", "w");
	for (int i = 0; i < XDIM; i++) fprintf(density_profile, "%.15f\n", n2[i]);
	fclose(density_profile);

	printf("Density profiles saved...\n\n");
} // end function getDensityProfile()


/**
 * @brief Function _getPressureProfile_ logs the pressure profile of a LB simulation for later plotting.
 * It also saves the associated theoretical pressure.
 */
void getPressureProfile() {
	FILE *pressure_profile;
	char pressure_profile_name[255];

	sprintf(pressure_profile_name,"pressure-profile_nA0%f_nB0%f", nA0, nB0);
	pressure_profile = openFile(pressure_profile_name, ".dat", "w");
	for (int i = 0; i < XDIM; i++) fprintf(pressure_profile, "%.15f\n", pressure[i]);
	fclose(pressure_profile);

	sprintf(pressure_profile_name,"pressure-profile-theoretical_nA0%f_nB0%f", nA0, nB0);
	pressure_profile = openFile(pressure_profile_name, ".dat", "w");
	for (int i = 0; i < XDIM; i++) fprintf(pressure_profile, "%.15f\n", theoreticalPressureArray[i]);
	fclose(pressure_profile);

	printf("Pressure profile saved...\n\n");
} // end function getPressureProfile()


/**
 * @brief Function _getChemicalPotentialProfile_ logs the chemical potential profiles for each component of a LB simulation for later plotting.
 * It also saves the associated theoretical chemical potentials.
 */
void getChemicalPotentialProfile() {
	FILE *chemical_potential_profile;
	char chemical_potential_profile_name[255];

	sprintf(chemical_potential_profile_name,"chemical-potential-profile-A_nA0%f_nB0%f", nA0, nB0);
	chemical_potential_profile = openFile(chemical_potential_profile_name, ".dat", "w");
	for (int i = 0; i < XDIM; i++) fprintf(chemical_potential_profile, "%.15f\n", mu1[i]);
	fclose(chemical_potential_profile);

	sprintf(chemical_potential_profile_name,"chemical-potential-profile-B_nA0%f_nB0%f", nA0, nB0);
	chemical_potential_profile = openFile(chemical_potential_profile_name, ".dat", "w");
	for (int i = 0; i < XDIM; i++) fprintf(chemical_potential_profile, "%.15f\n", mu2[i]);
	fclose(chemical_potential_profile);

	sprintf(chemical_potential_profile_name,"chemical-potential-profile-theoretical-A_nA0%f_nB0%f", nA0, nB0);
	chemical_potential_profile = openFile(chemical_potential_profile_name, ".dat", "w");
	for (int i = 0; i < XDIM; i++) fprintf(chemical_potential_profile, "%.15f\n", theoreticalMuAArray[i]);
	fclose(chemical_potential_profile);

	sprintf(chemical_potential_profile_name,"chemical-potential-profile-theoretical-B_nA0%f_nB0%f", nA0, nB0);
	chemical_potential_profile = openFile(chemical_potential_profile_name, ".dat", "w");
	for (int i = 0; i < XDIM; i++) fprintf(chemical_potential_profile, "%.15f\n", theoreticalMuBArray[i]);
	fclose(chemical_potential_profile);

	printf("Chemical potential profiles saved...\n\n");
} // end function getChemicalPotentialProfile()
