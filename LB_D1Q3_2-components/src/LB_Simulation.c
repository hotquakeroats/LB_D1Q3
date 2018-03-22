/**
 * @file LB_Simulation.c
 * @author Alexander Wagner (modified by Kent S. Ridl)
 * @date (last modified: 29 November 2017)
 *
 * The module _LB_Simulation.c_ is the main loop of the program. The execution starts here.
 */

#include "LB_D1Q3_2-components.h"


/**
 * @brief Function _main_ is the main loop of the program.  The execution starts here.
 *
 * This routine is responsible of executing the program and calling the Graphical User Interface (GUI). When the variable "done" is changed to a non-zero value
 * in the GUI the program ends. Otherwise the program loops through the iterations and displays the new graphics through the "DrawGraphs" routine. The GUI is
 * also able to change the program execution by calling functions, e.g. the init function in the Initializations sub-menu.  The program execution can be
 * accelerated by increasing the value of the "Repeat" variable, because less time is spent re-drawing the graphics.
 *
 * @note An additional control feature was added to automatically stop the simulation after a specified number of iterations.  The number of iterations is set
 * by the global parameter _phase_iterations_.  Note that if _phase_iterations_ is not an integer multiple of the "Repeat" variable, the program will automatically
 * stop with extra iterations.
 */
int main(int argc, char *argv[]){
	int newdata=1;
	int i;

	if (setDataDirectory()) {
		init();
		GUI();

		while (done==0){
			Events(newdata);
			DrawGraphs();

			// Modified to add the automatic stop after a specified number of phase_iterations
			if (next || !Pause || (run && (iterations < tmp_phase_iterations))){
				newdata=1;
				next=0;
				for (i=0;i<Repeat;i++){
					iteration();
				}
				if (iterations >= tmp_phase_iterations) {		// iterations is the total iteration count
					run = 0;									// tmp_phase_iterations is the next stopping point for iterations
					tmp_phase_iterations += phase_iterations;	// phase_iterations is the number of iterations to run before automatically stopping
				}
			}
			else sleep(1);
		}
	}

	return 0;

} // end function main()
