/*
 * LB_simulation.c
 *
 *  Created on: Jul 25, 2016
 *      Author: Kent Ridl
 */


#include "LB_D1Q3_3-components.h"


/**
This is the main loop of the program. The execution starts here.
This routine is responsible of executing the program and calling
the Graphical User Interface (GUI). When the variable "done" is changed
to a non-zero value in the GUI the program ends. Otherwise the
program loops through the iterations and displays the new graphics
through the "DrawGraphs" routine. The GUI is also able to change
the program execution by calling functions, e.g. the init function
in the Initializations sub-menu.
The program execution can be accelerated by increasing the value of
the "Repeat" variable, because less time is spent re-drawing the
graphics.
 */
int main(int argc, char *argv[]){
	int newdata=1;
	int i;

	init();
	GUI();

	total_phase_iterations = phase_iterations;

	while (done==0){
		Events(newdata);
		GetData();
		DrawGraphs();
		if (next || !Pause || (run && (iterations < total_phase_iterations))){
			newdata=1;
			next=0;
			for (i=0;i<Repeat;i++){
				iteration();
			}
			if (iterations == total_phase_iterations) {
				run = 0;
				total_phase_iterations += phase_iterations;
			}
		}
		else sleep(1);
	}

	return 0;

} // end function main()
