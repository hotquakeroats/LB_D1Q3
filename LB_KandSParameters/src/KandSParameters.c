
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "mygraph.h"


int next=0,Pause=1,done=0,Repeat=10,iterations;

double tcA = 0.4, tcB = 0.4, tcRatio = 1.0;
double ncA = 1.0, ncB = 1.0;
double pcA = 1.0, pcB = 1.0;

double aA = 0.1;
double aB = 0.1;
double aAB = 0.1;
double bA = 1./3.;
double bB = 1./3.;
double vdwInteractionFactor = 0.5;

double xi = 1.0;
double zeta = 1.0;
double bigLambda = 1.0;


void calculateVDWParameters() {
	pcA = 3.*tcA/8.; // determine pcA
	ncA = pcA / ((3./8.)*tcA); // fix ncA to be 1
	pcB = (3./8.)*tcB*ncB; // determine pcB

	aA = (27./64.)*(tcA*tcA/pcA);
	aB = (27./64.)*(tcB*tcB/pcB);
	aAB = sqrt(aA*aB) * vdwInteractionFactor;
	bA = tcA/(8.*pcA);
	bB = tcB/(8.*pcB);
}


void calculateKandSParameters() {
	calculateVDWParameters();

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

	tcRatio = tcA / tcB;
	xi = (b2 - b1) / (b1 + b2);
	zeta = (a2/(b2*b2) - a1/(b1*b1)) / (a1/(b1*b1) + a2/(b2*b2));
	bigLambda = (a1/(b1*b1) - 2.0*a12/(b1*b2) + a2/(b2*b2)) / (a1/(b1*b1) + a2/(b2*b2));
} // end function calculateKandSParameters()


void calculateMixtureParameters() {
	tcRatio = -ncB * (zeta + 1) / (zeta - 1);
	tcA = tcB * tcRatio;

	double nuSquared = (bigLambda-1)*(bigLambda-1) * (2.0 - (zeta+1)/(zeta-1) - (zeta-1)/(zeta+1)) / (4.0*ncB);
	vdwInteractionFactor = sqrt(nuSquared);

	calculateVDWParameters();
}


/**
This routine initializes the Grapical User Interface (GUI). First we
define the fields we want to visualize, and in the second part we
structure the menu where we can manipulate variables during the
simulation.
*/
void GUI() {
  StartMenu("K & S Parameters", 1);

	DefineFunction("calculateKandSParameters", &calculateKandSParameters);
	DefineDouble("* tc-A *", &tcA);
	DefineDouble("** tc-B **", &tcB);
	DefineDouble("* nc-B *", &ncB);
	DefineDouble("* VDW mix *", &vdwInteractionFactor);
	DefineDouble("tc ratio", &tcRatio);
	DefineDouble("nc-A", &ncA);
	DefineDouble("pc-A", &pcA);
	DefineDouble("pc-B", &pcB);
	DefineDouble("aA", &aA);
	DefineDouble("aB", &aB);
	DefineDouble("aAB", &aAB);
	DefineDouble("bA", &bA);
	DefineDouble("bB", &bB);
	DefineDouble("xi", &xi);
	DefineDouble("** zeta **", &zeta);
	DefineDouble("** Lambda **", &bigLambda);
	DefineFunction("calculateMixtureParameters", &calculateMixtureParameters);

    DefineBool("EXIT",&done);

  EndMenu();
  /* End of graphics menu */

} // end function GUI()


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

	calculateKandSParameters();
	GUI();

	while (done==0){
		Events(newdata);
		sleep(1);
	}

	return 0;
} // end function main()
