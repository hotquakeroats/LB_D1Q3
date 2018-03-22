/*
 * LB_collisions.c
 *
 *  Created on: Aug 5, 2016
 *      Author: clark
 */

#include "LB_D1Q3_3-components.h"


void setCollisionForcingNewChemicalPotentialGradient() {
	collision = collisionForcingNewChemicalPotentialGradient;
	printf("Switching to gradMu forcing - new...\n");
}

