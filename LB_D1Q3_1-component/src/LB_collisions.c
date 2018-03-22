/*
 * LB_collisions.c
 *
 *  Created on: Aug 5, 2016
 *      Author: clark
 */

#include "LB_D1Q3_1-component.h"


void collisionForcingBasic() {
	double a,A,tmp,phi,dn,ddn,dp,ddp;
	double omega = 1./oneOverTau;
	int i,ip,im;

	iterations++;
	for (i=0;i<XDIM;i++){
		n1[i]=f1_0[i]+f1_1[i]+f1_2[i];
	}

	calculatePressures();

	for (i=0;i<XDIM;i++){
		a=f1_1[i]-f1_2[i];
		A=a*a/n1[i]+n1[i]/3.;
		f1_0[i]+=omega*(n1[i]-A-f1_0[i])-Forces[i]*(2*a)/n1[i];
		f1_1[i]+=omega*(0.5*( a+A)-f1_1[i])+0.5*( Forces[i]+Forces[i]*(2*a)/n1[i]);
		f1_2[i]+=omega*(0.5*(-a+A)-f1_2[i])+0.5*(-Forces[i]+Forces[i]*(2*a)/n1[i]);
	}

} // end function collision()


void collisionForcingAlexander(){
	double a,A,tmp,phi,dn,ddn,dp,ddp;
	double omega = 1./oneOverTau;
	int i,ip,im;

	iterations++;
	for (i=0;i<XDIM;i++){
		n1[i]=f1_0[i]+f1_1[i]+f1_2[i];
	}

	calculatePressures();

	for (i=0;i<XDIM;i++){
		a=f1_1[i]-f1_2[i];
		A=a*a/n1[i]+n1[i]/3.;
		f1_0[i]+=omega*(n1[i]-A-f1_0[i])-Forces[i]*(2*a+Forces[i])/n1[i];
		f1_1[i]+=omega*(0.5*( a+A)-f1_1[i])+0.5*( Forces[i]+Forces[i]*(2*a+Forces[i])/n1[i]);
		f1_2[i]+=omega*(0.5*(-a+A)-f1_2[i])+0.5*(-Forces[i]+Forces[i]*(2*a+Forces[i])/n1[i]);
	}
}


void collisionForcingKyoto(){
	double a,A,tmp,phi,dn,dp,ddp;
	double omega = 1./oneOverTau;
	int i,ip,im;

	iterations++;
	for (i=0;i<XDIM;i++){
		n1[i]=f1_0[i]+f1_1[i]+f1_2[i];
	}
	for (i=0;i<XDIM;i++){
		ip=(i+1)%XDIM;
		im=(i+XDIM-1)%XDIM;
		dn=0.5*(n1[ip]-n1[im]);
		phi=(n1[i]-nc)/nc;

		pni[i]=Pni(n1[i],dn,ddni[i]);
		p[i] = P(n1[i],dn,ddni[i]);
	}

	calculatePressures();

	for (i=0;i<XDIM;i++){
		a=f1_1[i]-f1_2[i];
		A=a*a/n1[i]-0.25*Forces[i]*Forces[i]/n1[i]+1./12.*ddni[i]+n1[i]/3.;
		f1_0[i]+=omega*(n1[i]-A-f1_0[i])-Forces[i]*(2*a+Forces[i])/n1[i];
		f1_1[i]+=omega*(0.5*( a+A)-f1_1[i])+0.5*( Forces[i]+Forces[i]*(2*a+Forces[i])/n1[i]);
		f1_2[i]+=omega*(0.5*(-a+A)-f1_2[i])+0.5*(-Forces[i]+Forces[i]*(2*a+Forces[i])/n1[i]);
	}
}


void collisionIdeal(){
	double n,a,A,tmp;
	double omega = 1./oneOverTau;
	int i;

	iterations++;
	for (i=0;i<XDIM;i++){
		n1[i]=f1_0[i]+f1_1[i]+f1_2[i];
		a=0.5*(f1_1[i]-f1_2[i]);
		A=0.5*(4*a*a/n1[i]+n1[i]/3.);
		f1_0[i]+=omega*(n1[i]-2*A-f1_0[i]);
		f1_1[i]+=omega*(  a+A-f1_1[i]);
		f1_2[i]+=omega*( -a+A-f1_2[i]);
	}
}


void setCollisionForcingNewPressureGradient() {
	collision = collisionForcingNewPressureGradient;
	printf("gradP forcing - new...\n");
}


void setCollisionForcingNewChemicalPotentialGradient() {
	collision = collisionForcingNewChemicalPotentialGradient;
	printf("gradMu forcing - new...\n");
}


void setCollisionForcingBasic() {
	collision = collisionForcingBasic;
}


void setCollisionForcingAlexander() {
	collision = collisionForcingAlexander;
}


void setCollisionForcingKyoto() {
	collision = collisionForcingKyoto;
}


void setCollisionPressureMethod() {
	collision = collisionPressureMethod;
	printf("pressure method...\n");
}


void setCollisionIdeal() {
	collision = collisionIdeal;
}

