#include <stdio.h>
#include <math.h>

double a=-0.375,b=1./3.,T=0.333333;

double F(double N1,double V1){
	double V2=1-V1;
	double N2=1-N1;

	double F1;
	if ((N1*b>=1)||(N1<=0)||(V1<0)) F1=1e10;
	else F1= N1*T*log(N1/(V1-N1*b))+a*N1*N1/V1;
	double F2;
	if ((N2*b>=1)||(N2<=0)||(V2<0)) F2=1e10;
	else F2=N2*T*log(N2/(V2-N2*b))+a*N2*N2/V2;
	return F1+F2;
}

void Minimize(double *N1,double *V1){
	double FF[3][3];
	double D=0.1;

	while (D>0.000001){
		for (int i=0;i<3;i++)
			for (int j=0;j<3;j++){
				FF[i][j]=F(*N1+D*(i-1),*V1+D*(j-1));
			}
		int mini=1,minj=1;
		double minF=FF[1][1];
		for (int i=0;i<3;i++)
			for (int j=0;j<3;j++){
				if (FF[i][j]<minF){
					mini=i;
					minj=j;
					minF=FF[i][j];
				}
			}
		/*    printf("D=%g",D);
	  printf(" %i,%i \n",mini,minj);
    for (int i=0; i<3; i++){
      for (int j=0;j<3; j++)
	printf("%g",FF[i][j]);
      printf("\n");
    }
		 */
		if ((mini==1)&&(minj==1)) D/=2;
		else {
			*N1+=(mini-1)*D;
			*V1+=(minj-1)*D;
		}
	}
}

void minimize(){
	double N1=0.5,V1=0.5;
	FILE *phaseDiagram_rhoVsTemp_theoretical;
	phaseDiagram_rhoVsTemp_theoretical = fopen("/home/clark/school/Lattice Boltzmann/Maxwell construction/phaseDiagram_rhoVsTemp_theoretical.dat","w");

	for (T=1./3.;T>0;T-=0.00001){
		Minimize(&N1,&V1);
		printf("%e %e %e\n",T,N1/V1,(1-N1)/(1-V1));
		fprintf(phaseDiagram_rhoVsTemp_theoretical,"%.15f %.15f\n%.15f %.15f\n", N1/V1, T/(1./3.), (1-N1)/(1-V1), T/(1./3.));
	}

	fclose(phaseDiagram_rhoVsTemp_theoretical);
}
