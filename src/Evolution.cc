#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

#include "Evolution.hh"

extern int testProblem;
extern int D;
extern int Nc;
extern int Nv;
extern int N[3];
extern int NV[3];

extern double Cv;
extern double R;
extern double gma;
extern double K;


double TimeStep(){

	double dt = 128;
	return dt;
}


void Evolve(double* g, double* b, double* gbar, double* bbar, double* gbarp, double* bbarp, double* Sg, double* Sb, double* rho, double* rhov, double* rhoE, int effD, double dt, double Tf, double Tsim, double dtdump, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ){	

	dt = fmin(TimeStep(), dtdump);
	dt = fmin(dt, Tf-Tsim);

	Step1a(g, b, gbar, bbar, gbarp, bbarp, Sg, Sb, rho, rhov, rhoE, effD, dt, Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ);
	Step1b();
	Step1c();
	
	Step2a();
	Step2b();
	Step2c();
	
	Step3();
	
	Step4();
	
	Step5();
}


//Step 2
void Step1a(double* g, double* b, double* gbar, double* bbar, double* gbarp, double* bbarp, double* Sg, double* Sb, double* rho, double* rhov, double* rhoE, int effD, double dt, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ){

	double tau = 1;

	int Nx = N[0];
	int Ny = N[1];
	int Nz = N[2];

	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){
				for(int vx = 0; vx < NV[0]; vx++){
					for(int vy = 0; vy < NV[1]; vy++){
						for(int vz = 0; vz < NV[2]; vz++){

							int sidx = i + Nx*j + Nx*Ny*k; //spatial index
							int idx = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;
							//printf("idx = %d", idx);

							double U[3] = {rhov[D*sidx]/rho[sidx], rhov[D*sidx + 1]/rho[sidx], rhov[D*sidx + 2]/rho[sidx]};
							double T = Temperature(rhoE[sidx]/rho[sidx], sqrt(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]));


							double g_eq = geq(Co_X[vx], Co_Y[vy], Co_Z[vz], rho[sidx], U, T, Co_WX[vx], Co_WY[vy], Co_WZ[vz]);
							double b_eq = g_eq*(Co_X[vx]*Co_X[vx] + Co_Y[vy]*Co_Y[vy] + Co_Z[vz]*Co_Z[vz] + (3-effD+K)*R*T)/2;
							gbarp[idx] = (2*tau - dt/2.)/(2*tau)*g[idx] + dt/(4*tau)*g_eq + dt/4*Sg[sidx];
							bbarp[idx] = (2*tau - dt/2.)/(2*tau)*b[idx] + dt/(4*tau)*b_eq + dt/4*Sb[sidx];



						}
					}
				}
			}
		}
	}

}
void Step1b(){

}
void Step1c(){

}

//Step 2
void Step2a(){

}
void Step2b(){

}
void Step2c(){

}

//Step 3
void Step3(){

}

//Step 4
void Step4(){

}

//Step 5
void Step5(){

}