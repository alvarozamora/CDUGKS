#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

#include "Mesh.hh"
#include "Evolution.hh"

extern int testProblem;
extern int D;
extern int Nc;
extern int Nv;
extern int N[3];
extern int NV[3];
extern int effD;

extern double Cv;
extern double R;
extern double gma;
extern double K;


double TimeStep(){

	double dt = 128;
	return dt;
}


void Evolve(double* g, double* b, double* gbar, double* bbar, double* gbarp, double* bbarp, double* Sg, double* Sb, double* rho, double* rhov, double* rhoE, int effD, double dt, double Tf, double Tsim, double dtdump, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double* gsigma, double* bsigma, double* gsigma2, double* bsigma2, Cell* mesh){	

	dt = fmin(TimeStep(), dtdump);
	dt = fmin(dt, Tf-Tsim);

	Step1a(g, b, gbar, bbar, gbarp, bbarp, Sg, Sb, rho, rhov, rhoE, effD, dt, Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ);
	Step1b(gbarp, bbarp, effD, gsigma, bsigma, gsigma2, bsigma2, mesh);
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


							double mu = visc(T);

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
void Step1b(double* gbarp, double* bbarp, int effD, double* gsigma, double* bsigma, double* gsigma2, double* bsigma2, Cell* mesh){

	int Nx = N[0];
	int Ny = N[1];
	int Nz = N[2];

	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){
				for(int vx = 0; vx < NV[0]; vx++){
					for(int vy = 0; vy < NV[1]; vy++){
						for(int vz = 0; vz < NV[2]; vz++){


							//Compute Sigma
							int idx = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;

							for(int Dim = 0; Dim < effD; Dim++){
								int IL, IR, JL, JR, KL, KR;

								//Periodic Boundary Conditions
								if(Dim == 0){IL = (i - 1)%N[0]; IR = (i + 1)%N[0]; JL = j; JR = j; KL = k; KR = k;} 
								if(Dim == 1){JL = (j - 1)%N[1]; JR = (j + 1)%N[1]; IL = i; IR = i; KL = k; KR = k;}
								if(Dim == 2){KL = (k - 1)%N[2]; KR = (k + 1)%N[2]; IL = i; IR = i; JL = j; JR = k;}

								/*
								//Reflective Boundary Conditions
								if(Dim == 0){IL = (i - 1); IR = (i + 1);} 
								if(Dim == 1){JL = (j - 1); JR = (j + 1);}
								if(Dim == 2){KL = (k - 1); KR = (k + 1);}

								if(IL < 0){IL = 0;} if(IR == N[0]){IR = N[0] - 1;}
								if(JL < 0){JL = 0;} if(JR == N[1]){JR = N[1] - 1;}
								if(KL < 0){KL = 0;} if(KR == N[2]){KR = N[2] - 1;}
								*/

								int idxL = IL + Nx*JL + Nx*Ny*KL + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;
								int idxR = IR + Nx*JR + Nx*Ny*KR + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;


								double xL[3] = {mesh[idxL].x, mesh[idxL].y, mesh[idxL].z};
								double xR[3] = {mesh[idxR].x, mesh[idxR].y, mesh[idxR].z};
								double xC[3] = {mesh[idx].x,  mesh[idx].y,  mesh[idx].z};

								gsigma[effD*idx + Dim] = VanLeer(gbarp[idxL], gbarp[idx], gbarp[idxR], xL[Dim], xC[Dim], xR[Dim]);
								bsigma[effD*idx + Dim] = VanLeer(bbarp[idxL], bbarp[idx], bbarp[idxR], xL[Dim], xC[Dim], xR[Dim]);


							}



							//Compute gbarp/bbarp @ interfaces via interpolation 


						}
					}
				}
			}
		}
	}


}
void Step1c(){

	//Compute gbar/bbar @ t=n+1/2  with interface sigma
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