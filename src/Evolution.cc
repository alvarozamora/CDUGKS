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
extern double Pr;
extern double Vmax[3];


int Nx = N[0];
int Ny = N[1];
int Nz = N[2];

int debug = 1.0;

double TimeStep(double calcdt, double dumptime, double tend){

	
	double timestep = fmin(fmin(calcdt, dumptime), tend);
	return timestep;
}


int Evolve(double* g, double* b, double* gbar, double* bbar, double* gbarp, double* bbarp, double* Sg, double* Sb, double* rho, double* rhov, double* rhoE, double* dt, double Tf, double Tsim, double dtdump, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double* gsigma, double* bsigma, double* gsigma2, double* bsigma2, Cell* mesh, double* gbarpbound, double* bbarpbound, double* rhoh, double* rhovh, double* rhoEh, double* Tdump, int* BCs){	

	//Find timestep
	double CFL = 0.9; //safety factor
	double dxmin = 1.0/fmax(fmax(N[0],N[1]),N[2]); //smallest cell width 
	double umax = 2.0; // estimated maximum flow velocity , TODO calculate at each iteration for stronger problems
	double calcdt = CFL*dxmin/(1.0+sqrt(Vmax[0]*Vmax[0] + Vmax[1]*Vmax[1] + Vmax[2]*Vmax[2]));

	*dt = TimeStep(calcdt, dtdump-*Tdump, Tf-Tsim);
	
	int dump = (*dt < calcdt);


	//Evolution Cycle
	Step1a(g, b, gbar, bbar, gbarp, bbarp, Sg, Sb, rho, rhov, rhoE, *dt, Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ);
	Step1b(gbarp, bbarp, gsigma, bsigma, gsigma2, bsigma2, mesh, gbarpbound, bbarpbound, Co_X, Co_Y, Co_Z, BCs);
	Step1c(gbar, bbar, gbarpbound, bbarpbound, Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ, gsigma2, bsigma2, *dt, BCs);
	
	Step2a(gbar, bbar, Co_X, Co_Y, Co_Z, Co_WX, Co_WY, Co_WZ, *dt, rhoh, rhovh, rhoEh);
	Step2b(gbar, bbar, *dt, rhoh, rhovh, rhoEh, Co_X, Co_Y, Co_Z); //gbar, bbar are actually g/b at interface, recycling memory
	Step2c(gbar, bbar, Co_X, Co_Y, Co_Z, mesh, gbarp, bbarp, BCs); //gbar, bbar are actually g/b at interface, gbarp/bbarp are actually Fg/Fb -- Recycling memory
	
	Step3();
	
	Step4and5(rho, rhov, rhoE, *dt, mesh, gbarp, bbarp, Co_X, Co_Y, Co_Z, Co_WX, Co_WY, Co_WZ, g, b, BCs);
	

	return dump;
}


//Step 1: Phibar at interface
//Step 1a: Phibar at Cell Center.
void Step1a(double* g, double* b, double* gbar, double* bbar, double* gbarp, double* bbarp, double* Sg, double* Sb, double* rho, double* rhov, double* rhoE, double dt, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ){

	if(debug == 1){printf("Entering Step 1a\n");}

	double tg;
	double tb;
	

	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){
				int sidx = i + Nx*j + Nx*Ny*k; //spatial index

				

				for(int vx = 0; vx < NV[0]; vx++){
					for(int vy = 0; vy < NV[1]; vy++){
						for(int vz = 0; vz < NV[2]; vz++){

							
							int idx = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;

							double u = 0;
							for(int dim = 0; dim < effD; dim++){u += rhov[effD*sidx + dim]/rho[sidx]*rhov[effD*sidx + dim]/rho[sidx];} u = sqrt(u);
							double T = Temperature(rhoE[sidx]/rho[sidx], u);


							tg = visc(T)/rho[sidx]/R/T; // tau = mu/P, P = rho*R*T.
							tb = tg/Pr;

							//For Now...
							Sg[sidx] = 0.;
							Sb[sidx] = 0.;

							double c2 = 0;
							double Xi[3] = {Co_X[vx], Co_Y[vy], Co_Z[vz]};

							for(int dim = 0 ; dim < effD; dim++){ c2 += (Xi[dim]-rhov[effD*sidx + dim]/rho[sidx])*(Xi[dim]-rhov[effD*sidx + dim]/rho[sidx]);} //TODO: Potential BUG, did not double check algebra.

							double g_eq = geq(c2, rho[sidx], T);
							double b_eq = g_eq*(Co_X[vx]*Co_X[vx] + Co_Y[vy]*Co_Y[vy] + Co_Z[vz]*Co_Z[vz] + (3-effD+K)*R*T)/2;

							

							gbarp[idx] = (2*tg - dt/2.)/(2*tg)*g[idx] + dt/(4*tg)*g_eq + dt/4*Sg[sidx];
							bbarp[idx] = (2*tb - dt/2.)/(2*tb)*b[idx] + dt/(4*tb)*b_eq + dt/4*Sb[sidx];

							assert(gbarp[idx] == gbarp[idx]); // NaN Checker
							assert(bbarp[idx] == bbarp[idx]); // NaN Checker
						}
					}
				}
			}
		}
	}
}

//Step 1b: compute gradient of phibar to compute phibar at interface. compute phibar at interface.
void Step1b(double* gbarp, double* bbarp, double* gsigma, double* bsigma, double* gsigma2, double* bsigma2, Cell* mesh, double* gbarpbound, double* bbarpbound, double* Co_X, double* Co_Y, double* Co_Z, int* BCs){
	
	if(debug == 1){printf("Entering Step 1b\n");}

	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){

				int sidx = i + Nx*j + Nx*Ny*k;

				double xC[3] = {mesh[sidx].x,  mesh[sidx].y,  mesh[sidx].z};
				double sC[3] = {mesh[sidx].dx,  mesh[sidx].dy,  mesh[sidx].dz};

				for(int vx = 0; vx < NV[0]; vx++){
					for(int vy = 0; vy < NV[1]; vy++){
						for(int vz = 0; vz < NV[2]; vz++){


							//Compute Sigma
							int idx = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;

							

							for(int Dim = 0; Dim < effD; Dim++){
								int IL, IR, JL, JR, KL, KR;

								//Periodic Boundary Conditions
								if(Dim == 0 && BCs[0] == 0){IL = (i - 1 + N[0])%N[0]; IR = (i + 1)%N[0]; JL = j; JR = j; KL = k; KR = k;} 
								if(Dim == 1 && BCs[1] == 0){JL = (j - 1 + N[1])%N[1]; JR = (j + 1)%N[1]; IL = i; IR = i; KL = k; KR = k;}
								if(Dim == 2 && BCs[2] == 0){KL = (k - 1 + N[2])%N[2]; KR = (k + 1)%N[2]; IL = i; IR = i; JL = j; JR = j;}


								//Dirichlet Boundary Conditions
								if(Dim == 0 && BCs[0] == 1){IL = (i - 1); IR = (i + 1); if(IL < 0){IL = 0;} if(IR == N[0]){IR = N[0] - 1;} JL = j; JR = j; KL = k; KR = k;} 
								if(Dim == 1 && BCs[1] == 1){JL = (j - 1); JR = (j + 1); if(JL < 0){JL = 0;} if(JR == N[1]){JR = N[1] - 1;} IL = i; IR = i; KL = k; KR = k;}
								if(Dim == 2 && BCs[2] == 1){KL = (k - 1); KR = (k + 1); if(KL < 0){KL = 0;} if(KR == N[2]){KR = N[2] - 1;} IL = i; IR = i; JL = j; JR = j;}
								
								
								//Neumann Boundary Conditions
								if(Dim == 0 && BCs[0] == 2){IL = (i - 1); IR = (i + 1); if(IL < 0){IL = 0;} if(IR == N[0]){IR = N[0] - 1;} JL = j; JR = j; KL = k; KR = k;} 
								if(Dim == 1 && BCs[1] == 2){JL = (j - 1); JR = (j + 1); if(JL < 0){JL = 0;} if(JR == N[1]){JR = N[1] - 1;} IL = i; IR = i; KL = k; KR = k;}
								if(Dim == 2 && BCs[2] == 2){KL = (k - 1); KR = (k + 1); if(KL < 0){KL = 0;} if(KR == N[2]){KR = N[2] - 1;} IL = i; IR = i; JL = j; JR = j;}
								
								//Gather Left and Right Indices in real/vel space
								int idxL = IL + Nx*JL + Nx*Ny*KL + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;
								int idxR = IR + Nx*JR + Nx*Ny*KR + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;

								//Gather Left and Right Spacial Indices
								int sidxL = IL + Nx*JL + Nx*Ny*KL;
								int sidxR = IR + Nx*JR + Nx*Ny*KR;
								

								//Gather position of left/right cell centers
								double xL[3] = {mesh[sidxL].x, mesh[sidxL].y, mesh[sidxL].z};
								double xR[3] = {mesh[sidxR].x, mesh[sidxR].y, mesh[sidxR].z};
								
								//Gather cell size of left/right cells
								double sL[3] = {mesh[sidxL].dx, mesh[sidxL].dy, mesh[sidxL].dz};
								double sR[3] = {mesh[sidxR].dx, mesh[sidxR].dy, mesh[sidxR].dz};
								

								//Computing phisigma, at cell 
								gsigma[effD*idx + Dim] = VanLeer(gbarp[idxL], gbarp[idx], gbarp[idxR], xL[Dim], xC[Dim], xR[Dim]);
								bsigma[effD*idx + Dim] = VanLeer(bbarp[idxL], bbarp[idx], bbarp[idxR], xL[Dim], xC[Dim], xR[Dim]);

								//printf("Checking bsigma input[%d][%d]: bbarp[idxL] = %f, bbarp[idx] = %f, bbarp[idxR] = %f, xL[Dim] = %f, xC[Dim] = %f, xR[Dim] = %f\n", sidx, vx, bbarp[idxL], bbarp[idx], bbarp[idxR], xL[Dim], xC[Dim], xR[Dim]);

								

								//Computing phisigma, at interface
								for(int Dim2 = 0; Dim2 < effD; Dim2++){

									int IL2, IR2, JL2, JR2, KL2, KR2;

									//Periodic Boundary Conditions
									if(Dim2 == 0 && BCs[0] == 0){IL2 = (i - 1 + N[0])%N[0]; IR2 = (i + 1)%N[0]; JL2 = j; JR2 = j; KL2 = k; KR2 = k;} 
									if(Dim2 == 1 && BCs[1] == 0){JL2 = (j - 1 + N[1])%N[1]; JR2 = (j + 1)%N[1]; IL2 = i; IR2 = i; KL2 = k; KR2 = k;}
									if(Dim2 == 2 && BCs[2] == 0){KL2 = (k - 1 + N[2])%N[2]; KR2 = (k + 1)%N[2]; IL2 = i; IR2 = i; JL2 = j; JR2 = j;}
	
									
									//Dirichlet Boundary Conditions
									if(Dim2 == 0 && BCs[0] == 1){IL2 = (i - 1); IR2 = (i + 1); if(IL2 < 0){IL2 = 0;} if(IR2 == N[0]){IR2 = N[0] - 1;} JL2 = j; JR2 = j; KL2 = k; KR2 = k;} 
									if(Dim2 == 1 && BCs[1] == 1){JL2 = (j - 1); JR2 = (j + 1); if(JL2 < 0){JL2 = 0;} if(JR2 == N[1]){JR2 = N[1] - 1;} IL2 = i; IR2 = i; KL2 = k; KR2 = k;}
									if(Dim2 == 2 && BCs[2] == 1){KL2 = (k - 1); KR2 = (k + 1); if(KL2 < 0){KL2 = 0;} if(KR2 == N[2]){KR2 = N[2] - 1;} IL2 = i; IR2 = i; JL2 = j; JR2 = j;}


									//Neumann Boundary Conditions
									if(Dim2 == 0 && BCs[0] == 2){IL2 = (i - 1); IR2 = (i + 1); if(IL2 < 0){IL2 = 0;} if(IR2 == N[0]){IR2 = N[0] - 1;} JL2 = j; JR2 = j; KL2 = k; KR2 = k;} 
									if(Dim2 == 1 && BCs[1] == 2){JL2 = (j - 1); JR2 = (j + 1); if(JL2 < 0){JL2 = 0;} if(JR2 == N[1]){JR2 = N[1] - 1;} IL2 = i; IR2 = i; KL2 = k; KR2 = k;}
									if(Dim2 == 2 && BCs[2] == 2){KL2 = (k - 1); KR2 = (k + 1); if(KL2 < 0){KL2 = 0;} if(KR2 == N[2]){KR2 = N[2] - 1;} IL2 = i; IR2 = i; JL2 = j; JR2 = j;}


									int idxL2 = IL2 + Nx*JL2 + Nx*Ny*KL2 + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;
									int idxR2 = IR2 + Nx*JR2 + Nx*Ny*KR2 + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;

									int sidxL2 = IL2 + Nx*JL2 + Nx*Ny*KL2;
									int sidxR2 = IR2 + Nx*JR2 + Nx*Ny*KR2;
									

									double xL2[3] = {mesh[sidxL2].x, mesh[sidxL2].y, mesh[sidxL2].z};
									double xR2[3] = {mesh[sidxR2].x, mesh[sidxR2].y, mesh[sidxR2].z};
									double xC[3] = {mesh[sidx].x,  mesh[sidx].y,  mesh[sidx].z};

									//printf("sidx = %d, sidxL2 = %d, sidxR2 = %d\n", sidx, sidxL2, sidxR2);
									//printf("xL2 = {%f, %f, %f}\n", xL2[0], xL2[1], xL2[2]);
									//printf("xC = {%f, %f, %f}\n", xC[0], xC[1], xC[2]);

									double sL2[3] = {mesh[sidxL2].dx, mesh[sidxL2].dy, mesh[sidxL2].dz};
									double sR2[3] = {mesh[sidxR2].dx, mesh[sidxR2].dy, mesh[sidxR2].dz};
									double sC2[3] = {mesh[sidx].dx,  mesh[sidx].dy,  mesh[sidx].dz};


									//Dim 1 is vector component that is being interpolated.
									//Dim 2 is direction of interpolation.

									gsigma2[effD*effD*idx + effD*Dim + Dim2] = gsigma[effD*idx + Dim] + (sC2[Dim2]/2)*VanLeer(gsigma[effD*idxL2 + Dim], gsigma[effD*idx + Dim], gsigma[effD*idxR2 + Dim], xL2[Dim2], xC[Dim2], xR2[Dim2]);
									bsigma2[effD*effD*idx + effD*Dim + Dim2] = bsigma[effD*idx + Dim] + (sC2[Dim2]/2)*VanLeer(bsigma[effD*idxL2 + Dim], bsigma[effD*idx + Dim], bsigma[effD*idxR2 + Dim], xL2[Dim2], xC[Dim2], xR2[Dim2]);
									
									//printf("xL2[%d], xC[%d], xR2[%d] = %f, %f, %f\n", Dim2, Dim2, Dim2, xL2[Dim2], xC[Dim2], xR2[Dim2]);
									//printf("sidx = %d, gbarp[%d] = %f\n", sidx, idxR, gbarp[idxR]);
									//printf("gsigma = %f\n" ,gsigma[effD*idx + Dim]);
									//printf("gsigma2 = %f\n",gsigma2[effD*effD*idx + effD*Dim + Dim2]);
									
								}


								//Dot Product is just a single product when using rectangular mesh
								int interpidx = idx;
								double swap = 1.;
								     if(Co_X[vx] < 0 && Dim == 0){interpidx = idxR; swap = -1;}
								else if(Co_Y[vy] < 0 && Dim == 1){interpidx = idxR; swap = -1;}
								else if(Co_Z[vz] < 0 && Dim == 2){interpidx = idxR; swap = -1;}

								
								//TODO need to change sC to sR when swap, doesnt currently matter for sod bc dx_i = dx_0
								gbarpbound[effD*idx + Dim] = gbarp[interpidx] + swap*sC[Dim]/2*gsigma[effD*interpidx + Dim];
								bbarpbound[effD*idx + Dim] = bbarp[interpidx] + swap*sC[Dim]/2*bsigma[effD*interpidx + Dim];

								//printf("Checking bbarpbound[%d][%d]: bbarp = %f, bsigma = %f\n", sidx,vx, bbarp[interpidx], bsigma[effD*interpidx + Dim]);
								assert(gbarpbound[effD*idx + Dim] == gbarpbound[effD*idx + Dim]); // NaN Checker
								assert(bbarpbound[effD*idx + Dim] == bbarpbound[effD*idx + Dim]); // NaN Checker

								
							}
						}
					}
				}
			}
		}
	}

}


// Step 1c: Compute phibar at interface by interpolating w/ phisigma2, x-Xi*dt/2
void Step1c(double* gbar, double* bbar, double* gbarpbound, double*bbarpbound, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double* gsigma2, double* bsigma2, double dt, int* BCs){

	if(debug == 1){printf("Entering Step 1c\n");}

	//Compute gbar/bbar @ t=n+1/2  with interface sigma

	int Nx = N[0];
	int Ny = N[1];
	int Nz = N[2];

	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){
				for(int vx = 0; vx < NV[0]; vx++){
					for(int vy = 0; vy < NV[1]; vy++){
						for(int vz = 0; vz < NV[2]; vz++){
							int idx = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;
							for(int Dim = 0; Dim < effD; Dim++){

								double Xi[3] = {Co_X[vx], Co_Y[vy], Co_Z[vz]};

								//Phibar at Interface, at t = n+1/2
								gbar[effD*idx + Dim] = gbarpbound[effD*idx + Dim];
								bbar[effD*idx + Dim] = bbarpbound[effD*idx + Dim];


								for(int Dim2 = 0; Dim2 < effD; Dim2++){

									int IL2, IR2, JL2, JR2, KL2, KR2;
									//Periodic Boundary Conditions
									if(Dim2 == 0 && BCs[0] == 0){IL2 = (i - 1 + N[0])%N[0]; IR2 = (i + 1)%N[0]; JL2 = j; JR2 = j; KL2 = k; KR2 = k;} 
									if(Dim2 == 1 && BCs[1] == 0){JL2 = (j - 1 + N[1])%N[1]; JR2 = (j + 1)%N[1]; IL2 = i; IR2 = i; KL2 = k; KR2 = k;}
									if(Dim2 == 2 && BCs[2] == 0){KL2 = (k - 1 + N[2])%N[2]; KR2 = (k + 1)%N[2]; IL2 = i; IR2 = i; JL2 = j; JR2 = j;}
	
									
									//Dirichlet Boundary Conditions
									if(Dim2 == 0 && BCs[0] == 1){IL2 = (i - 1); IR2 = (i + 1); if(IL2 < 0){IL2 = 0;} if(IR2 == N[0]){IR2 = N[0] - 1;} JL2 = j; JR2 = j; KL2 = k; KR2 = k;} 
									if(Dim2 == 1 && BCs[1] == 1){JL2 = (j - 1); JR2 = (j + 1); if(JL2 < 0){JL2 = 0;} if(JR2 == N[1]){JR2 = N[1] - 1;} IL2 = i; IR2 = i; KL2 = k; KR2 = k;}
									if(Dim2 == 2 && BCs[2] == 1){KL2 = (k - 1); KR2 = (k + 1); if(KL2 < 0){KL2 = 0;} if(KR2 == N[2]){KR2 = N[2] - 1;} IL2 = i; IR2 = i; JL2 = j; JR2 = j;}


									//Neumann Boundary Conditions
									if(Dim2 == 0 && BCs[0] == 2){IL2 = (i - 1); IR2 = (i + 1); if(IL2 < 0){IL2 = 0;} if(IR2 == N[0]){IR2 = N[0] - 1;} JL2 = j; JR2 = j; KL2 = k; KR2 = k;} 
									if(Dim2 == 1 && BCs[1] == 2){JL2 = (j - 1); JR2 = (j + 1); if(JL2 < 0){JL2 = 0;} if(JR2 == N[1]){JR2 = N[1] - 1;} IL2 = i; IR2 = i; KL2 = k; KR2 = k;}
									if(Dim2 == 2 && BCs[2] == 2){KL2 = (k - 1); KR2 = (k + 1); if(KL2 < 0){KL2 = 0;} if(KR2 == N[2]){KR2 = N[2] - 1;} IL2 = i; IR2 = i; JL2 = j; JR2 = j;}


									int idxL2 = IL2 + Nx*JL2 + Nx*Ny*KL2 + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;
									int idxR2 = IR2 + Nx*JR2 + Nx*Ny*KR2 + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;

									int interpidx = idx;
									double swap = 1.;
								         if(Co_X[vx] < 0 && Dim2 == 0){interpidx = idxR2; swap = -1;}
									else if(Co_Y[vy] < 0 && Dim2 == 1){interpidx = idxR2; swap = -1;}
									else if(Co_Z[vz] < 0 && Dim2 == 2){interpidx = idxR2; swap = -1;}

									gbar[effD*idx + Dim] -= dt/2.0*Xi[Dim2]*gsigma2[effD*effD*interpidx + effD*Dim2 + Dim];
									bbar[effD*idx + Dim] -= dt/2.0*Xi[Dim2]*bsigma2[effD*effD*interpidx + effD*Dim2 + Dim];

									assert(gbar[effD*idx + Dim] == gbar[effD*idx + Dim]); // NaN Checker
									assert(bbar[effD*idx + Dim] == bbar[effD*idx + Dim]); // NaN Checker
								}	
							}
						}
					}
				}
			}
		}
	}
}

//Step 2: Microflux
//Step 2a: Interpolate W to interface.
void Step2a(double* gbar, double* bbar, double* Co_X, double* Co_Y, double* Co_Z, double* Co_WX, double* Co_WY, double* Co_WZ, double dt, double* rhoh, double* rhovh, double* rhoEh){
	
	if(debug == 1){printf("Entering Step 2a\n");}
	

	//Compute conserved variables W at t+1/2
	//First do density at boundary, density is needed for others.
	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){
				int sidx = i + Nx*j + Nx*Ny*k;
				for(int d = 0; d < effD; d++){

					rhoh[effD*sidx + d] = 0; 

					for(int vx = 0; vx < NV[0]; vx++){
						for(int vy = 0; vy < NV[1]; vy++){
							for(int vz = 0; vz < NV[2]; vz++){


								int idx = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;

								rhoh[effD*sidx + d] += Co_WX[vx]*Co_WY[vy]*Co_WZ[vz]*gbar[effD*idx + d]; 
							}
						}
					}
				}
			}
		}
	}

	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){
				int sidx = i + Nx*j + Nx*Ny*k;
				//Inialize E and momentum with source term
				
				//Dim is vector component that was interpolated
				//Dim2 is direction of interpolation (toward interface)
				for(int Dim2 = 0; Dim2 < effD; Dim2++){

					rhoEh[effD*sidx + Dim2] = 0.0;  // Initialize 0 before effD dot product.

					for(int Dim = 0; Dim < effD; Dim++){

						rhovh[effD*effD*sidx + effD*Dim + Dim2] = dt/2*rhoh[effD*sidx+Dim2]*0; //TODO: In future, replace 0 with acceleration field
						rhoEh[effD*sidx + Dim2] += dt/2.*rhoh[effD*sidx + Dim]*0; //TODO: In future replace 0 with u.dot(a), vel dot acc
					}

				}

				for(int vx = 0; vx < NV[0]; vx++){
					for(int vy = 0; vy < NV[1]; vy++){
						for(int vz = 0; vz < NV[2]; vz++){

							int idx = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;

							double U[3] = {Co_X[vx], Co_Y[vy], Co_Z[vz]};

							//Dim is vector component that was interpolated
							//Dim2 is direction of interpolation (toward interface)
							for(int Dim2 = 0; Dim2 < effD; Dim2++){
								for(int Dim = 0; Dim < effD; Dim++){

									rhovh[effD*effD*sidx + effD*Dim + Dim2] += Co_WX[vx]*Co_WY[vy]*Co_WZ[vz]*U[Dim2]*gbar[effD*idx + Dim]; 
									
								}

								rhoEh[effD*sidx + Dim2] += Co_WX[vx]*Co_WY[vy]*Co_WZ[vz]*bbar[effD*idx + Dim2]; 
							}	
						}
					}
				}
			}
		}
	}
}

//Step 2b: compute original phi at interface using gbar, W at interface
//Memory Recycling: phibar @ interface is used to store phi @ interface.
void Step2b(double* gbar, double* bbar, double dt, double* rhoh, double* rhovh, double* rhoEh, double* Co_X, double* Co_Y, double* Co_Z){

	if(debug == 1){printf("Entering Step 2b\n");}

	double tg;
	double tb;
	double u;
	double T;
	int idx;
	int sidx;
	double g_eq;
	double b_eq;

	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){
				sidx = i + Nx*j + Nx*Ny*k;
				for(int dim2 = 0; dim2 < effD; dim2++){ 

					u = 0;
					//Dim is vector component that was interpolated
					//Dim2 is direction of interpolation (toward interface)
					for(int dim = 0; dim < effD; dim++){u += rhovh[effD*effD*sidx + dim*effD + dim2]/rhoh[effD*sidx + dim2]*rhovh[effD*effD*sidx + dim*effD + dim2]/rhoh[effD*sidx + dim2];} u = sqrt(u);
					T = Temperature(rhoEh[effD*sidx+ dim2]/rhoh[effD*sidx + dim2], u);

					tg = visc(T)/rhoh[effD*sidx + dim2]/R/T;
					tb = tg/Pr;

					for(int vx = 0; vx < NV[0]; vx++){
						for(int vy = 0; vy < NV[1]; vy++){
							for(int vz = 0; vz < NV[2]; vz++){
								idx = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;


							
								double c2 = 0;
								double Xi[3] = {Co_X[vx], Co_Y[vy], Co_Z[vz]};

								for(int dim = 0 ; dim < effD; dim++){ c2 += (Xi[dim]-rhovh[effD*effD*sidx + effD*dim + dim2]/rhoh[effD*sidx + dim2])*(Xi[dim]-rhovh[effD*effD*sidx + effD*dim + dim2]/rhoh[effD*sidx + dim2]);} //TODO: Potential BUG, did not double check algebra or access dim order

								g_eq = geq(c2, rhoh[effD*sidx + dim2], T);
								b_eq = g_eq*(Co_X[vx]*Co_X[vx] + Co_Y[vy]*Co_Y[vy] + Co_Z[vz]*Co_Z[vz] + (3-effD+K)*R*T)/2;

								// this is actually the original distribution function, recycling memory from gbar
								gbar[effD*idx + dim2] = 2*tg/(2*tg + dt/2.)*gbar[effD*idx + dim2] + dt/(4*tg + dt)*g_eq + dt*tg/(4*tg + dt)*0; //TODO replace this last *0 with source term 
								bbar[effD*idx + dim2] = 2*tb/(2*tb + dt/2.)*bbar[effD*idx + dim2] + dt/(4*tb + dt)*b_eq + dt*tb/(4*tb + dt)*0; //TODO replace this last *0 with source term
							}
						}
					}
				}
			}
		}
	}
}

//Step 2c: Compute Microflux F at interface at half timestep using W/phi at interface.
void Step2c(double* gbar, double* bbar, double* Co_X, double* Co_Y, double* Co_Z, Cell* mesh, double* Fg, double* Fb, int* BCs){
	
	if(debug == 1){printf("Entering Step 2c\n");}


	int sidx;
	int idx;
	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){

				//Spatial Index
				sidx = i + Nx*j + Nx*Ny*k;
				int pos[3] = {i,j,k};


				//Area of Interface
				double A[3] = {mesh[sidx].dy*mesh[sidx].dz, mesh[sidx].dx*mesh[sidx].dz, mesh[sidx].dx*mesh[sidx].dy};

				for(int vx = 0; vx < NV[0]; vx++){
					for(int vy = 0; vy < NV[1]; vy++){
						for(int vz = 0; vz < NV[2]; vz++){
							idx = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;
								
							Fg[idx] = 0;
							Fb[idx] = 0;

							
							for(int dim = 0; dim < effD; dim++){
								int IL, JL, KL;

								//Boundary Conditions
								double right = 1.0;
								double left = 1.0;
								//Periodic Boundary Conditions
								if(dim == 0 && BCs[0] == 0){IL = (i - 1 + N[0])%N[0]; JL = j; KL = k;} 
								if(dim == 1 && BCs[1] == 0){JL = (j - 1 + N[1])%N[1]; IL = i; KL = k;}
								if(dim == 2 && BCs[2] == 0){KL = (k - 1 + N[2])%N[2]; IL = i; JL = j;}

								//Dirichlet Boundary Conditions
								if(dim == 0 && BCs[0] == 1){IL = i - 1; JL = j; KL = k; if(IL == -1){IL = 0; left = 0.; right = 0.;} if(i == N[0]-1){left = 0.; right = 0;}} 
								if(dim == 1 && BCs[1] == 1){JL = j - 1; IL = i; KL = k; if(JL == -1){JL = 0; left = 0.; right = 0.;} if(j == N[0]-1){left = 0.; right = 0;}}
								if(dim == 2 && BCs[2] == 1){KL = k - 1; IL = i; JL = j; if(KL == -1){KL = 0; left = 0.; right = 0.;} if(k == N[0]-1){left = 0.; right = 0;}}

								//Neumann Boundary Conditions (reflective?) //TODO ensure this is actually neumann
								if(dim == 0 && BCs[0] == 2){IL = i - 1; JL = j; KL = k; if(IL == -1){IL = 0; left = 0;} if(i == N[0]-1){right = 0;}} 
								if(dim == 1 && BCs[1] == 2){JL = j - 1; IL = i; KL = k; if(JL == -1){JL = 0; left = 0;} if(j == N[0]-1){right = 0;}}
								if(dim == 2 && BCs[2] == 2){KL = k - 1; IL = i; JL = j; if(KL == -1){KL = 0; left = 0;} if(k == N[0]-1){right = 0;}}

								int idxL = IL + Nx*JL + Nx*Ny*KL + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;

								double Xi[3] = {Co_X[vx], Co_Y[vy], Co_Z[vz]};

							
								Fg[idx] += Xi[dim]*A[dim]*(right*gbar[effD*idx + dim] - left*gbar[effD*idxL + dim]);
								Fb[idx] += Xi[dim]*A[dim]*(right*bbar[effD*idx + dim] - left*bbar[effD*idxL + dim]);
								

							} 
						}
					}
				}
			}
		}
	}
}

//Step 3: Source Terms
void Step3(){


}

//Step 4: Update Conservative Variables W at cell center at next timestep
//Step 5: Update Phi at cell center at next time step
void Step4and5(double* rho, double* rhov, double* rhoE, double dt, Cell* mesh, double* Fg, double* Fb, double* Co_X, double* Co_Y, double* Co_Z, double* Co_WX, double* Co_WY, double* Co_WZ, double* g, double* b, int* BCs){
	if(debug == 1){printf("Entering Step 4 & 5\n");}

	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){
	
				int sidx = i + Nx*j + Nx*Ny*k;

				double V = mesh[sidx].dx*mesh[sidx].dy*mesh[sidx].dz;

				double rhotest = 0;
				for(int vx = 0; vx < NV[0]; vx++){
					for(int vy = 0; vy < NV[1]; vy++){
						for(int vz = 0; vz < NV[2]; vz++){
							
							

							int idx = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;
							double Xi[3] = {Co_X[vx], Co_Y[vy], Co_Z[vz]};

							//Compute Flow velocity uo and temperature To before updating rho
							double uo = 0; //old flow velocity
							for(int dim = 0; dim < effD; dim++){uo += rhov[effD*sidx + dim]/rho[sidx]*rhov[effD*sidx + dim]/rho[sidx];} uo = sqrt(uo);
							double To = Temperature(rhoE[sidx]/rho[sidx], uo); //old temp
							//Compute old taus
							double tgo = visc(To)/rho[sidx]/R/To; // tau = mu/P, P = rho*R*T.
							double tbo = tgo/Pr; 

							//Compute old eq's
							double c2 = 0;
							for(int dim = 0 ; dim < effD; dim++){ c2 += (Xi[dim]-rhov[effD*sidx + dim]/rho[sidx])*(Xi[dim]-rhov[effD*sidx + dim]/rho[sidx]);} //TODO: Potential BUG, did not double check algebra.
							double g_eqo = geq(c2, rho[sidx],To);
							double b_eqo = g_eqo*(Co_X[vx]*Co_X[vx] + Co_Y[vy]*Co_Y[vy] + Co_Z[vz]*Co_Z[vz] + (3-effD+K)*R*To)/2;


							//Step 4: Update W at cell center
							rho[sidx] += -(dt/V*Fg[idx] + dt*0)*Co_WX[vx]*Co_WY[vy]*Co_WZ[vz]; //TODO replace 0 with source term.
							for(int dim = 0; dim < effD; dim++){
								rhov[effD*sidx + dim] += -dt/V*Fg[idx]*Xi[dim]*Co_WX[vx]*Co_WY[vy]*Co_WZ[vz];
							}

							rhoE[sidx] += -dt/V*Fb[idx]*Co_WX[vx]*Co_WY[vy]*Co_WZ[vz];

							//Step 5: Update Phi at cell center, need new eq and tau
							//Compute Flow velocity u and temperature T
							double u = 0; //old flow velocity
							for(int dim = 0; dim < effD; dim++){u += rhov[effD*sidx + dim]/rho[sidx]*rhov[effD*sidx + dim]/rho[sidx];} u = sqrt(u);
							double T = Temperature(rhoE[sidx]/rho[sidx], u); //old temp
							//Compute new taus
							double tg = visc(T)/rho[sidx]/R/T; // tau = mu/P, P = rho*R*T.
							double tb = tg/Pr; 

							//printf("taus are = {%f, %f}\n", tg, tb);

							//Compute new eq's
							c2 = 0; //reset c2 from before
							for(int dim = 0 ; dim < effD; dim++){ c2 += (Xi[dim]-rhov[effD*sidx + dim]/rho[sidx])*(Xi[dim]-rhov[effD*sidx + dim]/rho[sidx]);} //TODO: Potential BUG, did not double check algebra.
							double g_eq = geq(c2, rho[sidx], T);
							double b_eq = g_eq*(Co_X[vx]*Co_X[vx] + Co_Y[vy]*Co_Y[vy] + Co_Z[vz]*Co_Z[vz] + (3-effD+K)*R*To)/2;


							//Update phi
							double testold = g[idx];

							if((BCs[0] == 1 && i == 0) || (BCs[0] == 1 && i == N[0] - 1) ||
							   (BCs[1] == 1 && j == 0) || (BCs[1] == 1 && j == N[1] - 1) ||
							   (BCs[2] == 1 && k == 0) || (BCs[2] == 1 && k == N[2] - 1))
							{
								g[idx] = g[idx];
								b[idx] = b[idx];
							}else{
								g[idx] = (g[idx] + dt/2*(g_eq/tg + (g_eqo-g[idx])/tgo) - dt/V*Fg[idx] + dt*0)/(1+dt/2/tg); //TODO replace 0 with source term
								b[idx] = (b[idx] + dt/2*(b_eq/tb + (b_eqo-b[idx])/tbo) - dt/V*Fb[idx] + dt*0)/(1+dt/2/tb); //TODO replace 0 with source term
							}
						}
					}
				}
				if(debug == 1){printf("rho[%d] = %f, rhoE[%d] = %f\n", sidx, rho[sidx], sidx, rhoE[sidx]);}
				assert(rho[sidx] == rho[sidx]); // NaN checker
			}
		}
	}
}




