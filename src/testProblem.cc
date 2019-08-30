#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "testProblem.hh"


void TestProblem(int* N, int* NV, int* Nc, int* Nv, int* BCs, double* Vmin, double* Vmax, int testProblem, double* R, double* K, double* Cv, double* gma, double* w , double* ur, double* Tr, double* Pr, int* effD){

	//Sod Shock
	if(testProblem == 1){

		//Dimensionality
		*effD = 1; 

		//Resolution
		N[0]  = 256; N[1]  = 1; N[2]  = 1;  
		NV[0] = 256; NV[1] = 1; NV[2] = 1;
		*Nc = N[0]*N[1]*N[2];
		*Nv  = NV[0]*NV[1]*NV[2];
		Vmin[0] = -10; Vmin[1] = 0; Vmin[2] = 0;
		Vmax[0] =  10; Vmax[1] = 0; Vmax[2] = 0;


		// Boundary Conditions
		BCs[0] = 1; BCs[1] = 0; BCs[2] = 0;

		
		//Physical Parameters
		*R   = 0.5;           // Gas Constant
		*K   = 2.0;           //Internal DOF
		*Cv  = (3+(*K))*(*R)/2;     //Specific Heat
		*gma = ((*K)+5)/((*K)+3);   //gamma -- variable name taken
		*w   = 0.5;           //Viscosity exponent
		*ur  = 1e-4;          //Reference Visc
		*Tr  = 1.0;           //Reference Temp
		*Pr  = 2/3.;          //Prandtl Number
	
	}

	//Kelvin-Helmholtz
	else if(testProblem == 2){

		//Dimensionality
		*effD = 2; 

		//Resolution
		int num = 64;
		N[0]  = num; N[1]  = num; N[2]  = 1;  
		NV[0] = num; NV[1] = num; NV[2] = 1;
		*Nc = N[0]*N[1]*N[2];
		*Nv  = NV[0]*NV[1]*NV[2];
		Vmin[0] = -10; Vmin[1] = -10; Vmin[2] = 0;
		Vmax[0] =  10; Vmax[1] =  10; Vmax[2] = 0;

		// Boundary Conditions
		BCs[0] = 0; BCs[1] = 0; BCs[2] = 0;

		//Physical Parameters
		*R   = 0.5;           // Gas Constant
		*K   = 2.0;           //Internal DOF
		*Cv  = (3+(*K))*(*R)/2;     //Specific Heat
		*gma = ((*K)+5)/((*K)+3);   //gamma -- variable name taken
		*w   = 0.5;           //Viscosity exponent
		*ur  = 1e-4;          //Reference Visc
		*Tr  = 1.0;           //Reference Temp
		*Pr  = 2/3.;          //Prandtl Number
	}
}


void InitializeTestProblem(Cell* mesh, double* g, double* b, double* rho, double* rhov, double* rhoE, int testProblem, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double R, double K, double Cv, double gma, double w, double ur, double Tr, double Pr, int* N, int* NV, int effD){



	if (testProblem == 1){
		SodShock(mesh, g, b, rho, rhov, rhoE, Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ, R, K, Cv, gma, w, ur, Tr, Pr, N, NV, effD);
	}

	else if (testProblem == 2){
		KHI(mesh, g, b, rho, rhov, rhoE, Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ, R, K, Cv, gma, w, ur, Tr, Pr, N, NV, effD); 

	}

	else if (testProblem == 3){
		RTI(); //TODO: initialize
	}
}


void SodShock(Cell* mesh, double* g, double* b, double* rho, double* rhov, double* rhoE, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double R, double K, double Cv, double gma, double w, double ur, double Tr, double Pr, int* N, int* NV, int effD, double rhoL, double rhoR, double PL, double PR){
	
	int idx;
	int Nx = N[0];
	int Ny = N[1];
	int Nz = N[2];
	int Nvx = NV[0];
	int Nvy = NV[1];
	int Nvz = NV[2];


	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){

				idx = i + Nx*j + Nx*Ny*k; //spatial index

				//Left State
				if(mesh[idx].x <= 0.5){

					//Conserved Variables
					rho[idx] = rhoL; 
					rhoE[idx] = 0.0; //Initialize for Reduction
					for(int dim = 0; dim < effD; dim++){
						rhov[effD*idx + dim] = 0.0;
						rhoE[idx] += 0.5*rhov[effD*idx+dim]*rhov[effD*idx+dim]/rho[idx]; // 0.5* rhov**2/rho
					}
					rhoE[idx] += Cv*PL/R; // T = P/(rho*R), Ideal Gas 
					

				}

				//Right State
				if(mesh[idx].x > 0.5){

					//Conserved Variables
					rho[idx] = rhoR; 
					rhoE[idx] = 0.0; //Initialize for Reduction
					for(int dim = 0; dim < effD; dim++){
						rhov[effD*idx + dim] = 0.0;
						rhoE[idx] += 0.5*rhov[effD*idx+dim]*rhov[effD*idx+dim]/rho[idx]; // 0.5* rhov**2/rho
					}
					rhoE[idx] += Cv*PR/R; // T = P/(rho*R), Ideal Gas 

				}
				printf("Initialized W[%d]: rho = %f, rhov = %f, rhoE = %f\n", idx, rho[idx], rhov[idx], rhoE[idx]);

			}
		}
	}
				

	//Initialize g and b
	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){

				int sidx = i + Nx*j + Nx*Ny*k; //spatial index

				double u = 0;
				for(int dim = 0; dim < effD; dim++){u += rhov[effD*sidx + dim]/rho[sidx]*rhov[effD*sidx + dim]/rho[sidx];}
				u = sqrt(u);

				double T = Temperature(rhoE[sidx]/rho[sidx], u);

				//double rhotest = 0;
				//double rhovxtest = 0;
				//double rhoEtest = 0;

				for(int vx = 0; vx < NV[0]; vx++){
					for(int vy = 0; vy < NV[1]; vy++){
						for(int vz = 0; vz < NV[2]; vz++){

						
						idx = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;


						double c2 = 0;
						double Xi[3] = {Co_X[vx], Co_Y[vy], Co_Z[vz]};

						for(int dim = 0 ; dim < effD; dim++){ c2 += (Xi[dim]-rhov[effD*sidx + dim]/rho[sidx])*(Xi[dim]-rhov[effD*sidx + dim]/rho[sidx]);} //TODO: Potential BUG, did not double check algebra.
						g[idx] = geq(c2, rho[sidx], T);
						b[idx] = g[idx]*(Co_X[vx]*Co_X[vx] + Co_Y[vy]*Co_Y[vy] + Co_Z[vz]*Co_Z[vz] + (3-effD+K)*R*T)/2;


						//rhotest += g[idx]*Co_WX[vx]*Co_WY[vy]*Co_WZ[vz];//;
						//rhovxtest += g[idx]*Co_X[vx]*Co_WX[vx]*Co_WY[vy]*Co_WZ[vz];
						//rhoEtest += b[idx]*Co_WX[vx]*Co_WY[vy]*Co_WZ[vz];
						//printf("Initial g[%d] = %f\n", idx, g[idx]);
						//printf("Initial b[%d] = %f\n", idx, b[idx]);

						}
					}
				}

				//printf("Initialization W[%d] = {%f, %f, %f}\n", sidx, rhotest, rhovxtest, rhoEtest);
			}

		}

	}



		//printf("rho[idx] = %f, rhov[idx] = %f, rhoE[idx] = %f, idx = %d\n", rho[idx], rhov[D*idx], rhoE[idx], idx);


}

 void KHI(Cell* mesh, double* g, double* b, double* rho, double* rhov, double* rhoE, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double R, double K, double Cv, double gma, double w, double ur, double Tr, double Pr, int* N, int* NV, int effD, double rhoT, double rhoB, double PT, double PB, double vrel, double amp){

	int idx;
	int Nx = N[0];
	int Ny = N[1];
	int Nz = N[2];
	int Nvx = NV[0];
	int Nvy = NV[1];
	int Nvz = NV[2];

	double PI = 4*atan(1.0);

	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){

				idx = i + Nx*j + Nx*Ny*k; //spatial index

				//Bottom State
				if(mesh[idx].y <= 0.5){

					//Conserved Variables
					rho[idx] = rhoT; 
					rhoE[idx] = 0.0; //Initialize for Reduction
					for(int dim = 0; dim < effD; dim++){
						if(dim == 0){
							rhov[effD*idx + dim] = vrel/2.*rho[idx];
						}
						else if(dim == 1){
							rhov[effD*idx + dim] = amp*sin(2*PI*mesh[idx].x)*rho[idx];
						}
						rhoE[idx] += 0.5*rhov[effD*idx + dim]*rhov[effD*idx+dim]/rho[idx]; // 0.5* rhov**2/rho
					}
					rhoE[idx] += Cv*PT/R; // T = P/(rho*R), Ideal Gas 
					

				}

				//Right State
				if(mesh[idx].y > 0.5){

					//Conserved Variables
					rho[idx] = rhoB; 
					rhoE[idx] = 0.0; //Initialize for Reduction
					for(int dim = 0; dim < effD; dim++){
						if(dim == 0){
							rhov[effD*idx + dim] = -vrel/2.*rho[idx];
						}
                                                else if(dim == 1){
                                                        rhov[effD*idx + dim] = amp*sin(2*PI*mesh[idx].x)*rho[idx];
                                                }
				
						rhoE[idx] += 0.5*rhov[effD*idx+dim]*rhov[effD*idx+dim]/rho[idx]; // 0.5* rhov**2/rho
					}
					rhoE[idx] += Cv*PB/R; // T = P/(rho*R), Ideal Gas 

				}

				if(i == Nx/2) printf("Initialized central column, W[%d]: rho = %f, rhovx = %f, rhovy %f, rhoE = %f\n", idx, rho[idx], rhov[effD*idx + 0], rhov[effD*idx + 1], rhoE[idx]);

			}
		}
	}
				

	//Initialize g and b
	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){

				int sidx = i + Nx*j + Nx*Ny*k; //spatial index

				double u = 0;
				for(int dim = 0; dim < effD; dim++){u += rhov[effD*sidx + dim]/rho[sidx]*rhov[effD*sidx + dim]/rho[sidx];}
				u = sqrt(u);

				double T = Temperature(rhoE[sidx]/rho[sidx], u);

				//double rhotest = 0;
				//double rhovxtest = 0;
				//double rhoEtest = 0;

				for(int vx = 0; vx < NV[0]; vx++){
					for(int vy = 0; vy < NV[1]; vy++){
						for(int vz = 0; vz < NV[2]; vz++){

						
						idx = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;


						double c2 = 0;
						double Xi[3] = {Co_X[vx], Co_Y[vy], Co_Z[vz]};

						for(int dim = 0 ; dim < effD; dim++){ c2 += (Xi[dim]-rhov[effD*sidx + dim]/rho[sidx])*(Xi[dim]-rhov[effD*sidx + dim]/rho[sidx]);} //TODO: Potential BUG, did not double check algebra.
						g[idx] = geq(c2, rho[sidx], T);
						b[idx] = g[idx]*(Co_X[vx]*Co_X[vx] + Co_Y[vy]*Co_Y[vy] + Co_Z[vz]*Co_Z[vz] + (3-effD+K)*R*T)/2;


						//rhotest += g[idx]*Co_WX[vx]*Co_WY[vy]*Co_WZ[vz];//;
						//rhovxtest += g[idx]*Co_X[vx]*Co_WX[vx]*Co_WY[vy]*Co_WZ[vz];
						//rhoEtest += b[idx]*Co_WX[vx]*Co_WY[vy]*Co_WZ[vz];
						//printf("Initial g[%d] = %f\n", idx, g[idx]);
						//printf("Initial b[%d] = %f\n", idx, b[idx]);

						}
					}
				}

				//printf("Initialization W[%d] = {%f, %f, %f}\n", sidx, rhotest, rhovxtest, rhoEtest);
			}

		}

	}



		//printf("rho[idx] = %f, rhov[idx] = %f, rhoE[idx] = %f, idx = %d\n", rho[idx], rhov[D*idx], rhoE[idx], idx);


 	
 }

 void RTI(){

 	//TODO
 }

