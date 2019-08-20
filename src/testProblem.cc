#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "testProblem.hh"


void TestProblem(Cell* mesh, double* g, double* b, double* rho, double* rhov, double* rhoE, int testProblem, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ){

	if (testProblem == 1){
		SodShock(mesh, g, b, rho, rhov, rhoE, Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ);
	}

	else if (testProblem == 2){
		KHI(); //TODO: initialize
	}

	else if (testProblem == 3){
		RTI(); //TODO: initialize
	}
}


void SodShock(Cell* mesh,  double* g, double* b, double* rho, double* rhov, double* rhoE, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double rhoL, double rhoR, double PL, double PR){
	
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

				idx = i + Nx*j + Nx*Ny*k;

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
				//printf("idx = %d, rho[idx] = %f\n", idx, rho[idx]);

			}
		}
	}
				

	//Initialize g and b
	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){

				int sidx = i + Nx*j + Nx*Ny*k; //spatial index

				double U[3] = {rhov[effD*sidx]/rho[sidx], rhov[effD*sidx + 1]/rho[sidx], rhov[effD*sidx + 2]/rho[sidx]};

				double T = Temperature(rhoE[sidx]/rho[sidx], sqrt(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]));


				for(int vx = 0; vx < NV[0]; vx++){
					for(int vy = 0; vy < NV[1]; vy++){
						for(int vz = 0; vz < NV[2]; vz++){

						
						idx = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;


						g[idx] = geq(Co_X[vx], Co_Y[vy], Co_Z[vz], rho[sidx], U, T, Co_WX[vx], Co_WY[vy], Co_WZ[vz]);
						b[idx] = g[idx]*(Co_X[vx]*Co_X[vx] + Co_Y[vy]*Co_Y[vy] + Co_Z[vz]*Co_Z[vz] + (3-effD+K)*R*T)/2;


						//printf("Initial g[%d] = %f\n", idx, g[idx]);
						//printf("Initial b[%d] = %f\n", idx, b[idx]);

						}
					}
				}
			}

		}

	}



		//printf("rho[idx] = %f, rhov[idx] = %f, rhoE[idx] = %f, idx = %d\n", rho[idx], rhov[D*idx], rhoE[idx], idx);


}

 void KHI(){

 	//TODO
 }

 void RTI(){

 	//TODO
 }

