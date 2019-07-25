#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "testProblem.hh"
#include "Functions.hh"

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
			for(int k = 0; k < N[1]; k++){

				idx = i + Nx*j + Nx*Ny*k;

				//Left State
				if(mesh[idx].x <= 0.5){

					//Conserved Variables
					rho[idx] = rhoL; 
					rhoE[idx] = 0.0; //Initialize for Reduction
					for(int dim = 0; dim < D; dim++){
						rhov[D*idx + dim] = 0.0;
						rhoE[idx] += 0.5*rhov[D*idx+dim]*rhov[D*idx+dim]/rho[idx]; // 0.5* rhov**2/rho
					}
					rhoE[idx] += Cv*PL/R; // T = P/(rho*R), Ideal Gas 
					

				}

				//Right State
				if(mesh[idx].x > 0.5){

					//Conserved Variables
					rho[idx] = rhoR; 
					rhoE[idx] = 0.0; //Initialize for Reduction
					for(int dim = 0; dim < D; dim++){
						rhov[D*idx + dim] = 0.0;
						rhoE[idx] += 0.5*rhov[D*idx+dim]*rhov[D*idx+dim]/rho[idx]; // 0.5* rhov**2/rho
					}
					rhoE[idx] += Cv*PR/R; // T = P/(rho*R), Ideal Gas 

				}


				//Initialize g and b
				for(int i = 0; i < N[0]; i++){
					for(int j = 0; j < N[1]; j++){
						for(int k = 0; k < N[2]; k++){
							for(int vx = 0; vx < NV[0]; vx++){
								for(int vy = 0; vy < NV[1]; vy++){
									for(int vz = 0; vz < NV[2]; vz++){

									int sidx = i + Nx*j + Nx*Ny*k; //spatial index
									idx = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*vx + Nx*Ny*Nz*NV[0]*vy + Nx*Ny*Nz*NV[0]*NV[1]*vz;
									//printf("idx = %d", idx);

									double U[3] = {rhov[D*sidx]/rho[sidx], rhov[D*sidx + 1]/rho[sidx], rhov[D*sidx + 2]/rho[sidx]};
									double T = Temperature(rhoE[sidx]/rho[sidx], sqrt(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]));

									double effD = 1.0; //TODO: POTENTIAL BUG
									g[idx] = geq(vx, vy, vz, rho[sidx], U, T, Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ);
									b[idx] = g[idx]*(Co_X[vx]*Co_X[vx] + Co_Y[vy]*Co_Y[vy] + Co_Z[vz]*Co_Z[vz] + (3-effD+K)*R*T)/2;
									}
								}
							}
						}

					}

				}


			}
		}

		//printf("rho[idx] = %f, rhov[idx] = %f, rhoE[idx] = %f, idx = %d\n", rho[idx], rhov[D*idx], rhoE[idx], idx);

	}

}

 void KHI(){

 	//TODO
 }

 void RTI(){

 	//TODO
 }

