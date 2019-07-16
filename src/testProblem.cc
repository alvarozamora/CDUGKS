#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "testProblem.hh"


void TestProblem(Cell* mesh, double* g, double* b, double* rho, double* rhov, double* rhoE, int testProblem){

	if (testProblem == 1){
		SodShock(mesh, g, b, rho, rhov, rhoE); //TODO, set D=1 inside, initialize
	}

	else if (testProblem == 2){
		KHI(); //TODO, set D=2 inside, initialize
	}

	else if (testProblem == 3){
		RTI(); //TODO, set D=2 inside, initialize
	}
}


void SodShock(Cell* mesh,  double* g, double* b, double* rho, double* rhov, double* rhoE, double rhoL, double rhoR, double PL, double PR){
	
	int idx;
	int Nx = N[0];
	int Ny = N[1];
	int Nz = N[2];
	int Nvx = NV[0];
	int Nvy = NV[1];
	int Nvz = NV[2];


	printf("Cv = %f\n", Cv);

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

				//Rights State
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


			}
		}

	}

}

 void KHI(){

 	//TODO
 }

 void RTI(){

 	//TODO
 }

