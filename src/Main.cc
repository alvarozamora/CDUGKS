#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Main.hh"
#include "testProblem.hh"
#include "Functions.hh"
#include "Evolution.hh"

int main(){

	//Declare Physical Quantities
	printf("Declaring Variables\n");
	double g[Nc*Nv];   //g and b are reduced distrubution functions (Vel and E distribution)
	double b[Nc*Nv];
	double gbarp[Nc*Nv];   //gbarp and bbarp are reduced distrubution functions (Vel and E distribution)
	double bbarp[Nc*Nv];
	//At interface
	double gbar[Nc*Nv*effD];   //gbar and bbar are reduced distrubution functions (Vel and E distribution)
	double bbar[Nc*Nv*effD];
	double gbarpbound[Nc*Nv*effD];
	double bbarpbound[Nc*Nv*effD];

	double Sg[Nc]; //Source Terms
	double Sb[Nc];

	double rho[Nc];    //Conserved Variables at t
	double rhov[Nc*effD];
	double rhoE[Nc];
	double rhoh[Nc*effD];   //Conserved Variables at t + h, at interfaces
	double rhovh[Nc*effD*effD];
	double rhoEh[Nc*effD];

	double gsigma[Nc*Nv*effD]; //Gradients
	double bsigma[Nc*Nv*effD];
	double gsigma2[Nc*Nv*effD*effD];
	double bsigma2[Nc*Nv*effD*effD];



	//Newton-Cotes Quadrature
	printf("Setting up NC-Quadrature\n");
	double Co_X[NV[0]];   // Cotes points and weights
	double Co_WX[NV[0]]; 
	double Co_Y[NV[1]]; 
	double Co_WY[NV[1]]; 
	double Co_Z[NV[2]]; 
	double Co_WZ[NV[2]]; 
	Cotes(Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ);

	//Checking NC Weights on 128-cell Sod Problem
	//for(int i = 0; i < 128; i++){printf("Main.cc Co_X[%d] = %f\n", i, Co_X[i]);}
	//for(int i = 0; i < 128; i++){printf("Main.cc Co_WX[%d] = %f\n", i, Co_WX[i]);}

	//Generate Mesh: Grid Cell Centers and Sizes
	printf("Generating Mesh\n");
	int MeshType = 1; // 0 is UserDefinedMesh, 1 is RectangularMesh, 2 is Nested Rectangular Mesh.
	struct Cell mesh[N[0]*N[1]*N[2]];
	Mesh(N, mesh, MeshType);


	//Initialize Grid
	printf("Initializing Grid on Mesh\n");
	int testProblem = 1; //0 is None, 1 is Sod Shock, 2 is KHI, 3 is RTI.
	if (testProblem > 0){TestProblem(mesh, g, b, rho, rhov, rhoE, testProblem, Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ);}
	else{
		//TODO
	}

	//Evolve
	double Tsim = 0.;
	double dt = 128.;
	double Tf = 0.15;
	double dtdump = Tf/300;

	printf("Entering Evolution Loop\n");
	while(Tsim < Tf){
		Evolve(g, b, gbar, bbar, gbarp, bbarp, Sg, Sb, rho, rhov, rhoE, effD, dt, Tf, Tsim, dtdump, Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ, gsigma, bsigma, gsigma2, bsigma2, mesh, gbarpbound, bbarpbound,rhoh, rhovh, rhoEh);

		Tsim += dt;
	}

}