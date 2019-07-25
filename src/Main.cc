#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Main.hh"
#include "Mesh.hh"
#include "testProblem.hh"
#include "Functions.hh"
#include "Evolution.hh"

int main(){

	//Declare Physical Quantities
	double g[Nc*Nv];   //g and b are reduced distrubution functions (Vel and E distribution)
	double b[Nc*Nv];

	double rho[Nc];    //Conserved Variables
	double rhov[Nc*D];
	double rhoE[Nc];

	//Newton-Cotes Quadrature
	double Co_X[NV[0]]; double Co_WX[NV[0]]; double Co_Y[NV[1]]; double Co_WY[NV[1]]; double Co_Z[NV[2]]; double Co_WZ[NV[2]]; // Cotes points and weights
	Cotes(Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ);

	//Generate Mesh: Grid Cell Centers and Sizes
	int MeshType = 1; // 0 is UserDefinedMesh, 1 is RectangularMesh, 2 is Nested Rectangular Mesh.
	struct Cell mesh[N[0]*N[1]*N[2]];
	Mesh(N, mesh, MeshType);


	//Initialize Grid
	int testProblem = 1; //0 is None, 1 is Sod Shock, 2 is KHI, 3 is RTI.
	if (testProblem > 0){TestProblem(mesh, g, b, rho, rhov, rhoE, testProblem, Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ);}
	else{
		//TODO
	}

	//Evolve
	double Tsim = 0.;
	double dt = 128.;
	double Tf = 0.15;

	while(Tsim < Tf){
		Evolve();

		Tsim += dt;
	}

}