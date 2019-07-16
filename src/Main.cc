#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Main.hh"
#include "Mesh.hh"
#include "testProblem.hh"


int main(){

	//Declare Physical Quantities
	double g[Nc*Nv];   //g and b are reduced distrubution functions (Vel and E distribution)
	double b[Nc*Nv];

	double rho[Nc];    //Conserved Variables
	double rhov[Nc*D];
	double rhoE[Nc];


	//Generate Mesh: Grid Cell Centers and Sizes
	int MeshType = 1; // 0 is UserDefinedMesh, 1 is RectangularMesh, 2 is Nested Rectangular Mesh.
	struct Cell mesh[N[0]*N[1]*N[2]];
	Mesh(N, mesh, MeshType);


	//Initialize Grid
	int testProblem = 0; //0 is None, 1 is Sod Shock, 2 is KHI, 3 is RTI.
	if (testProblem > 0){TestProblem(mesh, g, b, rho, rhov, rhoE, testProblem);}
	else{
		//TODO
	}






}