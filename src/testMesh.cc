#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Mesh.hh"

#define D 3
#define Nv

int N[3] = {128,4,4};

int main(){

	int MeshType = 1; // 0 is UserDefinedMesh, 1 is RectangularMesh, 2 is Nested Rectangular Mesh.

	printf("Try Declaration\n");
	struct Cell mesh[N[0]*N[1]*N[2]];
	printf("Confirm Declaration\n");

	printf("Try Mesh\n");
	Mesh(N, mesh, MeshType);
	printf("Confirm Mesh\n");
}