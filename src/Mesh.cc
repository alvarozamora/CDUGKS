#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "Mesh.hh"

// This function generates the mesh structure. 
// Computes cell centers and cell sizes.
//
//

void Mesh(int* N, Cell* mesh, int MeshType){
	//User-Specified
	if (MeshType == 0){ 
		//TODO
	}
	//Rectangular
	else if (MeshType == 1){

		int Nx = N[0];
		int Ny = N[1];
		int Nz = N[2];

		double dX = 1.0/Nx;
		double dY = 1.0/Ny;
		double dZ = 1.0/Nz;

		int idx;
		for(int i = 0; i < Nx; i++){
			for(int j = 0; j < Ny; j++){
				for(int k = 0; k < Nz; k++){

					idx = i + Nx*j + Nx*Ny*k;

					mesh[idx].x = dX*(i + 1.0/2.0);
					mesh[idx].y = dY*(j + 1.0/2.0);
					mesh[idx].z = dZ*(k + 1.0/2.0);

					mesh[idx].dx = dX;
					mesh[idx].dy = dY;
					mesh[idx].dz = dZ;

				}
			}
		}

	}
	//Nested
	else if(MeshType == 2){
		//TODO
	}
}