#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "Main.hh"
#include "testProblem.hh"
#include "Functions.hh"
#include "Evolution.hh"

int main(){

	int testProblem = 2; //0 is None, 1 is Sod Shock, 2 is KHI, 3 is RTI.
	if (testProblem > 0){TestProblem(N, NV, &Nc, &Nv, BCs, Vmin, Vmax, testProblem, &R, &K, &Cv, &gma, &w , &ur, &Tr, &Pr, &effD);}
	else{
		//TODO Non Test Problems
	}


	printf("N = {%d, %d, %d}, effD = %d\n", N[0], N[1], N[2], effD);
	printf("NV = {%d, %d, %d}\n", NV[0], NV[1], NV[2]);
	printf("R = %f,  K = %f, Cv = %f,  g = %f\n", R, K, Cv, gma);
	printf("w = %f, ur = %f, Tr = %f, Pr = %f\n", w, ur, Tr, Pr);
	printf("Confirm Parameters, Enter/Return to Continue:\n");
	getchar();


	int numdoub = 0;
	numdoub += Nc*Nv;   //g and b are reduced distrubution functions (Vel and E distribution)
        numdoub += Nc*Nv;
        numdoub += Nc*Nv;   //gbarp and bbarp are reduced distrubution functions (Vel and E distribution)
        numdoub += Nc*Nv;

	numdoub += Nc*Nv*effD; // gbar/p and bbar/p are reduced distrubution functions (Vel and E distribution)
        numdoub += Nc*Nv*effD;
        numdoub += Nc*Nv*effD;
        numdoub += Nc*Nv*effD;


        numdoub += Nc; //Source Terms
        numdoub += Nc*Nv;

        numdoub += Nc;    //Conserved Variables at t
        numdoub += Nc*effD;
        numdoub += Nc;
        numdoub += Nc*effD;  //Conserved Variables at t + h, at interfaces
        numdoub += Nc*effD*effD;
        numdoub += Nc*effD;

        numdoub += Nc*Nv*effD;//Gradients
        numdoub += Nc*Nv*effD;
        numdoub += Nc*Nv*effD*effD;
        numdoub += Nc*Nv*effD*effD;

	printf("Total Number of Doubles = %d\n", numdoub);
	getchar();

	//Declare Physical Quantities
	printf("Declaring Variables\n");
	double* g = new double[Nc*Nv];   //g and b are reduced distrubution functions (Vel and E distribution)
	double* b = new double[Nc*Nv];
	double* gbarp = new double[Nc*Nv];   //gbarp and bbarp are reduced distrubution functions (Vel and E distribution)
	double* bbarp = new double[Nc*Nv];
	printf("Declared Reduced Distribution Functions\n");

	//At interface
	double* gbarpbound = new double[Nc*Nv*effD]; // gbar/p and bbar/p are reduced distrubution functions (Vel and E distribution)
	double* bbarpbound = new double[Nc*Nv*effD];
	double* gbar = new double[Nc*Nv*effD];
	double* bbar = new double[Nc*Nv*effD];
	printf("Declared Interface Reduced Distribution Functions\n");

	double* Sg = new double[Nc]; //Source Terms
	double* Sb = new double[Nc];
	printf("Declared Soure Terms\n");

	double* rho = new double[Nc];   //Conserved Variables at t
	double* rhov = new double[Nc*effD];
	double* rhoE = new double[Nc];
	double* rhoh = new double[Nc*effD]; //Conserved Variables at t + h, at interfaces
	double* rhovh = new double[Nc*effD*effD];
	double* rhoEh = new double[Nc*effD];
	printf("Declared Conserved Variable Arrays\n");

	double* gsigma = new double[Nc*Nv*effD]; // Gradients
	double* bsigma = new double[Nc*Nv*effD];
	printf("Declared Gradient Arrays\n");
	double* gsigma2 = new double[Nc*Nv*effD*effD];
	double* bsigma2 = new double[Nc*Nv*effD*effD];
	printf("Declared Second Gradient Arrays\n");


	//Flux
	//double Fg[Nc*Nv];
	//double Fb[Nc*Nv];

	//Newton-Cotes Quadrature
	printf("Setting up NC-Quadrature\n");
	double* Co_X = new double[NV[0]];   // Cotes points and weights
	double* Co_WX = new double[NV[0]];
	double* Co_Y = new double[NV[1]];
	double* Co_WY = new double[NV[1]];
	double* Co_Z = new double[NV[2]];
	double* Co_WZ = new double[NV[2]];
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
	if (testProblem > 0){

		InitializeTestProblem(mesh, g, b, rho, rhov, rhoE, testProblem, Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ, R, K, Cv, gma, w, ur, Tr, Pr, N, NV, effD);

		printf("Confirm Initial Conditions, Enter/Return to Continue:\n");
		getchar();
	}
	else{
		//TODO
	}
	
	

	//Evolve
	printf("Declaring Time Variables\n");
	double* Tsim = new double(0.);
	double* dt = new double;
	double* Tf = new double(0.15);
	double* Tdump = new double(0.0);
	double* dtdump = new double(*Tf/200.);
	printf("Declared Time Variables\n");

	if(testProblem == 1){
		*Tf = 0.15;
		*dtdump = *Tf/200.;
	}
	if(testProblem == 2){
		*Tf = 1.2;
		*Tf = 2.0;
		*dtdump = *Tf/400.;
	}



	int iter = 0;
	int dumpiter = 0;
	datadeal(mesh, rho, dumpiter, testProblem);
	printf("Entering Evolution Loop\n");
	while(*Tsim < *Tf){ // && iter < itermax
		iter++;

		int dump = Evolve(g, b, gbar, bbar, gbarp, bbarp, Sg, Sb, rho, rhov, rhoE, dt, *Tf, *Tsim, *dtdump, Co_X, Co_WX, Co_Y, Co_WY, Co_Z, Co_WZ, gsigma, bsigma, gsigma2, bsigma2, mesh, gbarpbound, bbarpbound, rhoh, rhovh, rhoEh, Tdump, BCs, R, K, Cv, gma, w, ur, Tr, Pr, N, NV, effD, Vmax);

		*Tsim += *dt;
		*Tdump += *dt;

		if(dump == 1){
			*Tdump = 0.0;
			dumpiter++;
			datadeal(mesh, rho, dumpiter, testProblem);
		}
		printf("iteration = %d, timestep = %f, Tsim = %f, Tdump = %f, dtdump = %f, dumpiter = %d\n", iter, *dt, *Tsim, *Tdump, *dtdump, dumpiter);
	}

	//show data
	for(int i = 0; i < N[0]; i++){
		for(int j = 0; j < N[1]; j++){
			for(int k = 0; k < N[2]; k++){
				int sidx = i + N[0]*j + N[0]*N[1]*k;
				printf("rho[%d] = %1.10f\n", sidx, rho[sidx]);
			}
		}
	}


}



void datadeal(Cell* mesh, double* rho, int iter, int testProblem){

	FILE *fp;

	if (iter == 0){
		fp=fopen("Data/x.txt","w");
		for (int i=0; i<N[0]*N[1]*N[2]; i++) fprintf(fp,"%e\n", mesh[i].x);
		fclose(fp);

		if(testProblem > 0){
			fp=fopen("Data/index.txt","w");
			fprintf(fp, "%d", testProblem);
			fclose(fp);
		}
	}

	// Allocates storage
	char *rhofile = (char*)malloc(17 * sizeof(char));
	sprintf(rhofile, "Data/rho%04d.txt", iter);

	printf(rhofile); printf("\n");
	fp=fopen(rhofile ,"w");
	for (int i=0; i<N[0]*N[1]*N[2]; i++) fprintf(fp,"%f\n",rho[i]);
	fclose(fp);

	

	/*

	fp=fopen("u.dat","w");
	for (i=imin; i<=imax; i++) fprintf(fp,"%e ",ux[i]);
	fclose(fp);

	fp=fopen("T.dat","w");
	for (i=imin; i<=imax; i++) fprintf(fp,"%e ",T[i]);
	fclose(fp);
	*/
}
