#ifndef MESH_HH
#define MESH_HH

extern int D;
extern int Nv;

struct Cell{

	double x;
	double y;
	double z;

	double dx;
	double dy;
	double dz;
};



void Mesh(int* N, Cell* mesh, int MeshType);
#endif
