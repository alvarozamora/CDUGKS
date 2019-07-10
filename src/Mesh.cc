#include <iostream>
#define D 1
#define B 1
#define Nv 128

struct cell
{
	double X[D];    //Position of Cell Center
	double dX[D];   //Cell Size
	double I[B][D]; //Interface Area Vectors
};

int mesh()
{
	if(UserDefinedMesh){


	}
	else if(UniformMesh){

	int N = 128;




	std::cout << "a = "<<a << std::endl;


	}
	else if(NestedMesh){
		
	}
}