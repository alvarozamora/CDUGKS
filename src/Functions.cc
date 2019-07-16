#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Functions.hh"


double Temperature(double E, double u){

	double T = (gamma - 1)/R * (E - 0.5*u*u);

	return T
}