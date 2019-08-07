#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH


double Temperature(double E, double u);

double geq(double vx, double vy, double vz, double rho, double* U, double T, double wx, double wy, double wz);

double visc(double T);

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
double VanLeer(double L, double C, double R, double xL, double xC, double xR);



#endif