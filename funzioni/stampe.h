
typedef double Real;
#include <fstream>
void stampamat(Real *a, int n, int m, std::ofstream *prt = 0);
void stampa(int d, Real t, Real *u, void (*sol)(Real *, Real), std::ofstream *prt = 0);
void stampafloat(int d, float t, float *u, std::ofstream *prt = 0);
void stampa(int d, Real t, long double *u, std::ofstream *prt = 0);
