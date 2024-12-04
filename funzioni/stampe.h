using namespace std;
typedef double Real;
#include <fstream>
void stampamat(Real *a,int n,int m,ofstream *prt=0);
void stampa(int d,Real t,Real *u,void (*sol)(Real *,Real),ofstream *prt=0);
void stampafloat(int d,float t,float *u,ofstream *prt=0);
void stampa(int d,Real t,long double *u,ofstream *prt=0);
