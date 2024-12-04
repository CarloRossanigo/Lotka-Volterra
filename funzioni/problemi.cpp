#include <math.h>
#include <iostream>
#include "problemi.h"
using namespace std;
typedef double Real;
const Real pi=4.0*atan(1.0);
const Real A=1.0;
const Real C=A;
const Real B=0.1;
const Real D=0.2;
const Real G=0.3;
// ********************************************************
//
//
////


void F_lotkav(Real *f,Real t,Real *u)
{

    f[0]=A*u[0]-B*u[0]*u[1]-G*u[0];
    f[1]=-C*u[1]+D*u[0]*u[1]-G*u[1];
}


void jeffe_lotkav(Real *jf,Real t,Real *u)
{
    jf[0]=A-B*u[1]-G;
    jf[1]=-B*u[0];
    jf[2]=+D*u[1];
    jf[3]=-C+D*u[0]-G;
}
void dati_lotkav(Real *v,Real *T,Real *U)
{
    *v=0.0;
    *T=20.0;
    U[0]=16.0;
    U[1]=12.0;


}


