#include <iostream>
#include <fstream>
#include "stampe.h"
using namespace std;
typedef double Real;
void stampamat(Real *a, int n, int m, std::ofstream *prt)
{
    if (prt == 0)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                cout << *(a + j + i * m) << "\t";
            }
            cout << "\n";
        }
        cout << " " << endl;
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                *prt << *(a + j + i * m) << "\t";
            }
            *prt << "\n";
        }
        *prt << " " << endl;
    }
}

void stampa(int d, Real t, Real *u, void (*sol)(Real *, Real), std::ofstream *prt)
{

    if (prt == 0)
    {

        cout << t << "\t";
        for (int i = 0; i < d; i++)
            cout << u[i] << "\t";
        if (sol != 0)
        {
            Real s[d];
            sol(s, t);
            for (int i = 0; i < d; i++)
                cout << s[i] << "\t";
        }

        cout << endl;
    }
    else
    {

        *prt << t << "\t";
        for (int i = 0; i < d; i++)
            *prt << u[i] << "\t";
        if (sol != 0)
        {
            Real s[d];
            sol(s, t);
            for (int i = 0; i < d; i++)
                *prt << s[i] << "\t";
        }
        *prt << endl;
    }
}

void stampafloat(int d, float t, float *u, std::ofstream *prt)
{
    if (prt == 0)
    { // stampa a video
        cout << t << "\t";
        for (int i = 0; i < d; i++)
        {
            cout << u[i] << "\t";
        }
        cout << endl;
    }
    else
    { // stampa su file al puntatore
        *prt << t << "\t";
        for (int i = 0; i < d; i++)
        {
            *prt << u[i] << "\t";
        }
        *prt << endl;
    }
}

//! long double
void stampa(int d, Real t, long double *u, std::ofstream *prt)
{
    if (prt == 0)
    {
        cout << t << "\t";
        for (int i = 0; i < d; i++)
        {
            cout << u[i] << "\t";
        }
        cout << endl;
    }
    else
    {
        *prt << t << "\t";
        for (int i = 0; i < d; i++)
        {
            *prt << u[i] << "\t";
        }
        *prt << endl;
    }
}
