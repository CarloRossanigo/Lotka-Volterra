##include <iostream>
#include <math.h>
#include <fstream>
//using MacOSÉ if you use another OS  you may change "/" to "\"
#include "../funzioni/problemi.h" //#include "..\funzioni\problemi.h"
#include "../funzioni/mat_vet.h" // #include "..\funzioni\mat_vet.h"
#include "../funzioni/Runge_Kutta.h" //#include "..\funzioni\Runge_Kutta.h"
#include "../funzioni/stampe.h" // #include "..\funzioni\stampe.h"

using namespace std;
typedef double Real;
const Real pi = 4.0 * atan(1.0); // Definizione di pi greco (Pi definition)
// Puntatori per la funzione del campo vettoriale e altre operazioni (Pointers for vector field function and other operations)
void (*effe)(Real*, Real, Real*);                           // funzione del campo vettoriale dell'equazione differenziale (Vector field function of differential equation)
void (*dati)(Real*, Real*, Real*);                          // dati iniziali del problema (Initial data for the problem)
void (*butcher)(int, Real*, Real*, Real*, Real*, int*);    // coefficienti della tabella butcher (Butcher table coefficients)
void (*sol)(Real* y, Real t);                               // funzione soluzione (Solution function)

// Variabili comuni (Common variables)
Real t, T, h, hnew;                                        // variabili temporali (time variables)
unsigned long valutazioni_f;                                // variabile per il conteggio delle valutazioni di f (variable for counting evaluations of f)
const int dim = 2;                                          // dimensione del problema (dimension of the problem)
const Real p = 8.0;                                         // ordine del metodo (method order)
Real u[dim];                                                // variabile per la soluzione (solution variable)

// Impostazione della tabella di Butcher (Setting up the Butcher table)
const int ns = 13;  // numero di stadi (number of stages)
Real b[ns];         // coefficienti b della tabella di Butcher (coefficients b of the Butcher table)
Real c[ns];         // coefficienti c della tabella di Butcher (coefficients c of the Butcher table)
Real bc[ns];        // coefficienti bc della tabella di Butcher (coefficients bc of the Butcher table)
Real A[ns][ns];     // matrice A della tabella di Butcher (matrix A of the Butcher table)
Real K[ns][dim] = {{0}};  // matrice K per gli stadi (matrix K for the stages)
Real Z[dim];        // vettore per il risultato (vector for the result)

int main() {
    effe = F_lotkav;                                       // funzione del campo vettoriale (set vector field function)
    dati = dati_lotkav;                                     // dati iniziali del problema (initial problem data)
    sol = nullptr;                                          // funzione soluzione (solution function)
    butcher = DP87;                                         // impostazione del metodo Runge-Kutta (setting the Runge-Kutta method)
    butcher(ns, A[0], b, c, bc, 0);                         // impostazione della matrice di Butcher (setting the Butcher matrix)

    // Impostiamo il file su cui stampare (Setting up the output file)
    char n_file[21] = {0};  // nome del file (file name)
    cout << "inserire nome del file di stampa(max 20 caratteri):";  // richiede il nome del file (asks for file name)
    cin >> n_file;
    ofstream prt(n_file);  // apertura del file per la scrittura (open file for writing)
    prt.precision(14);      // impostazione della precisione di stampa (set print precision)
    cout.precision(14);     // impostazione della precisione di stampa (set print precision)

    // Input dei dati (Input data)
    dati(&t, &T, u); // richiede i dati iniziali (requests initial data)

    Real cpy[dim];  // variabile per copiare i dati (variable for copying data)
    copia(cpy, u, dim);  // copia i dati iniziali (copies initial data)

    unsigned long N;
    cout << "inserisci numero di passi N  = ";  // chiede il numero di passi (asks for the number of steps)
    cin >> N;
    h = (T - t) / Real(N);  // calcola il passo (calculates step size)

    Real toll = 1e-8;  // tolleranza per il passo adattivo (tolerance for adaptive step)
    Real r = 0.5;      // costante per il passo adattivo (constant for adaptive step)
    Real errmax = -1.0; // variabile per errore massimo (variable for maximum error)
    unsigned long n = 0;  // numero di passi (number of steps)

    for (int i = 0; i < ns; i++) {  // calcola la differenza b - bc per il metodo adattivo (calculates b - bc for adaptive method)
        bc[i] = b[i] - bc[i];  // calcolo della differenza b - bc (calculates the difference b - bc)
    }
    valutazioni_f = 0;  // reset valutazioni di f (reset evaluations of f)

    int yes = 0;  // flag per interrompere il ciclo (flag to break the loop)

    while (t < T && yes == 0) {  // ciclo principale (main loop)

        // Costruisco K per ogni stadio (Building K for each stage)
        for (int i = 0; i < ns; i++) {  // per ogni stadio (for each stage)
            Real aux[dim] = {0};  // array ausiliario (auxiliary array)
            step(u, h, A[i], K[0], i, dim, aux);  // calcola il passo (calculates step)
            effe(K[i], t + h * c[i], aux);  // aggiorna il campo vettoriale (update vector field)
        }

        // Calcolo l'errore e aggiorno h per il passo adattivo (Calculate error and update h for adaptive step)
        Real err[dim], Err;
        step(0, 1, bc, K[0], ns, dim, err);  // moltiplica i coefficienti bc con K (multiplies bc coefficients with K)
        Err = h * norm_2(err, dim);  // calcola l'errore (calculates error)
        hnew = pow(r * toll / Err, 1.0 / p) * h;  // aggiorna h in base all'errore (updates h based on error)

        valutazioni_f += 2;  // aggiorna il contatore delle valutazioni di f (updates the function evaluations counter)

        if (Err <= toll) {  // se l'errore  minore della tolleranza (if the error is smaller than the tolerance)
            Real Aux[dim] = {0};  // variabile ausiliaria (auxiliary variable)
            step(u, h, b, K[0], ns, dim, u);  // calcola la soluzione (calculate the solution)
            t += h;  // aggiorna il tempo (update time)
            h = min(hnew, T - t);  // aggiorna il passo (update the step size)
            n++;  // incrementa il numero di passi (increment step count)
            stampa(dim, t, u, sol, &prt);  // stampa i risultati (print the results)
            if (n == N) {  // se si  raggiunto il numero massimo di passi (if maximum number of steps is reached)
                yes = 1;  // interrompe il ciclo (break the loop)
            }
        } else {
            h = hnew;  // aggiorna h (update h)
        }
    }

    // Al termine del ciclo (At the end of the loop)
    for (int i = 0; i < dim; i++) {
        cpy[i] -= u[i];  // calcola la differenza finale (calculates the final difference)
    }

    cout << "n step= " << n << endl;  // numero di passi (number of steps)
    cout << "valutazioni f = " << valutazioni_f << endl;  // numero di valutazioni di f (number of function evaluations)
    // cout << "errore = " << norm_2(cpy, dim) << endl;  // errore finale (final error)
    return 0;
}
